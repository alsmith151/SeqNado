#!/usr/bin/env python3
from __future__ import annotations

import json
import os
import pathlib
import shutil
import subprocess
import sys
from importlib import resources
from typing import List, Optional, Any

import typer
import click
from loguru import logger
import pandera.pandas as pa
from pandas.api.types import is_bool_dtype, is_integer_dtype, is_float_dtype
import pandas as pd
from seqnado import Assay


# Optional: prettier tracebacks/console with rich
try:
    from rich import print  # noqa: F401
    from rich.traceback import install as _rich_tb_install

    _rich_tb_install(show_locals=False)
except Exception:
    pass

try:
    from rich.console import Console
    from rich.text import Text
    _RICH_CONSOLE = Console(force_terminal=True)
except Exception:  # rich not installed or no TTY
    _RICH_CONSOLE = None
    Text = None  # type: ignore

app = typer.Typer(
    add_completion=True,
    no_args_is_help=True,
    help="""
[bold]SeqNado CLI[/bold]

Initialize your environment, build configs, create design files, and run pipelines.
Use --help on any subcommand for details.
""",
)


# ------------------------------- Utilities ---------------------------------- #


def _pkg_files(pkg: str) -> pathlib.Path:
    """Return a concrete filesystem path to a package's resources directory."""
    return pathlib.Path(resources.files(pkg))  # PEP 451; zip-safe


def _read_json(path: pathlib.Path) -> dict:
    return json.loads(path.read_text())


def _write_json(path: pathlib.Path, data: dict) -> None:
    path.write_text(json.dumps(data, indent=4))
    try:
        os.chmod(path, 0o600)
    except Exception:
        pass


def _snakemake_available() -> bool:
    return shutil.which("snakemake") is not None


def _configure_logging(verbose: bool) -> None:
    logger.remove()
    logger.add(sys.stderr, level="DEBUG" if verbose else "INFO", colorize=True)


def _find_fastqs(hints: list[str]) -> list[pathlib.Path]:
    for loc in hints:
        files = sorted(pathlib.Path(loc).glob("*.fastq.gz"))
        if files:
            return files
    return []


# ---- Pandera schema helpers ------------------------------------------------- #


def _extract_candidate_defaults_from_schema(
    model: type[pa.DataFrameModel],
    assay: "Assay",
) -> dict[str, dict[str, Any]]:
    """
    Build a metadata dict for columns that we might add with defaults.

    - Ignores: r2, r1_control, r2_control, ip, control
    - Considers a column a candidate if it has a schema default OR is nullable

    Returns:
      {
        col_name: {
          "default": <value or None>,
          "nullable": <bool>,
          "description": <str or None>,
          "dtype": <pandas dtype-ish>,
          "categories": <list or None>   # for categoricals
        },
        ...
      }
    """
    schema = model.to_schema()
    out: dict[str, dict[str, Any]] = {}

    to_ignore = {"r2", "r1_control", "r2_control", "ip", "control"}
    for name, col in schema.columns.items():
        if name in to_ignore:
            continue

        # Pandera versions differ; use getattr defensively
        default = getattr(col, "default", None)
        nullable = bool(getattr(col, "nullable", False))
        if not (default is not None or nullable):
            continue

        description = getattr(col, "description", None)
        dtype = getattr(col, "dtype", None)

        # Extract categorical choices if possible
        categories = None
        try:
            if isinstance(dtype, pd.CategoricalDtype):
                categories = list(dtype.categories)
            elif getattr(dtype, "categories", None) is not None:
                categories = list(dtype.categories)
            elif str(dtype) == "category" and name == "assay":
                # Fallback for assay column if dtype is generic 'category'
                try:
                    categories = [a.value for a in type(assay)]
                except Exception:
                    categories = None
        except Exception:
            categories = None

        out[name] = {
            "default": default,
            "nullable": nullable,
            "description": description,
            "dtype": dtype,
            "categories": categories,
        }

    return out



def _style_name_with_rich(name: str, style: str = "bold cyan") -> str:
    """Return `name` styled with Rich as an ANSI string, or plain name on fallback."""
    if _RICH_CONSOLE is None or Text is None:
        return name
    with _RICH_CONSOLE.capture() as cap:
        _RICH_CONSOLE.print(Text(name, style=style), end="")
    return cap.get()


def _format_col_hint(name: str, meta: dict[str, Any]) -> str:
    """
    Build a hint like:
      <colored name>: <description> · choices=[...]
    """
    name_colored = _style_name_with_rich(name)

    parts: list[str] = []
    if meta.get("description"):
        parts.append(str(meta["description"]))

    cats = meta.get("categories")
    if cats:
        parts.append("choices=[" + ", ".join(map(str, cats)) + "]")

    return f"{name_colored}: " + " · ".join(parts) if parts else name_colored


def _coerce_value_to_dtype(
    value: str, dtype: Any, categories: Optional[List[Any]]
) -> Any:
    """
    Best-effort conversion from string input to the schema dtype.
    Keeps string if unsure. Enforces categorical choices if provided.
    """
    if categories is not None:
        if value == "" or value in categories:
            return value
        raise ValueError(f"Value must be one of: {', '.join(map(str, categories))}")

    if dtype is None:
        return value

    try:
        if is_bool_dtype(dtype):
            low = value.strip().lower()
            if low in {"true", "t", "yes", "y", "1"}:
                return True
            if low in {"false", "f", "no", "n", "0"}:
                return False
            raise ValueError("Enter a boolean (y/n, true/false, 1/0).")
        if is_integer_dtype(dtype):
            return int(value)
        if is_float_dtype(dtype):
            return float(value)
    except Exception as e:
        raise ValueError(str(e))

    return value


def _apply_interactive_defaults(
    df: pd.DataFrame,
    candidates: dict[str, dict[str, Any]],
    interactive: bool,
    accept_all_defaults: bool,
) -> pd.DataFrame:
    """
    For any candidate column not present in df, offer to add it.
    - If `accept_all_defaults` is True, auto-add only when a schema default exists.
    - If `interactive` is False, do nothing (safe for CI/batch).
    Prompts display Pandera description/dtype/nullable/choices.
    """
    missing = [c for c in candidates.keys() if c not in df.columns]
    if not missing:
        return df

    for col in missing:
        meta = candidates[col]
        default = meta.get("default", None)
        nullable = bool(meta.get("nullable", False))
        categories = meta.get("categories")
        dtype = meta.get("dtype")

        hint = _format_col_hint(col, meta)

        if accept_all_defaults and default is not None:
            logger.info(f"Adding '{col}' with schema default={default!r}  ({hint})")
            df[col] = default
            continue
        elif accept_all_defaults:
            # No schema default → skip in accept-all mode
            continue

        if not interactive:
            # Non-interactive and no mass-apply flag → skip
            continue

        # Interactive:
        if default is not None:
            add_col = typer.confirm(
                f"⚠️ Column '{col}' is missing.\n{hint}\nAdd it with default or enter 'n' to skip?",
                default=True,
            )
            if add_col:
                df[col] = default
        else:
            # No schema default; still nullable—offer to set a uniform value
            add_col = typer.confirm(
                f"⚠️ Optional column '{col}' is missing.\n{hint}\nAdd it with a uniform value or enter 'n' to skip?",
                default=False,
            )
            if add_col:
                while True:
                    raw = typer.prompt(
                        f"Enter a default value for '{col}' "
                        f"(press Enter for empty{f'; choices: {categories}' if categories else ''})",
                        default="",
                    )
                    # Empty string allowed for nullable columns
                    if raw == "" and nullable:
                        df[col] = raw
                        break
                    try:
                        coerced = _coerce_value_to_dtype(raw, dtype, categories)
                        df[col] = coerced
                        break
                    except ValueError as e:
                        typer.echo(f"[invalid] {e}")

    return df


# ------------------------------ Autocomplete -------------------------------- #


def _assay_names() -> list[str]:
    # Import lazily to keep CLI startup snappy
    from seqnado.inputs import Assay

    return list(Assay.all_assay_clean_names())


def assay_autocomplete(_: str) -> list[str]:
    return _assay_names()


def fastq_autocomplete(incomplete: str) -> list[str]:
    p = pathlib.Path(incomplete or ".")
    base = p.parent if p.name else p
    pattern = (p.name or "*") + ("*.fastq.gz" if not p.suffixes else "")
    return [str(x) for x in base.glob(pattern) if x.is_file()]


# --------------------------------- init ------------------------------------- #


@app.command(
    help="""
Initialize SeqNado user environment.

- Logs the current Conda environment if active (optional).
- Runs packaged Apptainer/Singularity init (if `apptainer` is on PATH).
- Ensures ~/.config/seqnado/genome_config.json exists (template or preset).
"""
)
def init(
    preset: bool = typer.Option(
        False, help="Use packaged preset genomes instead of the editable template."
    ),
    dry_run: bool = typer.Option(
        False, help="Show actions without writing files or running scripts."
    ),
    verbose: bool = typer.Option(
        False, "--verbose", "-v", help="Increase logging verbosity."
    ),
) -> None:
    _configure_logging(verbose)

    conda_env = os.environ.get("CONDA_DEFAULT_ENV")
    logger.info(f"Conda environment: {conda_env or 'none'}")

    # Apptainer/Singularity bootstrap
    if shutil.which("apptainer"):
        init_script = _pkg_files("seqnado").joinpath("init.sh")
        if init_script.exists():
            logger.info(f"Configuring Apptainer/Singularity via {init_script}")
            if not dry_run:
                try:
                    subprocess.run(["bash", str(init_script)], check=True)
                except subprocess.CalledProcessError as e:
                    logger.warning(f"Apptainer init script failed (continuing): {e}")
                except Exception as e:
                    logger.warning(f"Skipping Apptainer init due to error: {e}")
            else:
                logger.info(f"[dry-run] Would execute: bash {init_script}")
        else:
            logger.warning("Apptainer init script not found in package; skipping.")
    else:
        logger.info("Apptainer not found on PATH; skipping container setup.")

    # Genome config
    cfg_dir = pathlib.Path.home() / ".config" / "seqnado"
    cfg_dir.mkdir(parents=True, exist_ok=True)
    genome_config = cfg_dir / "genome_config.json"

    data_pkg = "seqnado.data"
    template_name = "preset_genomes.json" if preset else "genomes_template.json"
    template_path = _pkg_files(data_pkg).joinpath(template_name)

    if genome_config.exists():
        logger.info(f"Found genome config: {genome_config}")
        try:
            genome_data = _read_json(genome_config)
        except Exception as e:
            logger.warning(
                f"Could not parse existing genome config ({e}); leaving as-is."
            )
            genome_data = None

        if isinstance(genome_data, dict):
            has_placeholders = any(
                isinstance(v, str) and ("PATH" in v or v == "NA")
                for g in genome_data.values()
                if isinstance(g, dict)
                for v in g.values()
            )
            if has_placeholders:
                logger.info(
                    "Genome config contains placeholder paths—update before running pipelines."
                )
    else:
        if not template_path.exists():
            logger.error(
                "Packaged genome templates missing—installation may be corrupted."
            )
            raise typer.Exit(code=1)

        if dry_run:
            logger.info(
                f"[dry-run] Would create {genome_config} from {template_path.name}"
            )
        else:
            try:
                template = _read_json(template_path)
                _write_json(genome_config, template)
                logger.info(
                    "Created genome config "
                    + (
                        "from preset genomes."
                        if preset
                        else "template (please update paths)."
                    )
                )
            except Exception as e:
                logger.error(f"Failed to write genome config: {e}")
                raise typer.Exit(code=1)

    logger.success("Initialization complete.")


# -------------------------------- config ------------------------------------ #


@app.command(help="Build a workflow configuration YAML for the selected ASSAY.")
def config(
    assay: str = typer.Argument(
        ..., metavar="ASSAY", autocompletion=assay_autocomplete
    ),
    dont_make_directories: bool = typer.Option(
        False, help="Do not create the output project directory or fastq subdir."
    ),
    all_options: bool = typer.Option(
        False, help="Render all options (even if not used by the workflow)."
    ),
    output: Optional[pathlib.Path] = typer.Option(
        None, "-o", "--output", help="Explicit path for the rendered config file."
    ),
    verbose: bool = typer.Option(
        False, "--verbose", "-v", help="Increase logging verbosity."
    ),
) -> None:
    _configure_logging(verbose)

    from importlib.metadata import version as _pkg_version
    from datetime import date
    from seqnado.inputs import Assay
    from seqnado.config.user_input import build_workflow_config, render_config

    if assay not in Assay.all_assay_clean_names():
        allowed = ", ".join(Assay.all_assay_clean_names())
        logger.error(f"Unknown assay '{assay}'. Allowed: {allowed}")
        raise typer.Exit(code=2)

    seqnado_version = _pkg_version("seqnado")
    assay_obj = Assay.from_clean_name(assay)
    workflow_config = build_workflow_config(assay_obj, seqnado_version)
    if not workflow_config:
        logger.error("Failed to build workflow configuration.")
        raise typer.Exit(code=1)

    if dont_make_directories:
        config_output = output or pathlib.Path(f"config_{assay_obj.clean_name}.yaml")
        config_output.parent.mkdir(parents=True, exist_ok=True)
    else:
        dirname = f"{date.today().isoformat()}_{assay_obj.value}_{workflow_config.project.name}"
        outdir = pathlib.Path(dirname)
        (outdir / "fastqs").mkdir(parents=True, exist_ok=True)
        logger.info(f"Created output directory: {outdir / 'fastqs'}")
        config_output = output or (outdir / f"config_{assay_obj.clean_name}.yaml")

    tpl = _pkg_files("seqnado.data").joinpath("config_template.jinja")
    if not tpl.exists():
        logger.error("Packaged config template missing—installation may be corrupted.")
        raise typer.Exit(code=1)

    try:
        render_config(
            template=tpl,
            workflow_config=workflow_config,
            outfile=config_output,
            all_options=all_options,
        )
    except Exception as e:
        logger.error(f"Failed to render config: {e}")
        raise typer.Exit(code=1)

    logger.success(f"Wrote config → {config_output}")


# -------------------------------- design ------------------------------------ #


@app.command(help="Generate a SeqNado design CSV from FASTQ files for ASSAY.")
def design(
    assay: str = typer.Argument(
        ...,
        metavar="ASSAY",
        autocompletion=assay_autocomplete,
        show_choices=True,
        click_type=click.Choice(Assay.all_assay_clean_names(), case_sensitive=False),
        help="Assay type. Options: " + ", ".join(Assay.all_assay_clean_names()),
    ),
    files: List[pathlib.Path] = typer.Argument(
        None, metavar="[FASTQ ...]", autocompletion=fastq_autocomplete
    ),
    output: pathlib.Path = typer.Option(
        pathlib.Path("design.csv"), "-o", "--output", help="Output CSV filename."
    ),
    merge: bool = typer.Option(
        False, "--merge", help="Add a 'merge' column (for non-RNA assays)."
    ),
    auto_discover: bool = typer.Option(
        True,
        "--auto-discover/--no-auto-discover",
        help="Search common folders if none provided.",
    ),
    interactive: bool = typer.Option(
        True,
        "--interactive/--no-interactive",
        help="Interactively offer to add missing columns using schema defaults.",
    ),
    accept_all_defaults: bool = typer.Option(
        False,
        "--accept-all-defaults",
        help="Non-interactive: auto-add only columns that have a schema default.",
    ),
    verbose: bool = typer.Option(
        False, "--verbose", "-v", help="Increase logging verbosity."
    ),
) -> None:
    _configure_logging(verbose)

    from seqnado.inputs import FastqCollection, FastqCollectionForIP
    from seqnado.inputs.validation import DesignDataFrame

    if assay not in Assay.all_assay_clean_names():
        allowed = ", ".join(Assay.all_assay_clean_names())
        logger.error(f"Unknown assay '{assay}'. Allowed: {allowed}")
        raise typer.Exit(code=2)

    fastq_paths: list[pathlib.Path] = []
    for p in files or []:
        if p.suffixes[-2:] == [".fastq", ".gz"]:
            fastq_paths.append(p)

    if not fastq_paths and auto_discover:
        hints = [".", "fastqs", "fastq", "data", "data/fastqs"]
        logger.info(f"No FASTQs provided; searching {hints}")
        fastq_paths = _find_fastqs(hints)

    if not fastq_paths:
        logger.error("No FASTQ files provided or found.")
        typer.echo(
            "Tip: provide paths explicitly or place *.fastq.gz files in one of: ./, fastqs/, fastq/, data/, data/fastqs/",
            err=True,
        )
        raise typer.Exit(code=1)

    _assay = Assay.from_clean_name(assay)
    logger.info(f"Found {len(fastq_paths)} FASTQ files for assay '{_assay.value}'")

    if _assay in {"chip", "cat"}:
        design_obj = FastqCollectionForIP.from_fastq_files(
            assay=_assay, files=fastq_paths
        )
    else:
        design_obj = FastqCollection.from_fastq_files(assay=_assay, files=fastq_paths)

    df = design_obj.to_dataframe().sort_values("sample_id")

    # This adds default columns if missing, either interactively or via flags
    # based on the schema in DesignDataFrame
    # Gives the user a chance to fill in missing metadata e.g. consensus_group
    # Can ignore if running in batch/non-interactive mode
    schema_candidates = _extract_candidate_defaults_from_schema(DesignDataFrame, _assay)
    df = _apply_interactive_defaults(
        df,
        schema_candidates,
        interactive=interactive,
        accept_all_defaults=accept_all_defaults,
    )

    output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output, index=False)
    logger.success(f"Design file saved → {output}")


# -------------------------------- pipeline ---------------------------------- #
# Allow pass-through of *unknown* options to Snakemake via ctx.args
@app.command(
    help="Run the data processing pipeline for ASSAY (Snakemake under the hood).",
    context_settings={"allow_extra_args": True, "ignore_unknown_options": True},
)
def pipeline(
    ctx: typer.Context,
    assay: Optional[str] = typer.Argument(
        None, metavar="ASSAY", autocompletion=assay_autocomplete
    ),
    config_file: Optional[pathlib.Path] = typer.Option(
        None,
        "--configfile",
        help="Path to a SeqNado config YAML (default: config_ASSAY.yaml).",
    ),
    show_version: bool = typer.Option(
        False, "--version", help="Print SeqNado version and exit."
    ),
    preset: str = typer.Option(
        "lc",
        "--preset",
        help="Snakemake job profile preset.",
        show_choices=True,
        case_sensitive=False,
    ),
    clean_symlinks: bool = typer.Option(
        False, help="Remove symlinks created by previous runs."
    ),
    scale_resources: float = typer.Option(
        1.0, "-s", "--scale-resources", help="Scale memory/time (env: SCALE_RESOURCES)."
    ),
    verbose: bool = typer.Option(
        False, "--verbose", "-v", help="Increase logging verbosity."
    ),
    queue: Optional[str] = typer.Option(
        None, "-q", "--queue", help="Slurm queue/partition for the `ss` preset."
    ),
    print_cmd: bool = typer.Option(
        False, "--print-cmd", help="Print the Snakemake command before running it."
    ),
) -> None:
    _configure_logging(verbose)

    if show_version:
        from importlib.metadata import version as _pkg_version

        logger.info(f"SeqNado version {_pkg_version('seqnado')}")
        raise typer.Exit(code=0)

    if not _snakemake_available():
        logger.error(
            "`snakemake` not found on PATH. Install/activate the environment that provides it."
        )
        raise typer.Exit(code=127)

    if not assay:
        logger.error(
            "No assay provided. Use `seqnado pipeline ASSAY` or `--version` to print the version."
        )
        raise typer.Exit(code=2)

    # Lazy import to keep --help snappy
    from seqnado.helpers import extract_cores_from_options

    # Typer/Click hands us all unknown flags/args via ctx.args; pass them to Snakemake
    extra_args = list(ctx.args)
    cleaned_opts, cores = extract_cores_from_options(extra_args)

    # Scale resources via env var
    os.environ["SCALE_RESOURCES"] = str(scale_resources)

    # Clean old symlinks if requested
    if clean_symlinks:
        target = pathlib.Path("seqnado_output/fastqs")
        logger.info(f"Cleaning symlinks in {target} ...")
        for link in target.glob("*"):
            if link.is_symlink():
                link.unlink(missing_ok=True)

    # Build Snakemake command
    pkg_root = _pkg_files("seqnado")
    snakefile = pkg_root / "workflow" / "Snakefile"
    if not snakefile.exists():
        logger.error(f"Snakefile for assay '{assay}' not found: {snakefile}")
        raise typer.Exit(code=1)

    cmd: list[str] = [
        "snakemake",
        "-c",
        str(cores),
        "--snakefile",
        str(snakefile),
        "--show-failed-logs",
        "--configfile",
        str(config_file or pathlib.Path(f"config_{assay}.yaml")),
        *cleaned_opts,
    ]

    if preset.lower() == "ss":
        profile = (
            pkg_root / "workflow" / "envs" / "profiles" / "profile_slurm_singularity"
        )
        cmd += [
            "--profile",
            str(profile),
            "--default-resources",
            f"slurm_partition={queue or 'short'}",
        ]
    elif preset.lower() == "ls":
        profile = pkg_root / "workflow" / "envs" / "profiles" / "profile_singularity"
        cmd += ["--profile", str(profile)]
    # lc -> no extra flags (handled by environment/snakefile/profile conventions)

    # Print the logo (best-effort)
    logo_path = pkg_root / "data" / "logo.txt"
    if logo_path.exists():
        try:
            print(logo_path.read_text())
        except Exception:
            pass

    # Resolve symlinked PWD for Singularity bind compatibility
    cwd = str(pathlib.Path(".").resolve())
    os.chdir(cwd)
    os.environ["PWD"] = cwd

    if print_cmd:
        logger.info("Snakemake command:\n$ " + " ".join(map(str, cmd)))

    completed = subprocess.run(cmd, cwd=cwd)
    raise typer.Exit(code=completed.returncode)


# -------------------------------- Entrypoint -------------------------------- #

if __name__ == "__main__":
    app()
