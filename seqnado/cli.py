#!/usr/bin/env python3
from __future__ import annotations

import json
import os
import shlex
import shutil
import subprocess
import sys
import tempfile
from datetime import date, datetime
from importlib import resources
from pathlib import Path
from typing import Any, List, Optional, TYPE_CHECKING

import typer
from loguru import logger
import click

# Optional: prettier tracebacks/console with rich if available
try:
    from rich.traceback import install as _rich_tb_install

    _rich_tb_install(show_locals=False)
except Exception:
    pass

# Rich console capture may or may not be available — keep an opt-in safe reference
try:
    from rich.console import Console  # type: ignore
    from rich.text import Text  # type: ignore

    _RICH_CONSOLE = Console(force_terminal=True)
except Exception:
    _RICH_CONSOLE = None
    Text = None  # type: ignore

if TYPE_CHECKING:
    # Available to type-checkers only; does not execute at runtime.
    import pandas as pd

from seqnado import Assay

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


def _pkg_traversable(pkg: str):
    """
    Return the importlib.resources Traversable for a package.
    Use resources.as_file(...) when you need a filesystem Path for a particular file.
    """
    return resources.files(pkg)


def _read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_json(path: Path, data: dict) -> None:
    """
    Atomically write JSON to `path` using a temporary file + os.replace.
    Attempt to set secure permissions but ignore if not supported.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = Path(tempfile.mktemp(dir=str(path.parent)))
    try:
        tmp.write_text(json.dumps(data, indent=4), encoding="utf-8")
        os.replace(str(tmp), str(path))
        try:
            os.chmod(path, 0o600)
        except Exception:
            # best-effort — ignore permission errors
            logger.debug("Unable to set file permissions on %s", path)
    finally:
        if tmp.exists():
            try:
                tmp.unlink(missing_ok=True)
            except Exception:
                # best-effort cleanup
                pass


def _snakemake_available() -> bool:
    return shutil.which("snakemake") is not None


def _configure_logging(verbose: bool) -> None:
    logger.remove()
    logger.add(sys.stderr, level="DEBUG" if verbose else "INFO", colorize=True)


def _find_fastqs(hints: List[str]) -> List[Path]:
    """
    Search the provided hint directories for *.fastq.gz files.
    Skip hints that don't exist and return a sorted list.
    """
    out: List[Path] = []
    for loc in hints:
        p = Path(loc)
        if not p.exists():
            continue
        out.extend(sorted(p.glob("*.fastq.gz")))
    return out


def _style_name_with_rich(name: str, style: str = "bold cyan") -> str:
    """Return `name` styled with Rich as an ANSI string, or plain name on fallback."""
    if _RICH_CONSOLE is None or Text is None:
        return name
    with _RICH_CONSOLE.capture() as cap:
        _RICH_CONSOLE.print(Text(name, style=style), end="")
    return cap.get()


# ---- Pandera schema helpers ------------------------------------------------- #
# (these functions are defensive and robust against Pandera API differences)


def _extract_candidate_defaults_from_schema(
    model: type, assay: Any
) -> dict[str, dict[str, Any]]:
    """
    Build a metadata dict for columns that we might add with defaults.

    - Ignores: r2, r1_control, r2_control, ip, control
    - Considers a column a candidate if it has a schema default OR is nullable
    """
    try:
        schema = model.to_schema()  # Pandera DataFrameModel -> Schema
    except Exception as e:
        logger.debug("Could not convert model to schema: %s", e)
        return {}

    out: dict[str, dict[str, Any]] = {}
    to_ignore = {"r2", "r1_control", "r2_control", "ip", "control"}

    for name, col in schema.columns.items():
        if name in to_ignore:
            continue

        default = getattr(col, "default", None)
        nullable = bool(getattr(col, "nullable", False))
        if not (default is not None or nullable):
            continue

        description = getattr(col, "description", None)
        dtype = getattr(col, "dtype", None)

        categories = None
        try:
            import pandas as pd

            if isinstance(dtype, pd.CategoricalDtype):
                categories = list(dtype.categories)
            elif getattr(dtype, "categories", None) is not None:
                categories = list(dtype.categories)
            elif str(dtype) == "category" and name == "assay":
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


def _format_col_hint(name: str, meta: dict[str, Any]) -> str:
    """
    Build a hint like: <colored name>: <description> · choices=[...]
    """
    name_colored = _style_name_with_rich(name)

    parts: List[str] = []
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
    import pandas as pd
    from pandas.api import types as _pdt

    if categories is not None:
        # Allow empty/nullable
        if value == "" or value in categories:
            return value
        raise ValueError(f"Value must be one of: {', '.join(map(str, categories))}")

    if dtype is None:
        return value

    try:
        pd_dtype = pd.api.types.pandas_dtype(dtype)
    except Exception:
        pd_dtype = dtype

    # boolean handling
    if _pdt.is_bool_dtype(pd_dtype):
        low = value.strip().lower()
        if low in {"true", "t", "yes", "y", "1"}:
            return True
        if low in {"false", "f", "no", "n", "0"}:
            return False
        raise ValueError("Enter a boolean (y/n, true/false, 1/0).")

    # integer
    if _pdt.is_integer_dtype(pd_dtype):
        try:
            return int(value)
        except Exception as e:
            raise ValueError(str(e))

    # float
    if _pdt.is_float_dtype(pd_dtype):
        try:
            return float(value)
        except Exception as e:
            raise ValueError(str(e))

    # fallback to string
    return value


def _apply_interactive_defaults(
    df_in: pd.DataFrame,
    candidates: dict[str, dict[str, Any]],
    interactive: bool,
    accept_all_defaults: bool,
) -> pd.DataFrame:
    """
    For any candidate column not present in df_in, optionally add it.

    - Returns a new DataFrame (does not mutate the input).
    - If `accept_all_defaults` is True, auto-add only when a schema default exists.
    - If `interactive` is False, do nothing (safe for CI/batch).
    """
    import pandas as pd

    df = df_in.copy()
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
            try:
                df[col] = pd.Series([default] * len(df), index=df.index)
            except Exception:
                df[col] = default
            continue
        elif accept_all_defaults:
            continue

        if not interactive:
            continue

        # Interactive path
        if default is not None:
            add_col = typer.confirm(
                f"⚠️ Column '{col}' is missing.\n{hint}\nAdd it with default or enter 'n' to skip?",
                default=True,
            )
            if add_col:
                try:
                    df[col] = pd.Series([default] * len(df), index=df.index)
                except Exception:
                    df[col] = default
        else:
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
                    if raw == "" and nullable:
                        df[col] = pd.Series([pd.NA] * len(df), index=df.index)
                        break
                    try:
                        coerced = _coerce_value_to_dtype(raw, dtype, categories)
                        try:
                            df[col] = pd.Series([coerced] * len(df), index=df.index)
                        except Exception:
                            df[col] = coerced
                        break
                    except ValueError as e:
                        typer.echo(f"[invalid] {e}")

    return df


# ------------------------------ Autocomplete -------------------------------- #


def _assay_names() -> List[str]:
    # Import lazily to keep CLI startup snappy
    from seqnado.inputs import Assay  # local import

    return list(Assay.all_assay_clean_names())


def assay_autocomplete(_: str) -> List[str]:
    return _assay_names()


def fastq_autocomplete(incomplete: str) -> List[str]:
    # Provide some file suggestions without expensive globbing
    p = Path(incomplete or ".")
    base = p.parent if p.name else p
    pattern = p.name or "*"
    try:
        candidates = sorted(base.glob(pattern + "*.fastq.gz"))
    except Exception:
        candidates = []
    return [str(x) for x in candidates if x.is_file()][:40]


def _get_profile_name(fn: Path) -> Optional[str]:
    name = fn.name
    if not name.startswith("profile_"):
        return None

    profile_parts = name.split("_")
    if len(profile_parts) < 2:
        return None
    initials = "".join(part[0] for part in profile_parts[1:] if part)
    return initials


def _preset_profiles() -> List[str]:
    profiles = [
        f.name
        for f in _pkg_traversable("seqnado.workflow.envs.profiles").iterdir()
        if f.is_dir()
    ]

    # At the moment these referenced by two unique letters based on the name after profile_ e.g. 'slurm_sigularity' -> 'ss'
    profiles = {
        _get_profile_name(Path(p)): p for p in profiles if p.startswith("profile_")
    }
    return profiles


def _profile_autocomplete() -> List[str]:
    return list(_preset_profiles().keys())


# --------------------------------- init ------------------------------------- #


@app.command(
    help="""
Initialize SeqNado user environment.

- Logs the current Conda environment if active (optional).
- Runs packaged Apptainer/Singularity init (if `apptainer` on PATH).
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
        init_script_trav = _pkg_traversable("seqnado").joinpath("init.sh")
        try:
            with resources.as_file(init_script_trav) as init_script:
                if init_script.exists():
                    logger.info(f"Configuring Apptainer/Singularity via {init_script}")
                    if not dry_run:
                        try:
                            subprocess.run(["bash", str(init_script)], check=True)
                        except subprocess.CalledProcessError as e:
                            logger.warning(
                                f"Apptainer init script failed (continuing): {e}"
                            )
                        except Exception as e:
                            logger.warning(f"Skipping Apptainer init due to error: {e}")
                    else:
                        logger.info(f"[dry-run] Would execute: bash {init_script}")
                else:
                    logger.warning(
                        "Apptainer init script not found in package; skipping."
                    )
        except Exception as e:
            logger.warning("Could not access package init script: %s", e)
    else:
        logger.info("Apptainer not found on PATH; skipping container setup.")

    # Genome config
    cfg_dir = Path.home().joinpath(".config", "seqnado")
    cfg_dir.mkdir(parents=True, exist_ok=True)
    genome_config = cfg_dir.joinpath("genome_config.json")

    data_pkg = "seqnado.data"
    template_genomes = "preset_genomes.json"
    template_config = "genomes_template.json"
    template_name = template_genomes if preset else template_config
    template_trav = _pkg_traversable(data_pkg).joinpath(template_name)
    template_config_trav = _pkg_traversable(data_pkg).joinpath(template_config)

    
    # move to ~/.config/snakemake.
    profile_target_dir = Path.home().joinpath(".config", "snakemake")
    profile_target_dir.mkdir(parents=True, exist_ok=True)
    profiles = _preset_profiles()
    for profile in profiles.values():
        profile_src_trav = _pkg_traversable("seqnado.workflow.envs.profiles").joinpath(
            profile
        )
        profile_dest = profile_target_dir.joinpath(profile)
        if profile_dest.exists():
            logger.info(f"Snakemake profile already exists: {profile_dest}")
            continue
        try:
            with resources.as_file(profile_src_trav) as profile_src:
                if dry_run:
                    logger.info(
                        f"[dry-run] Would copy Snakemake profile {profile_src.name} to {profile_dest}"
                    )
                else:
                    shutil.copytree(profile_src, profile_dest)
                    logger.info(f"Copied Snakemake profile to {profile_dest}")
        except Exception as e:
            logger.error(f"Failed to copy Snakemake profile {profile}: {e}")

    if genome_config.exists():
        logger.info(f"Found genome config: {genome_config}")
        try:
            genome_data = _read_json(genome_config)
        except Exception as e:
            logger.warning(
                f"Could not parse existing genome config ({e}); leaving as-is."
            )
            genome_data = None

        try:
            with resources.as_file(template_config_trav) as template_path:
                template_data = _read_json(Path(template_path))
            if genome_data == template_data:
                logger.warning(
                    "Genome config matches the template exactly. Please update paths as needed."
                )
            else:
                logger.info("Genome config appears to be customized; leaving as-is.")
        except Exception as e:
            logger.debug("Could not read packaged template to compare: %s", e)
    else:
        try:
            with resources.as_file(template_trav) as template_path:
                if not Path(template_path).exists():
                    logger.error(
                        "Packaged genome templates missing—installation may be corrupted."
                    )
                    raise typer.Exit(code=1)
                if dry_run:
                    logger.info(
                        f"[dry-run] Would create {genome_config} from {template_path.name}"
                    )
                else:
                    template = _read_json(Path(template_path))
                    _write_json(genome_config, template)
                    logger.info(
                        "Created genome config "
                        + (
                            "from preset genomes."
                            if preset
                            else "template (please update paths)."
                        )
                    )
        except typer.Exit:
            raise
        except Exception as e:
            logger.error(f"Failed to write genome config: {e}")
            raise typer.Exit(code=1)

    logger.success("Initialization complete.")


# -------------------------------- utils ----------------------------------- #


@app.command()
def genomes(
    subcommand: str = typer.Argument(..., help="Subcommand: list | edit | build"),
    assay: str = typer.Argument(
        "atac",
        metavar="ASSAY",
        autocompletion=assay_autocomplete,
        show_choices=True,
        help="Assay type. Options: " + ", ".join(_assay_names()),
    ),
    fasta: Optional[Path] = typer.Option(
        None, "--fasta", "-f", help="Input FASTA (required for build)"
    ),
    name: Optional[str] = typer.Option(
        None, "--name", "-n", help="Genome name (prefix) for built genome"
    ),
    outdir: Path = typer.Option(
        Path.cwd() / "genome_build", "--outdir", "-o", help="Output directory for build"
    ),
) -> None:
    """
    Subcommands for genomes:
      - list  : show packaged and user genome presets
      - edit  : open user genome config in $EDITOR (creates if missing)
      - build : build a genome from a fasta using a Snakemake pipeline
    """
    # Import locally for snappy startup
    from seqnado.config import load_genome_configs  # local import

    sub = subcommand.lower().strip()
    if sub == "list":
        try:
            cfg = load_genome_configs(assay=assay)
        except Exception as e:
            logger.error(f"Failed to load genome config: {e}")
            raise typer.Exit(code=1)

        if not cfg:
            logger.warning("No genome config found.")
            raise typer.Exit(code=0)

        for name, details in cfg.items():
            typer.echo(f"[bold]{name}[/bold]")
            try:
                items = details.dict()
            except Exception:
                items = dict(details)
            for k, v in items.items():
                typer.echo(f"  {k}: {v or '[not set]'}")
            typer.echo("")
        raise typer.Exit(code=0)

    elif sub == "edit":
        env = os.environ.get("SEQNADO_GENOME_CONFIG")
        cfg_path = (
            Path(env)
            if env
            else Path.home().joinpath(".config", "seqnado", "genome_config.json")
        )
        if not cfg_path.exists():
            logger.error(
                f"Genome config not found: {cfg_path} (try `seqnado init` first)"
            )
            raise typer.Exit(code=1)

        editor = (
            os.environ.get("EDITOR")
            or os.environ.get("VISUAL")
            or ("notepad" if sys.platform == "win32" else "nano")
        )
        editor_cmd = shlex.split(editor)
        if shutil.which(editor_cmd[0]) is None:
            logger.error(
                f"Editor '{editor_cmd[0]}' not found. Please set $EDITOR to your preferred editor."
            )
            raise typer.Exit(code=4)

        try:
            subprocess.check_call(editor_cmd + [str(cfg_path)])
        except subprocess.CalledProcessError as e:
            logger.error(f"Editor returned non-zero exit status: {e}")
            raise typer.Exit(code=3)
        except FileNotFoundError:
            logger.error(
                f"Editor '{editor_cmd[0]}' not found. Please set $EDITOR to your preferred editor."
            )
            raise typer.Exit(code=4)
        except Exception as e:
            logger.error(f"Failed to launch editor: {e}")
            raise typer.Exit(code=5)

        typer.echo(f"Edited genome config: [italic]{cfg_path}[/italic]")
        raise typer.Exit(code=0)

    elif sub == "build":
        if not fasta:
            logger.error("The --fasta option is required for 'build'.")
            raise typer.Exit(code=2)

        fasta = Path(fasta).expanduser().resolve()
        if not fasta.exists():
            logger.error(f"FASTA not found: {fasta}")
            raise typer.Exit(code=3)

        outdir = Path(outdir).expanduser().resolve()
        outdir.mkdir(parents=True, exist_ok=True)

        raise NotImplementedError("Genome build not yet implemented.")

    else:
        logger.error("Unknown subcommand for genomes: %s", sub)
        raise typer.Exit(code=2)


# -------------------------------- config ------------------------------------ #


@app.command(help="Build a workflow configuration YAML for the selected ASSAY.")
def config(
    assay: str = typer.Argument(
        ...,
        metavar="ASSAY",
        autocompletion=assay_autocomplete,
        show_choices=True,
        help=", ".join(_assay_names()),
    ),
    make_dirs: bool = typer.Option(
        True,
        "--make-dirs/--no-make-dirs",
        help="Create/don't create the output project directory or fastq subdir.",
    ),
    render_options: bool = typer.Option(
        False, help="Render all options (even if not used by the workflow)."
    ),
    output: Optional[Path] = typer.Option(
        None, "-o", "--output", help="Explicit path for the rendered config file."
    ),
    verbose: bool = typer.Option(
        False, "--verbose", "-v", help="Increase logging verbosity."
    ),
    interactive: bool = typer.Option(
        True,
        "--interactive/--no-interactive",
        help="Interactively prompt for config values",
    ),
) -> None:
    _configure_logging(verbose)

    # Local imports to keep help fast
    from importlib.metadata import version as _pkg_version
    from seqnado.inputs import Assay
    from seqnado.config import (
        build_workflow_config,
        render_config,
        build_default_workflow_config,
    )

    if assay not in Assay.all_assay_clean_names():
        allowed = ", ".join(Assay.all_assay_clean_names())
        logger.error(f"Unknown assay '{assay}'. Allowed: {allowed}")
        raise typer.Exit(code=2)

    seqnado_version = _pkg_version("seqnado")
    assay_obj = Assay.from_clean_name(assay)

    if not interactive:
        logger.info("Running in non-interactive mode; using defaults where possible.")
        workflow_config = build_default_workflow_config(assay_obj)
        render_options = True
    else:
        workflow_config = build_workflow_config(assay_obj, seqnado_version)
        if not workflow_config:
            logger.error("Failed to build workflow configuration.")
            raise typer.Exit(code=1)

    if not make_dirs:
        config_output = output or Path(f"config_{assay_obj.clean_name}.yaml")
        config_output.parent.mkdir(parents=True, exist_ok=True)
    else:
        dirname = f"{date.today().isoformat()}_{assay_obj.clean_name}_{workflow_config.project.name}"
        outdir = Path(dirname)
        (outdir / "fastqs").mkdir(parents=True, exist_ok=True)
        logger.info(f"Created output directory: {outdir / 'fastqs'}")
        config_output = output or (outdir / f"config_{assay_obj.clean_name}.yaml")

    tpl_trav = _pkg_traversable("seqnado.data").joinpath("config_template.jinja")
    try:
        with resources.as_file(tpl_trav) as tpl_path:
            if not Path(tpl_path).exists():
                logger.error(
                    "Packaged config template missing—installation may be corrupted."
                )
                raise typer.Exit(code=1)
            try:
                render_config(
                    template=Path(tpl_path),
                    workflow_config=workflow_config,
                    outfile=config_output,
                    all_options=render_options,
                )
            except Exception as e:
                logger.error(f"Failed to render config: {e}")
                raise typer.Exit(code=1)
    except typer.Exit:
        raise
    except Exception as e:
        logger.error("Failed to locate or access packaged template: %s", e)
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
        help="Assay type. Options: " + ", ".join(Assay.all_assay_clean_names()),
    ),
    files: List[Path] = typer.Argument(
        None, metavar="[FASTQ ...]", autocompletion=fastq_autocomplete
    ),
    output: Path = typer.Option(
        Path("metadata.csv"), "-o", "--output", help="Output CSV filename."
    ),
    group_by: bool = typer.Option(
        False, "--group-by", help="Group samples by a regular expression or a column."
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

    # Local imports
    from seqnado.inputs import FastqCollection, FastqCollectionForIP
    from seqnado.inputs.validation import DesignDataFrame
    from seqnado.inputs import Assay as AssayEnum
    import pandas as pd

    if assay not in AssayEnum.all_assay_clean_names():
        allowed = ", ".join(AssayEnum.all_assay_clean_names())
        logger.error(f"Unknown assay '{assay}'. Allowed: {allowed}")
        raise typer.Exit(code=2)

    fastq_paths: List[Path] = []
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

    _assay = AssayEnum.from_clean_name(assay)
    logger.info(f"Found {len(fastq_paths)} FASTQ files for assay '{_assay.value}'")

    if _assay in {AssayEnum.CHIP, AssayEnum.CAT}:
        design_obj = FastqCollectionForIP.from_fastq_files(
            assay=_assay, files=fastq_paths
        )
    else:
        design_obj = FastqCollection.from_fastq_files(assay=_assay, files=fastq_paths)

    df = design_obj.to_dataframe().sort_values("sample_id")

    schema_candidates = _extract_candidate_defaults_from_schema(DesignDataFrame, _assay)
    df = _apply_interactive_defaults(
        df,
        schema_candidates,
        interactive=interactive,
        accept_all_defaults=accept_all_defaults,
    )

    if group_by:
        if isinstance(group_by, str) and group_by in df.columns:
            df["consensus_group"] = df[group_by].astype(str)
            logger.info(
                f"Grouped samples by column '{group_by}' into 'consensus_group'."
            )
        else:
            # Treat group_by as a regex pattern to extract from sample_id
            try:
                if _assay in [AssayEnum.CHIP, AssayEnum.CAT]:
                    samples = df["sample_id"] + df["ip"]
                else:
                    samples = df["sample_id"]

                df["consensus_group"] = samples.str.extract(group_by, expand=False)
                if df["consensus_group"].isnull().all():
                    raise ValueError(
                        f"No matches found with the provided regex '{group_by}'"
                    )

                df["consensus_group"] = df["consensus_group"].fillna("unknown")
                logger.info(
                    f"Grouped samples by regex '{group_by}' into 'consensus_group'."
                )
            except Exception as e:
                logger.error(f"Failed to group by '{group_by}': {e}")
                raise typer.Exit(code=3)

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
    config_file: Optional[Path] = typer.Option(
        None,
        "--configfile",
        help="Path to a SeqNado config YAML (default: config_<ASSAY>.yaml).",
    ),
    show_version: bool = typer.Option(
        False, "--version", help="Print SeqNado version and exit."
    ),
    preset: str = typer.Option(
        "le",
        "--preset",
        click_type=click.Choice(_profile_autocomplete(), case_sensitive=False),
        help="Snakemake job profile preset.",
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
        'short', "-q", "--queue", help="Slurm queue/partition for the `ss` preset."
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

    # Local import for helper
    from seqnado.helpers import extract_cores_from_options

    extra_args = list(ctx.args)
    cleaned_opts, cores = extract_cores_from_options(extra_args)

    os.environ["SCALE_RESOURCES"] = str(scale_resources)

    if clean_symlinks:
        target = Path("seqnado_output/fastqs")
        logger.info(f"Cleaning symlinks in {target} ...")
        for link in target.glob("*"):
            if link.is_symlink():
                link.unlink(missing_ok=True)

    pkg_root_trav = _pkg_traversable("seqnado")
    snake_trav = pkg_root_trav.joinpath("workflow").joinpath("Snakefile")

    if not config_file:
        config_file = Path(f"config_{assay}.yaml")
        if not config_file.exists():
            logger.error(
                f"No config file provided and default not found: {config_file}"
            )
            raise typer.Exit(code=1)

    # Use resources.as_file to ensure the packaged Snakefile is available as a path
    try:
        with resources.as_file(snake_trav) as snakefile_path:
            if not Path(snakefile_path).exists():
                logger.error(
                    f"Snakefile for assay '{assay}' not found: {snakefile_path}"
                )
                raise typer.Exit(code=1)

            cmd: List[str] = [
                "snakemake",
                "-c",
                str(cores),
                "--snakefile",
                str(snakefile_path),
                "--show-failed-logs",
                "--configfile",
                str(config_file),
                *cleaned_opts,
            ]

            if preset:
                profiles = _preset_profiles()
                profile_dir = profiles.get(preset.lower())
                cmd += ["--profile", profile_dir]
                logger.info(
                    f"Using Snakemake profile preset '{preset}' -> {profile_dir}"
                )

            if queue and preset.startswith("s"):
                cmd += ['--default-resources', f'slurm_partition={queue}']

            logo_trav = pkg_root_trav.joinpath("data").joinpath("logo.txt")
            try:
                with resources.as_file(logo_trav) as lp:
                    try:
                        print(Path(lp).read_text())
                    except Exception:
                        pass
            except Exception:
                pass

            cwd = str(Path(".").resolve())
            os.chdir(cwd)
            os.environ["PWD"] = cwd

            if print_cmd:
                logger.info("Snakemake command:\n$ " + " ".join(map(str, cmd)))

            completed = subprocess.run(cmd, cwd=cwd)
            raise typer.Exit(code=completed.returncode)
    except typer.Exit:
        raise
    except Exception as e:
        logger.exception("Failed to run snakemake: %s", e)
        raise typer.Exit(code=1)


# -------------------------------- Entrypoint --------------------------------
if __name__ == "__main__":
    app()
