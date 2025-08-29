import shlex
from typing import Optional, Literal, Any, Annotated
from enum import Enum
from pathlib import Path

from pydantic import BaseModel, Field, field_validator, model_validator
from pydantic.functional_serializers import PlainSerializer

from seqnado import Assay
from .mixins import PathValidatorMixin


# =============================================================================
# Core Configuration Classes
# =============================================================================

class CommandLineArguments(BaseModel):
    """Base class for CLI options with validation and filtering capabilities."""
    
    value: str = Field(default="", description="CLI options string")
    exclude: set[str] = Field(
        default_factory=set, description="Options to exclude from the final command"
    )

    @field_validator("value", mode="before")
    @classmethod
    def validate_options_syntax(cls, v: str) -> str:
        if not isinstance(v, str):
            raise ValueError("Value must be a string")
        try:
            shlex.split(v)
        except ValueError as e:
            raise ValueError(f"Invalid CLI option string: {e}")
        if any(c in v for c in ["\n", "\r", "\x00"]):
            raise ValueError("Value string contains unsafe characters")
        return v

    @model_validator(mode="after")
    def ensure_leading_space(self) -> "CommandLineArguments":
        if self.value and not self.value.startswith(" "):
            self.value = " " + self.value
        return self

    def _filter_excluded_options(self, options_str: str) -> str:
        if not self.exclude:
            return options_str
        try:
            tokens = shlex.split(options_str.strip())
            filtered_tokens = []
            i = 0
            while i < len(tokens):
                token = tokens[i]
                excluded = False
                for exclude_pattern in self.exclude:
                    if token == exclude_pattern or token.startswith(
                        exclude_pattern + "="
                    ):
                        excluded = True
                        if (
                            token == exclude_pattern
                            and i + 1 < len(tokens)
                            and not tokens[i + 1].startswith("-")
                        ):
                            i += 1
                        break
                if not excluded:
                    filtered_tokens.append(token)
                i += 1
            return " " + shlex.join(filtered_tokens) if filtered_tokens else ""
        except ValueError:
            return options_str

    @property
    def option_string_filtered(self) -> str:
        return self._filter_excluded_options(self.value)

    @property
    def option_string_raw(self) -> str:
        return self.value
    
    def __str__(self) -> str:
        return self._filter_excluded_options(self.value).strip()


def serialize_options_to_string(options: CommandLineArguments) -> str:
    """Custom serializer that returns only the filtered option string."""
    return options.option_string_filtered.strip()


# Create the annotated Options type with custom serialization
CommandLineArgumentsType = Annotated[
    CommandLineArguments,
    PlainSerializer(serialize_options_to_string, return_type=str)
]


class ToolConfig(BaseModel):
    """Base configuration for a tool with threads and CLI options."""
    
    threads: int = Field(default=1, ge=1, description="Number of threads to use")
    command_line_arguments: CommandLineArgumentsType = Field(default_factory=CommandLineArguments, description="CLI options")

    @field_validator("command_line_arguments", mode="before")
    @classmethod
    def validate_command_line_arguments(cls, v):
        if isinstance(v, CommandLineArguments):
            return v
        if isinstance(v, str):
            return CommandLineArguments(value=v)
        if isinstance(v, dict):
            return CommandLineArguments(**v)
        raise TypeError("command_line_arguments must be a string, dict, or CommandLineArguments instance")


# =============================================================================
# Enums
# =============================================================================

class CutadaptMode(Enum):
    DEFAULT = "default"
    CRISPR = "crispr"


class MacsMode(Enum):
    ATAC = "atac"
    CHIP_BROAD = "chip_broad"
    CHIP_TF = "chip_tf"
    GENERIC = "generic"


# =============================================================================
# Simple Tool Configurations (single ToolConfig)
# =============================================================================

class Cutadapt(ToolConfig):
    """Cutadapt adapter trimming tool configuration."""
    
    mode: CutadaptMode = Field(default=CutadaptMode.CRISPR)

    @model_validator(mode="before")
    @classmethod
    def set_defaults(cls, data: dict) -> dict:
        data = dict(data)
        data.setdefault("threads", 4)
        data.setdefault(
            "options",
            CommandLineArguments(value="-g 'ACACCG' --cut 0 -l 20 -m 20 -M 20 --discard-untrimmed"),
        )
        return data


class LanceotronMCC(ToolConfig):
    """Lanceotron MCC peak calling tool configuration."""
    
    @model_validator(mode="before")
    @classmethod
    def set_defaults(cls, data: dict) -> dict:
        data = dict(data)
        data.setdefault("threads", 8)
        data.setdefault(
            "options",
            CommandLineArguments(
                value="--max-peak-size 2000 --no-mcc-filter-enriched-regions --no-mcc-filter-background"
            ),
        )
        return data


class Methyldackel(ToolConfig):
    """Methyldackel methylation analysis tool configuration."""
    
    @model_validator(mode="before")
    @classmethod
    def set_defaults(cls, data: dict) -> dict:
        data = dict(data)
        data.setdefault("threads", 8)
        data.setdefault("options", CommandLineArguments(value=""))
        return data


# =============================================================================
# Complex Tool Configurations (multiple sub-tools)
# =============================================================================

class Bowtie2(BaseModel):
    """Bowtie2 alignment tool configuration."""
    
    align: ToolConfig = Field(
        default_factory=lambda: ToolConfig(
            threads=8, command_line_arguments=CommandLineArguments(value="--very-sensitive")
        ),
        description="Alignment configuration"
    )
    index: ToolConfig = Field(
        default_factory=lambda: ToolConfig(
            threads=8, command_line_arguments=CommandLineArguments(value="--threads {threads} --quiet")
        ),
        description="Index building configuration"
    )


class Samtools(BaseModel):
    """Samtools suite configuration."""
    
    sort: ToolConfig = Field(
        default_factory=lambda: ToolConfig(
            threads=8, command_line_arguments=CommandLineArguments(value="-@ {threads} -m 2G")
        ),
        description="Sort configuration"
    )
    index: ToolConfig = Field(
        default_factory=lambda: ToolConfig(threads=8, command_line_arguments=CommandLineArguments(value="")),
        description="Index configuration"
    )
    view: ToolConfig = Field(
        default_factory=lambda: ToolConfig(threads=8, command_line_arguments=CommandLineArguments(value="-bS")),
        description="View configuration"
    )


class Deeptools(BaseModel):
    """Deeptools suite configuration."""
    
    bam_coverage: ToolConfig = Field(
        default_factory=lambda: ToolConfig(
            threads=8, command_line_arguments=CommandLineArguments(
                value="--binSize 10 --normalizeUsing RPKM", 
                exclude={"-e", "--extendReads"}
            )
        ),
        description="BAM coverage configuration"
    )
    compute_matrix: ToolConfig = Field(
        default_factory=lambda: ToolConfig(threads=8, command_line_arguments=CommandLineArguments(value="")),
        description="Compute matrix configuration"
    )
    plot_heatmap: ToolConfig = Field(
        default_factory=lambda: ToolConfig(
            threads=8,
            command_line_arguments=CommandLineArguments(
                value="--colorMap RdYlBu --whatToShow 'heatmap and colorbar' --boxAroundHeatmaps no"
            ),
        ),
        description="Plot heatmap configuration"
    )



class Homer(BaseModel):
    """HOMER motif analysis suite configuration."""
    
    make_tag_directory: ToolConfig = Field(
        default_factory=lambda: ToolConfig(command_line_arguments=CommandLineArguments(value="")),
        description="Tag directory creation configuration"
    )
    make_bigwig: ToolConfig = Field(
        default_factory=lambda: ToolConfig(command_line_arguments=CommandLineArguments(value="")),
        description="BigWig creation configuration"
    )
    annotate_peaks: ToolConfig = Field(
        default_factory=lambda: ToolConfig(command_line_arguments=CommandLineArguments(value="")),
        description="Peak annotation configuration"
    )
    find_peaks: ToolConfig = Field(
        default_factory=lambda: ToolConfig(command_line_arguments=CommandLineArguments(value="")),
        description="Peak finding configuration"
    )
    find_motifs_genome: ToolConfig = Field(
        default_factory=lambda: ToolConfig(command_line_arguments=CommandLineArguments(value="")),
        description="Motif finding configuration"
    )


class Star(BaseModel):
    """STAR RNA-seq aligner configuration."""
    
    align: ToolConfig = Field(
        default_factory=lambda: ToolConfig(
            threads=8,
            command_line_arguments=CommandLineArguments(
                value="--quantMode TranscriptomeSAM GeneCounts --outSAMunmapped Within --outSAMattributes Standard"
            ),
        ),
        description="Alignment configuration"
    )
    index: ToolConfig = Field(
        default_factory=lambda: ToolConfig(
            threads=8,
            command_line_arguments=CommandLineArguments(value="--runThreadN {threads} --genomeSAindexNbases 14"),
        ),
        description="Index building configuration"
    )


class BcfTools(BaseModel):
    """BCFtools variant calling suite configuration."""
    
    call: ToolConfig = Field(
        default_factory=lambda: ToolConfig(
            threads=8, command_line_arguments=CommandLineArguments(value="-c -v -f GQ,DP")
        ),
        description="Variant calling configuration"
    )
    consensus: ToolConfig = Field(
        default_factory=lambda: ToolConfig(threads=8, command_line_arguments=CommandLineArguments(value="-s -c")),
        description="Consensus calling configuration"
    )

class Salmon(BaseModel):
    """Salmon quantification tool configuration."""

    quant: ToolConfig = Field(
        default_factory=lambda: ToolConfig(
            threads=8, command_line_arguments=CommandLineArguments(value="-p {threads} -o {output.dir}")
        ),
        description="Quantification configuration"
    )


# =============================================================================
# Single-purpose Tool Configurations
# =============================================================================

class Bamnado(BaseModel):
    """Bamnado coverage analysis tool configuration."""
    
    bam_coverage: ToolConfig = Field(
        default_factory=lambda: ToolConfig(
            threads=8, command_line_arguments=CommandLineArguments(value="--bin-size 10 --norm-method rpkm")
        ),
        description="BAM coverage analysis configuration"
    )


class Subread(BaseModel):
    """Subread feature counting tool configuration."""
    
    feature_counts: ToolConfig = Field(
        default_factory=lambda: ToolConfig(
            threads=16, command_line_arguments=CommandLineArguments(value="-p --countReadPairs")
        ),
        description="Feature counting configuration"
    )


class Lanceotron(BaseModel):
    """Lanceotron peak calling tool configuration."""
    
    call_peak: ToolConfig = Field(
        default_factory=lambda: ToolConfig(threads=8, command_line_arguments=CommandLineArguments(value="-c 0.5")),
        description="Peak calling configuration"
    )


class Picard(BaseModel):
    """Picard tools configuration."""
    
    mark_duplicates: ToolConfig = Field(
        default_factory=lambda: ToolConfig(
            threads=8, command_line_arguments=CommandLineArguments(value="REMOVE_DUPLICATES=true")
        ),
        description="Mark duplicates configuration"
    )


class Trimgalore(BaseModel):
    """Trim Galore adapter trimming tool configuration."""
    
    trim: ToolConfig = Field(
        default_factory=lambda: ToolConfig(
            threads=4, command_line_arguments=CommandLineArguments(value="--2colour 20")
        ),
        description="Trimming configuration"
    )


# =============================================================================
# Special Configuration Classes
# =============================================================================

class Macs(BaseModel):
    """MACS peak calling tool configuration with mode-specific defaults."""
    
    version: Literal[2, 3] = Field(default=3, description="MACS version to use")
    mode: MacsMode = Field(default=MacsMode.GENERIC, description="Analysis mode")
    call_peak: Optional[ToolConfig] = Field(default=None, description="Peak calling configuration")

    def _get_options_for_mode(self) -> str:
        """Get default options based on the analysis mode."""
        match self.mode:
            case MacsMode.ATAC:
                return "--nomodel --shift -100 --extsize 200"
            case MacsMode.CHIP_BROAD:
                return "--broad --nomodel --shift -100 --extsize 200"
            case MacsMode.CHIP_TF:
                return "--nomodel --shift -100 --extsize 200"
            case _:
                return "-f BAMPE"

    @model_validator(mode="after")
    def set_default_call_peak(self) -> "Macs":
        if self.call_peak is None:
            self.call_peak = ToolConfig(
                command_line_arguments=CommandLineArguments(value=self._get_options_for_mode())
            )
        return self


class Seacr(BaseModel):
    """SEACR peak calling tool configuration (non-ToolConfig based)."""
    
    threshold: float = Field(default=0.01, ge=0, le=1, description="Peak calling threshold")
    normalization: Optional[str] = Field(default="non", description="Normalization method")
    stringency: Literal["stringent", "relaxed"] = Field(
        default="stringent", description="Peak calling stringency"
    )


class FastqScreen(BaseModel):
    """FastqScreen tool configuration."""
    
    threads: int = Field(default=4, description="Number of threads to use")
    config: Path = Field(default=Path("fastq_screen_config.yaml"), description="Path to the FastqScreen config file")
    command_line_arguments: CommandLineArguments = Field(
        default_factory=lambda: CommandLineArguments()
    )


# =============================================================================
# Main Configuration Class
# =============================================================================

def get_assay_specific_tools(assay: Assay) -> list[type[BaseModel]]:
    """Get the list of tool classes appropriate for a given assay type."""
    generic_dna_tools = [
        Bowtie2,
        Bamnado,
        Deeptools,
        Homer,
        Lanceotron,
        Macs,
        Picard,
        Samtools,
        Trimgalore,
        Subread,
    ]
    
    tools = {
        Assay.ATAC: generic_dna_tools,
        Assay.CHIP: generic_dna_tools,
        Assay.CAT: generic_dna_tools + [Seacr],
        Assay.RNA: [Star, Deeptools, Samtools, Trimgalore, Salmon],
        Assay.SNP: [Bowtie2, Trimgalore, Samtools, BcfTools],
        Assay.METH: [Bowtie2, Methyldackel, Samtools, Picard],
        Assay.CRISPR: [Cutadapt, Bowtie2, Subread, Samtools],
    }
    return tools.get(assay, [])


class ThirdPartyToolsConfig(BaseModel):
    """Configuration for all third-party bioinformatics tools."""
    
    # Alignment tools
    bowtie2: Optional[Bowtie2] = Field(default=None, description="Bowtie2 aligner configuration")
    star: Optional[Star] = Field(default=None, description="STAR RNA-seq aligner configuration")

    # Processing tools
    samtools: Optional[Samtools] = Field(default=None, description="Samtools suite configuration")
    picard: Optional[Picard] = Field(default=None, description="Picard tools configuration")
    
    # Trimming tools
    cutadapt: Optional[Cutadapt] = Field(default=None, description="Cutadapt trimming configuration")
    trimgalore: Optional[Trimgalore] = Field(default=None, description="Trim Galore configuration")
    
    # Peak calling tools
    macs: Optional[Macs] = Field(default=None, description="MACS peak caller configuration")
    lanceotron: Optional[Lanceotron] = Field(default=None, description="Lanceotron peak caller configuration")
    lanceotron_mcc: Optional[LanceotronMCC] = Field(default=None, description="Lanceotron MCC configuration")
    seacr: Optional[Seacr] = Field(default=None, description="SEACR peak caller configuration")
    
    # Analysis tools
    deeptools: Optional[Deeptools] = Field(default=None, description="Deeptools suite configuration")
    homer: Optional[Homer] = Field(default=None, description="HOMER suite configuration")
    bamnado: Optional[Bamnado] = Field(default=None, description="Bamnado coverage analysis configuration")
    subread: Optional[Subread] = Field(default=None, description="Subread feature counting configuration")
    salmon: Optional[Salmon] = Field(default=None, description="Salmon quantification configuration")

    
    # Specialized tools
    methyldackel: Optional[Methyldackel] = Field(default=None, description="Methyldackel methylation analysis configuration")
    bcftools: Optional[BcfTools] = Field(default=None, description="BCFtools variant calling suite configuration")

    @classmethod
    def for_assay(cls, assay: Assay, **overrides) -> "ThirdPartyToolsConfig":
        """
        Factory method to create a configuration with defaults based on assay type.
        
        Args:
            assay: The assay type to configure tools for
            **overrides: Override configurations for specific tools
            
        Returns:
            ThirdPartyToolsConfig: Configured instance with assay-appropriate defaults
            
        Raises:
            ValueError: If no tools are configured for the given assay type
        """
        assay_tools = get_assay_specific_tools(assay)
        if not assay_tools:
            raise ValueError(f"No tools configured for assay {assay.value}")
        
        # Create default instances for the assay's tools
        defaults = {}
        for tool_class in assay_tools:
            tool_name = tool_class.__name__.lower()
            defaults[tool_name] = tool_class()
        
        # Merge with any user overrides
        defaults.update(overrides)
        
        return cls(**defaults)

    def get_configured_tools(self) -> dict[str, BaseModel]:
        """Get a dictionary of all configured (non-None) tools."""
        return {
            name: tool for name, tool in self.model_dump(exclude_none=True).items()
            if tool is not None
        }

    def has_tool(self, tool_name: str) -> bool:
        """Check if a specific tool is configured."""
        return getattr(self, tool_name, None) is not None