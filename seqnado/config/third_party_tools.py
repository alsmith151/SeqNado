
import shlex
from pydantic import BaseModel, Field, field_validator, model_validator
from typing import Optional, Literal, Any
from enum import Enum

from seqnado import Assay

class Options(BaseModel):
    """Pydantic model for validating CLI options."""
    
    value: str = Field(default="", description="CLI options string")
    exclude: set[str] = Field(default_factory=set, description="Options to exclude from the final command")

    @field_validator("value", mode="before")
    @classmethod
    def validate_options_syntax(cls, v: str) -> str:
        if not isinstance(v, str):
            raise ValueError("Value must be a string")
        
        try:
            # Check that it's safely tokenizable like a CLI
            shlex.split(v)
        except ValueError as e:
            raise ValueError(f"Invalid CLI option string: {e}")
    
        # Optional: Reject unsafe characters
        if any(c in v for c in ["\n", "\r", "\x00"]):
            raise ValueError("Value string contains unsafe characters")

        return v

    @model_validator(mode="after")
    def ensure_leading_space(self) -> "Options":
        """Ensure options have a leading space if not empty."""
        if self.value and not self.value.startswith(" "):
            self.value = " " + self.value
        return self
    
    def _filter_excluded_options(self, options_str: str) -> str:
        """Remove excluded options from the options string."""
        if not self.exclude:
            return options_str
            
        try:
            tokens = shlex.split(options_str.strip())
            filtered_tokens = []
            i = 0
            
            while i < len(tokens):
                token = tokens[i]
                
                # Check if this token is an excluded option
                excluded = False
                for exclude_pattern in self.exclude:
                    if token == exclude_pattern or token.startswith(exclude_pattern + "="):
                        excluded = True
                        # If it's a flag with separate value, skip the next token too
                        if (token == exclude_pattern and 
                            i + 1 < len(tokens) and 
                            not tokens[i + 1].startswith("-")):
                            i += 1  # Skip the value token
                        break
                
                if not excluded:
                    filtered_tokens.append(token)
                
                i += 1
            
            return " " + " ".join(filtered_tokens) if filtered_tokens else ""
        except ValueError:
            # If parsing fails, return original
            return options_str
    
    def __str__(self) -> str:
        """Return options without leading space for display, with exclusions applied."""
        filtered = self._filter_excluded_options(self.value)
        return filtered.strip()
    
    @property
    def raw(self) -> str:
        """Return raw options string with leading space for CLI usage, with exclusions applied."""
        return self._filter_excluded_options(self.value)
    
    @property
    def raw_unfiltered(self) -> str:
        """Return raw options string without applying exclusions."""
        return self.value
    
    def model_dump_string(self) -> str:
        """Serialize to string for easy CLI usage, with exclusions applied."""
        return self.raw
    
    @classmethod
    def __get_pydantic_json_schema__(cls, core_schema, handler):
        """Make this serialize as a string in JSON schema."""
        json_schema = handler(core_schema)
        json_schema.update(type="string")
        return json_schema
    
    @classmethod
    def __get_pydantic_core_schema__(cls, source_type, handler):
        """Define core schema for serialization."""
        from pydantic_core import core_schema
        
        def serialize_options(instance, _info):
            if isinstance(instance, Options):
                return instance.raw  # Use filtered version for serialization
            return str(instance)
            
        def validate_options(value):
            if isinstance(value, Options):
                return value
            return cls(value=value)
        
        return core_schema.no_info_after_validator_function(
            validate_options,
            core_schema.str_schema(),
            serialization=core_schema.to_string_ser_schema(
                when_used='json'
            )
        )

class ToolConfig(BaseModel):
    """Base model for CLI-based tools."""
    threads: int = Field(default=1, ge=1)
    options: Options = Field(default_factory=Options)


    @field_validator("options", mode="before")
    @classmethod
    def ensure_options(cls, v):
        if isinstance(v, Options):
            return v
        if isinstance(v, str):
            return Options(value=v)
        if isinstance(v, dict):
            return Options(**v)
        raise TypeError("options must be a string, dict, or Options instance")

    @model_validator(mode="after")
    def coerce_options(self) -> "ToolConfig":
        if not isinstance(self.options, Options):
            if isinstance(self.options, str):
                self.options = Options(value=self.options)
            elif isinstance(self.options, dict):
                self.options = Options(**self.options)
            else:
                raise TypeError("options must be a string, dict, or Options instance")
        return self
    




class Bowtie2(BaseModel):
    """Configuration for Bowtie2 tool."""
    align: ToolConfig = Field(
        default_factory=lambda: ToolConfig(
            threads=8,
            options=Options(value="--very-sensitive")
        )
    )
    index: ToolConfig = Field(
        default_factory=lambda: ToolConfig(
            threads=8,
            options=Options(value="--threads {threads} --quiet")
        )
    )

class Subread(BaseModel):
    """Configuration for Subread tool."""
    feature_counts: ToolConfig = Field(
        default_factory=lambda: ToolConfig(
            threads=16,
            options=Options(value="-p --countReadPairs")
        )
    )
    

class Bamnado(BaseModel):
    """Configuration for Bamnado tool."""
    bam_coverage: ToolConfig = Field(
        default_factory=lambda: ToolConfig(
            threads=8,
            options=Options(value="--bin-size 10 --norm-method rpkm")
        )
    )

class Deeptools(BaseModel):
    """Configuration for Deeptools tool."""
    bam_coverage: ToolConfig = Field(
        default_factory=lambda: ToolConfig(
            threads=8,
            options=Options(value="--binSize 10 --normalizeUsing RPKM")
        )
    )
    compute_matrix: ToolConfig = Field(
        default_factory=lambda: ToolConfig(
            threads=8,
            options=Options(value="")
        )
    )
    plot_heatmap: ToolConfig = Field(
        default_factory=lambda: ToolConfig(
            threads=8,
            options=Options(value="--colorMap RdYlBu --whatToShow 'heatmap and colorbar'")
        )
    )


class Homer(BaseModel):
    make_tag_directory: ToolConfig = Field(
        default_factory=lambda: ToolConfig(
            options=Options(value="")
        )
    )
    make_bigwig: ToolConfig = Field(
        default_factory=lambda: ToolConfig(
            options=Options(value="")
        )
    )
    annotate_peaks: ToolConfig = Field(
        default_factory=lambda: ToolConfig(
            options=Options(value="")
        )
    )
    find_peaks: ToolConfig = Field(
        default_factory=lambda: ToolConfig(
            options=Options(value="")
        )
    )

class Lanceotron(BaseModel):
    """Configuration for Lanceotron tool."""
    call_peak: ToolConfig = Field(
        default_factory=lambda: ToolConfig(
            threads=8,
            options=Options(value="-c 0.5")
        )
    )

class MacsMode(Enum):
    """Enum for MACS modes."""
    ATAC = "atac"
    CHIP_BROAD = "chip_broad"
    CHIP_TF = "chip_tf"
    GENERIC = "generic"



class Macs(BaseModel):
    """Configuration for MACS tool."""
    version: Literal[2, 3] = 3
    mode: MacsMode = MacsMode.GENERIC
    call_peak: ToolConfig = Field(default=None)

    def _get_options_for_mode(self) -> str:
        """Get options string based on the mode."""
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
        """Set default call_peak configuration based on mode if not provided."""
        if self.call_peak is None:
            self.call_peak = ToolConfig(
                options=Options(value=self._get_options_for_mode())
            )
        return self

class Picard(BaseModel):
    """Configuration for Picard tool."""
    mark_duplicates: ToolConfig = Field(
        default_factory=lambda: ToolConfig(
            threads=8,
            options=Options(value="REMOVE_DUPLICATES=true")
        )
    )

class Samtools(BaseModel):
    """Configuration for Samtools tool."""
    sort: ToolConfig = Field(
        default_factory=lambda: ToolConfig(
            threads=8,
            options=Options(value="-@ {threads} -m 2G")
        )
    )
    index: ToolConfig = Field(
        default_factory=lambda: ToolConfig(
            threads=8,
            options=Options(value="")
        )
    )
    view: ToolConfig = Field(
        default_factory=lambda: ToolConfig(
            threads=8,
            options=Options(value="-bS")
        )
    )

class Trimgalore(BaseModel):
    """Configuration for Trim Galore tool."""
    trim: ToolConfig = Field(
        default_factory=lambda: ToolConfig(
            threads=4,
            options=Options(value="--2colour 20")
        )
    )

class Seacr(BaseModel):
    """Configuration for Seacr tool."""
    threshold: float = Field(default=0.01)
    normalization: str | None = Field(
        default="non",
        description="Normalization method to use, e.g., 'RPKM', 'CPM', etc."
    )
    stringency: Literal["stringent", "relaxed"] = Field(
        default="stringent",
        description="Stringency level to use, e.g., 'stringent', 'relaxed', etc."
    )


class CutadaptMode(Enum):
    """Enum for Cutadapt modes."""
    DEFAULT = "default"
    CRISPR = "crispr"


class Cutadapt(ToolConfig):
    """Configuration for Cutadapt tool."""
    mode: CutadaptMode = Field(default=CutadaptMode.CRISPR)
    
    def __init__(self, **data):
        if 'options' not in data:
            data['options'] = Options(value="-g 'ACACCG' --cut 0 -l 20 -m 20 -M 20 --discard-untrimmed")
        if 'threads' not in data:
            data['threads'] = 4
        super().__init__(**data)


class LanceotronMCC(ToolConfig):
    """Configuration for Lanceotron MCC tool."""
    
    def __init__(self, **data):
        if 'options' not in data:
            data['options'] = Options(value="--max-peak-size 2000 --no-mcc-filter-enriched-regions --no-mcc-filter-background")
        if 'threads' not in data:
            data['threads'] = 8
        super().__init__(**data)


class Methyldackel(ToolConfig):
    """Configuration for Methyldackel tool."""
    
    def __init__(self, **data):
        if 'options' not in data:
            data['options'] = Options(value="")
        if 'threads' not in data:
            data['threads'] = 8
        super().__init__(**data)

class Star(BaseModel):
    """Configuration for STAR aligner."""
    align: ToolConfig = Field(
        default_factory=lambda: ToolConfig(
            threads=8,
            options=Options(value="--quantMode TranscriptomeSAM GeneCounts --outSAMunmapped Within --outSAMattributes Standard")
        )
    )
    index: ToolConfig = Field(
        default_factory=lambda: ToolConfig(
            threads=8,
            options=Options(value="--runThreadN {threads} --genomeSAindexNbases 14")
        )
    )


class BcfTools(BaseModel):
    """Configuration for BCFtools."""
    call: ToolConfig = Field(
        default_factory=lambda: ToolConfig(
            threads=8,
            options=Options(value="-c -v -f GQ,DP")
        )
    )
    consensus: ToolConfig = Field(
        default_factory=lambda: ToolConfig(
            threads=8,
            options=Options(value="-s -c")
        )
    )



def get_assay_specific_tools(assay: Assay) -> list[BaseModel]:
    """Get the specific tools configuration based on the assay type."""

    generic_dna_tools = [
        Bowtie2, Bamnado, Deeptools, Homer, Lanceotron, Macs, Picard, Samtools, Trimgalore, Subread
    ]


    tools = {
        Assay.ATAC: generic_dna_tools,
        Assay.CHIP: generic_dna_tools,
        Assay.CAT: generic_dna_tools + [Seacr],
        Assay.RNA: [Star, Deeptools, Samtools, Trimgalore],
        Assay.SNP: [Bowtie2, Trimgalore, Samtools, BcfTools],
        Assay.METH: [Bowtie2, Methyldackel, Samtools, Picard],
        Assay.CRISPR: [Cutadapt, Bowtie2, Subread, Samtools],
    }

    return tools.get(assay)

class ThirdPartyToolsConfig(BaseModel):
    """Configuration for third-party tools used in SeqNado."""

    assay: Assay
    bowtie2: Bowtie2 | None = None
    bamnado: Bamnado | None = None
    deeptools: Deeptools | None = None
    homer: Homer | None = None
    lanceotron: Lanceotron | None = None
    macs: Macs | None = None
    picard: Picard | None = None
    samtools: Samtools | None = None
    subread: Subread | None = None
    trimgalore: Trimgalore | None = None
    seacr: Seacr | None = None
    cutadapt: Cutadapt | None = None
    lanceotronmcc: LanceotronMCC | None = None
    methyldackel: Methyldackel | None = None
    star: Star | None = None
    bcftools: BcfTools | None = None


    def model_post_init(self, context: Any) -> None:
        """Set default tools based on assay type if not provided."""

        assay_tools = get_assay_specific_tools(self.assay)
        if assay_tools is None:
            raise ValueError(f"No tools configured for assay {self.assay.value}")
        
        for tool_class in assay_tools:
            tool_name = tool_class.__name__.lower()
            if getattr(self, tool_name) is None:
                setattr(self, tool_name, tool_class())



        










