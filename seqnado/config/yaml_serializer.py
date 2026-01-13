"""YAML serialization utilities with support for computed field comments."""

from typing import Any, TextIO
from io import StringIO
from pathlib import Path

from ruamel.yaml import YAML
from pydantic import BaseModel

from seqnado.config.mixins import ComputedFieldOverrideMixin


def dump_config_to_yaml(
    config: BaseModel,
    stream: TextIO | Path | None = None,
    add_computed_field_comments: bool = True
) -> str | None:
    """Serialize a pydantic config model to YAML with optional computed field comments.
    
    Args:
        config: The pydantic model instance to serialize
        stream: Optional file path or text stream to write to. If None, returns YAML string.
        add_computed_field_comments: If True, add comments to computed fields
        
    Returns:
        YAML string if stream is None, otherwise None (writes to stream)
        
    Example:
        >>> config = SeqnadoConfig(...)
        >>> yaml_str = dump_config_to_yaml(config)
        >>> # Or write to file:
        >>> dump_config_to_yaml(config, Path("config.yml"))
    """
    yaml = YAML()
    yaml.default_flow_style = False
    yaml.width = 4096  # Prevent line wrapping
    yaml.preserve_quotes = True
    
    # Get the model dict, excluding private _explicit_overrides field
    config_dict = config.model_dump(
        mode="json", 
        exclude_none=False,
        exclude={'_explicit_overrides'}
    )
    
    # Get computed field names if model uses ComputedFieldOverrideMixin
    computed_fields = set()
    if isinstance(config, ComputedFieldOverrideMixin):
        computed_fields = getattr(config.__class__, '__computed_fields__', set())
    
    # Add comments to computed fields
    if add_computed_field_comments and computed_fields:
        # Convert to CommentedMap for comment support
        from ruamel.yaml.comments import CommentedMap
        
        def add_comments_recursive(data: Any, parent_key: str = "") -> Any:
            """Recursively add comments to computed fields in nested structures."""
            if isinstance(data, dict):
                commented = CommentedMap(data)
                for key, value in data.items():
                    # Check if this key is a computed field
                    if key in computed_fields:
                        # Add comment before the key
                        commented.yaml_set_comment_before_after_key(
                            key,
                            before="computed - override if needed",
                            indent=0
                        )
                    # Recursively process nested dicts/lists
                    if isinstance(value, (dict, list)):
                        commented[key] = add_comments_recursive(value, key)
                return commented
            elif isinstance(data, list):
                return [add_comments_recursive(item, parent_key) for item in data]
            else:
                return data
        
        config_dict = add_comments_recursive(config_dict)
    
    # Handle output
    if stream is None:
        # Return as string
        string_stream = StringIO()
        yaml.dump(config_dict, string_stream)
        return string_stream.getvalue()
    elif isinstance(stream, Path):
        # Write to file path
        with open(stream, 'w') as f:
            yaml.dump(config_dict, f)
        return None
    else:
        # Write to provided stream
        yaml.dump(config_dict, stream)
        return None


def load_config_from_yaml(path: Path) -> dict[str, Any]:
    """Load a YAML config file into a dictionary.
    
    Args:
        path: Path to the YAML file
        
    Returns:
        Dictionary containing the config data
    """
    yaml = YAML()
    with open(path, 'r') as f:
        return yaml.load(f)
