#!/usr/bin/env python3
"""
HBV Pipeline Configuration Parser

Parses YAML configuration files and exports values as environment variables
or shell-compatible format for use by bash scripts.

Usage:
    # Get single value
    python config_parser.py config.yaml --get project.work_dir

    # Export all as shell variables
    python config_parser.py config.yaml --export

    # Validate configuration
    python config_parser.py config.yaml --validate

    # As Python module
    from scripts.utils.config_parser import load_config, get_nested
    config = load_config("config.yaml")
    work_dir = get_nested(config, "project.work_dir")
"""

import argparse
import os
import sys
from pathlib import Path
from typing import Any, Dict, Optional

try:
    import yaml
except ImportError:
    print("Error: PyYAML is required. Install with: pip install pyyaml", file=sys.stderr)
    sys.exit(1)


def load_config(config_path: str) -> Dict[str, Any]:
    """
    Load configuration from YAML file.
    
    Args:
        config_path: Path to YAML configuration file
        
    Returns:
        Dictionary containing configuration values
        
    Raises:
        FileNotFoundError: If config file doesn't exist
        yaml.YAMLError: If YAML parsing fails
    """
    path = Path(config_path)
    if not path.exists():
        raise FileNotFoundError(f"Configuration file not found: {config_path}")
    
    with open(path, 'r') as f:
        config = yaml.safe_load(f)
    
    if config is None:
        return {}
    
    return config


def get_nested(config: Dict[str, Any], key_path: str, default: Any = None) -> Any:
    """
    Get a nested value from config using dot notation.
    
    Args:
        config: Configuration dictionary
        key_path: Dot-separated path (e.g., "project.work_dir")
        default: Default value if key not found
        
    Returns:
        Value at the specified path, or default if not found
        
    Examples:
        >>> config = {"project": {"work_dir": "/path/to/work"}}
        >>> get_nested(config, "project.work_dir")
        '/path/to/work'
        >>> get_nested(config, "project.missing", "default")
        'default'
    """
    keys = key_path.split('.')
    value = config
    
    for key in keys:
        if isinstance(value, dict) and key in value:
            value = value[key]
        else:
            return default
    
    return value


def flatten_config(config: Dict[str, Any], prefix: str = "") -> Dict[str, str]:
    """
    Flatten nested config into flat dictionary with dot-notation keys.
    
    Args:
        config: Nested configuration dictionary
        prefix: Current key prefix (used in recursion)
        
    Returns:
        Flat dictionary with dot-notation keys and string values
        
    Examples:
        >>> config = {"project": {"work_dir": "/path"}}
        >>> flatten_config(config)
        {'project.work_dir': '/path'}
    """
    flat = {}
    
    for key, value in config.items():
        full_key = f"{prefix}.{key}" if prefix else key
        
        if isinstance(value, dict):
            flat.update(flatten_config(value, full_key))
        else:
            # Convert to string for shell compatibility
            if value is None:
                flat[full_key] = ""
            elif isinstance(value, bool):
                flat[full_key] = "true" if value else "false"
            else:
                flat[full_key] = str(value)
    
    return flat


def to_shell_var_name(key_path: str) -> str:
    """
    Convert dot-notation key to shell variable name.
    
    Args:
        key_path: Dot-notation key (e.g., "project.work_dir")
        
    Returns:
        Shell-compatible variable name (e.g., "HBV_PROJECT_WORK_DIR")
        
    Examples:
        >>> to_shell_var_name("project.work_dir")
        'HBV_PROJECT_WORK_DIR'
    """
    return "HBV_" + key_path.upper().replace(".", "_").replace("-", "_")


def export_as_shell(config: Dict[str, Any]) -> str:
    """
    Export config as shell variable assignments.
    
    Args:
        config: Configuration dictionary
        
    Returns:
        String of shell export statements
    """
    flat = flatten_config(config)
    lines = []
    
    for key, value in sorted(flat.items()):
        var_name = to_shell_var_name(key)
        # Escape single quotes in value
        escaped_value = str(value).replace("'", "'\"'\"'")
        lines.append(f"export {var_name}='{escaped_value}'")
    
    return "\n".join(lines)


def validate_config(config: Dict[str, Any]) -> tuple[bool, list[str]]:
    """
    Validate configuration for required fields.
    
    Args:
        config: Configuration dictionary
        
    Returns:
        Tuple of (is_valid, list of error messages)
    """
    errors = []
    
    # Required paths
    required_paths = [
        ("project.work_dir", "Working directory"),
        ("reference.hbv_unified", "HBV unified reference"),
    ]
    
    for key_path, description in required_paths:
        value = get_nested(config, key_path)
        if not value or value == "/path/to/your/hbv_project":
            errors.append(f"Missing or placeholder: {description} ({key_path})")
    
    # Validate numeric ranges
    threads = get_nested(config, "resources.threads", 1)
    try:
        if int(threads) < 1:
            errors.append(f"resources.threads must be >= 1, got {threads}")
    except (ValueError, TypeError):
        errors.append(f"resources.threads must be an integer, got {threads}")
    
    # Check reference files exist (optional validation)
    ref_paths = [
        "reference.hbv_panel",
        "reference.hbv_unified",
        "reference.human_ref",
    ]
    
    for key_path in ref_paths:
        path = get_nested(config, key_path)
        if path and not path.startswith("/path/to"):
            if not Path(path).exists():
                errors.append(f"Reference file not found: {path} ({key_path})")
    
    return len(errors) == 0, errors


def print_config_summary(config: Dict[str, Any]) -> None:
    """Print a human-readable config summary."""
    print("=" * 60)
    print("HBV Pipeline Configuration Summary")
    print("=" * 60)
    
    sections = [
        ("Project", [
            ("project.work_dir", "Work Directory"),
            ("project.raw_fastq_dir", "Raw FASTQ Directory"),
            ("project.results_dir", "Results Directory"),
        ]),
        ("Reference", [
            ("reference.hbv_unified", "HBV Reference"),
            ("reference.human_ref", "Human Reference"),
        ]),
        ("Resources", [
            ("resources.threads", "Threads"),
            ("resources.threads_per_job", "Threads/Job"),
        ]),
        ("QC", [
            ("qc.min_read_length", "Min Read Length"),
            ("qc.max_read_length", "Max Read Length"),
            ("qc.min_quality", "Min Quality"),
        ]),
    ]
    
    for section_name, fields in sections:
        print(f"\n{section_name}:")
        for key_path, label in fields:
            value = get_nested(config, key_path, "not set")
            print(f"  {label}: {value}")
    
    print("\n" + "=" * 60)


def main():
    parser = argparse.ArgumentParser(
        description="HBV Pipeline Configuration Parser",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    parser.add_argument(
        "config",
        help="Path to YAML configuration file"
    )
    
    parser.add_argument(
        "--get",
        metavar="KEY",
        help="Get single value using dot notation (e.g., project.work_dir)"
    )
    
    parser.add_argument(
        "--export",
        action="store_true",
        help="Export all config as shell variable assignments"
    )
    
    parser.add_argument(
        "--validate",
        action="store_true",
        help="Validate configuration and report errors"
    )
    
    parser.add_argument(
        "--summary",
        action="store_true",
        help="Print configuration summary"
    )
    
    parser.add_argument(
        "--json",
        action="store_true",
        help="Output as JSON (for --get with complex values)"
    )
    
    args = parser.parse_args()
    
    try:
        config = load_config(args.config)
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except yaml.YAMLError as e:
        print(f"Error parsing YAML: {e}", file=sys.stderr)
        sys.exit(1)
    
    if args.get:
        value = get_nested(config, args.get)
        if value is None:
            print(f"Key not found: {args.get}", file=sys.stderr)
            sys.exit(1)
        if args.json:
            import json
            print(json.dumps(value))
        else:
            print(value)
    
    elif args.export:
        print(export_as_shell(config))
    
    elif args.validate:
        is_valid, errors = validate_config(config)
        if is_valid:
            print("Configuration is valid!")
            sys.exit(0)
        else:
            print("Configuration errors:", file=sys.stderr)
            for error in errors:
                print(f"  - {error}", file=sys.stderr)
            sys.exit(1)
    
    elif args.summary:
        print_config_summary(config)
    
    else:
        # Default: print summary
        print_config_summary(config)


if __name__ == "__main__":
    main()
