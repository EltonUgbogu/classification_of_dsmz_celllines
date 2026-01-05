#!/usr/bin/env python3
"""
Independence check script: validates that all input paths in config.yaml
point to locations under ~/dsmz/pipeline/data/

Outputs are allowed to be outside data/, but all inputs must be under data/.
"""

import yaml
import os
import sys
from pathlib import Path


def collect_strings(x):
    """Recursively collect all string values from nested dict/list structures."""
    if isinstance(x, dict):
        for v in x.values():
            yield from collect_strings(v)
    elif isinstance(x, list):
        for v in x:
            yield from collect_strings(v)
    elif isinstance(x, str):
        yield x


def is_output_path(path_key, path_value):
    """Determine if a path is an output (allowed outside data/) or an input."""
    output_indicators = [
        'unsup_root',
        'hvg_dir',
        'method_outdir',
        'agnostic_cluster_root',
        'tumour_nh_input_root',
        'hvg_expr',
        'tumour_nh_expr',
        'tumour_nh_map_rds',
        'hvg_joint_scores_tsv',
        'hvg_final_gene_list',
    ]
    return any(indicator in path_key.lower() for indicator in output_indicators)


def is_code_path(path_key, path_value):
    """Determine if a path points to code/scripts (not data inputs)."""
    code_indicators = [
        'base_functions',
        'functions_dir',
        '/R/',
        '/scripts/',
        '.R',
    ]
    return any(indicator in path_key.lower() or indicator in path_value.lower() 
               for indicator in code_indicators)


def check_config_independence(config_path):
    """Check that all input paths are under data/ directory."""
    pipeline_root = Path(config_path).parent.parent
    data_root = pipeline_root / "data"
    
    with open(config_path, 'r') as f:
        cfg = yaml.safe_load(f)
    
    bad_inputs = []
    all_paths = []
    
    for prof, prof_data in cfg.get("profiles", {}).items():
        paths_dict = prof_data.get("paths", {})
        
        for path_key, path_value in paths_dict.items():
            if not isinstance(path_value, str):
                continue
                
            all_paths.append((prof, path_key, path_value))
            
            # Skip output paths (they're allowed outside data/)
            if is_output_path(path_key, path_value):
                continue
            
            # Skip code/script paths (they're not data inputs)
            if is_code_path(path_key, path_value):
                continue
            
            # Check if it's an absolute path
            if path_value.startswith("/"):
                # Check if it's under data/ directory
                if "/dsmz/pipeline/data/" not in path_value:
                    bad_inputs.append((prof, path_key, path_value))
    
    # Print results
    print(f"Total paths checked: {len(all_paths)}")
    print(f"Non-data absolute input paths: {len(bad_inputs)}")
    
    if bad_inputs:
        print("\n❌ Found input paths outside data/ directory:")
        for prof, key, path in bad_inputs:
            print(f"  [{prof}] {key}: {path}")
        return False
    else:
        print("\n✅ All input paths are under data/ directory!")
        return True


if __name__ == "__main__":
    config_path = Path(__file__).parent.parent / "config" / "config.yaml"
    
    if not config_path.exists():
        print(f"Error: Config file not found at {config_path}", file=sys.stderr)
        sys.exit(1)
    
    success = check_config_independence(config_path)
    sys.exit(0 if success else 1)

