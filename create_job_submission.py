#!/usr/bin/env python3
"""Render or submit a DNAnexus workflow run from a YAML config.

Example config:

workflow: "/workflows/fine_mapping/fine_mapping"
project: "weinstock lab"
destination: "/analysis/fine_mapping/"
name: "susie-finemap-run"
priority: "normal"
cost_limit: 20
tags:
  - fine_mapping
  - susie
inputs:
  summary_stats_parquet: "file-xxxx"
  step2_chunk_manifest: "file-yyyy"
  covariates_file: "file-zzzz"
  plink2_binary: "file-aaaa"
  p_threshold: 5e-8
  window_kb: 500
  max_causal: 10
  min_abs_corr: 0.5
  coverage: 0.95
  max_iter: 500
  tol: 0.001
  refine: true
  estimate_residual_variance: true
  threads: 4
  r_docker: "rocker/r-ver:4.3.0"
"""

import argparse
import json
import shlex
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

import yaml


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create or launch a DNAnexus job for fine_mapping.WDL from YAML config."
    )
    parser.add_argument("config", help="Path to YAML config file")
    parser.add_argument(
        "--run",
        action="store_true",
        help="Execute the dx run command instead of printing it",
    )
    return parser.parse_args()


def load_config(path: str) -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    if not isinstance(config, dict):
        raise ValueError("Config file must contain a YAML mapping at the top level.")
    return config


def require_mapping(config: Dict[str, Any], key: str) -> Dict[str, Any]:
    value = config.get(key)
    if not isinstance(value, dict):
        raise ValueError(f"Config key '{key}' must be a mapping.")
    return value


def normalize_dx_target(value: str, project: Optional[str]) -> str:
    if not isinstance(value, str):
        raise ValueError("Workflow and destination values must be strings.")
    if not project or ":" in value or value.startswith("project-"):
        return value
    if value.startswith("/"):
        return f"{project}:{value}"
    return value


def format_input_value(value: Any) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (int, float, str)):
        return str(value)
    if isinstance(value, (list, dict)):
        return json.dumps(value)
    raise ValueError(f"Unsupported input value type: {type(value).__name__}")


def build_dx_command(config: Dict[str, Any]) -> List[str]:
    required_top_level = ["workflow", "destination", "inputs"]
    missing = [key for key in required_top_level if key not in config]
    if missing:
        raise ValueError(f"Missing required config keys: {', '.join(missing)}")

    inputs = require_mapping(config, "inputs")
    if not inputs:
        raise ValueError("Config key 'inputs' cannot be empty.")

    project = config.get("project")
    input_prefix = config.get("input_prefix", "stage-common")

    workflow = normalize_dx_target(str(config["workflow"]), project)
    destination = normalize_dx_target(str(config["destination"]), project)

    cmd = ["dx", "run", workflow]

    for name, value in inputs.items():
        cmd.append(f"-i{input_prefix}.{name}={format_input_value(value)}")

    if "name" in config:
        cmd.extend(["--name", str(config["name"])])
    if "priority" in config:
        cmd.extend(["--priority", str(config["priority"])])
    if "cost_limit" in config:
        cmd.extend(["--cost-limit", str(config["cost_limit"])])

    tags = config.get("tags", [])
    if isinstance(tags, str):
        tags = [tags]
    if not isinstance(tags, list):
        raise ValueError("Config key 'tags' must be a string or list of strings.")
    for tag in tags:
        cmd.extend(["--tag", str(tag)])

    cmd.extend(["--destination", destination])

    if config.get("brief", True):
        cmd.append("--brief")
    if config.get("yes", True):
        cmd.append("-y")

    extra_args = config.get("extra_dx_run_args", [])
    if isinstance(extra_args, str):
        extra_args = [extra_args]
    if not isinstance(extra_args, list):
        raise ValueError("Config key 'extra_dx_run_args' must be a string or list of strings.")
    cmd.extend(str(arg) for arg in extra_args)

    return cmd


def shell_join(cmd: List[str]) -> str:
    return " ".join(shlex.quote(part) for part in cmd)


def main() -> int:
    args = parse_args()
    config_path = Path(args.config)
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")

    config = load_config(str(config_path))
    cmd = build_dx_command(config)

    if not args.run:
        print(shell_join(cmd))
        return 0

    print(shell_join(cmd), file=sys.stderr)
    subprocess.run(cmd, check=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
