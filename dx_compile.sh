#!/usr/bin/env bash
set -euo pipefail

dxCompiler="${DX_COMPILER_JAR:-dxCompiler-2.15.0.jar}"
project="${DX_PROJECT:-weinstock lab}"
base_folder="${DX_FOLDER:-/workflows}"

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ ! -f "$dxCompiler" && -f "${repo_dir}/${dxCompiler}" ]]; then
    dxCompiler="${repo_dir}/${dxCompiler}"
fi

if [[ ! -f "$dxCompiler" ]]; then
    echo "dxCompiler jar not found: ${dxCompiler}" >&2
    echo "Set DX_COMPILER_JAR or place the jar next to this script." >&2
    exit 1
fi

shopt -s nullglob

if [[ "$#" -gt 0 ]]; then
    wdls=("$@")
else
    wdls=("${repo_dir}"/*.WDL)
fi

if [[ "${#wdls[@]}" -eq 0 ]]; then
    echo "No WDL files found to compile." >&2
    exit 1
fi

for wdl in "${wdls[@]}"; do
    if [[ ! -f "$wdl" ]]; then
        echo "Missing WDL: ${wdl}" >&2
        exit 1
    fi

    workflow_name="$(basename "${wdl}" .WDL)"
    dest_folder="${base_folder%/}/${workflow_name}/"

    echo "Compiling ${wdl} -> ${project}:${dest_folder}"
    java -jar "$dxCompiler" compile "$wdl" \
        -project "${project}" \
        -folder "${dest_folder}" \
        -f
done
