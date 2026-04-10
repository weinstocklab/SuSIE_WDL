#!/usr/bin/env bash
set -euo pipefail

uv run dx run workflow-J7209P8Jy2bf473XvX6ZBXZ3 \
    -istage-common.summary_stats_parquet=file-J6xQ4K8JJ6JJzQqBky2Kyzjp \
    -istage-common.step2_chunk_manifest=file-J62FjK0Jy2bVPbyQFKP3p886 \
    -istage-common.covariates_file=file-J6vq5X8Jy2bygQz9068JYzxJ \
    -istage-common.plink2_binary=file-Gzq7v9jJy2bQ9gF5vJXj1Xp3 \
    -istage-common.min_maf=0.001 \
    --name fine-mapping-gte0-residualized-passengers \
    --priority normal \
    --cost-limit 20 \
    --tag fine_mapping \
    --destination 'weinstock lab:/somatic/fine_mapping_results/gte0_residualized_passengers/' \
    --brief \
    -y
