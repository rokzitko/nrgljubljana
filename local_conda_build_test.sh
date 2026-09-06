#!/usr/bin/env bash
# Note: run this from a clean environment (e.g. do "module purge")

set -euo pipefail

repo_root="$(git -C "$(dirname -- "${BASH_SOURCE[0]}")" rev-parse --show-toplevel)"

# Initialize Conda for this non-interactive shell.
eval "$(conda shell.bash hook)"
conda activate base
export CONDA_CHANNEL_PRIORITY=strict
conda install --override-channels -c conda-forge conda-build

# Build this checkout's committed HEAD; uncommitted source changes are excluded.
export NRGLJUBLJANA_CONDA_GIT_URL="file://${repo_root}"
NRGLJUBLJANA_CONDA_GIT_REV="$(git -C "${repo_root}" rev-parse HEAD)"
export NRGLJUBLJANA_CONDA_GIT_REV

variants="{nrgljubljana_build_tests: ['ON'], nrgljubljana_build_jobs: ['2'], nrgljubljana_test_jobs: ['2']}"

conda render --override-channels -c conda-forge \
  --variants "$variants" "${repo_root}/recipe"

conda build --override-channels -c conda-forge --no-anaconda-upload \
  --variants "$variants" "${repo_root}/recipe"
