#!/usr/bin/env bash
set -euo pipefail

benchmarks=(
  02_star_constant
  03_star_linear
  04_star_sine
  05_star_radial_decay
  06_star_localized
)

for benchmark in "${benchmarks[@]}"; do
  echo
  echo "=== Running ${benchmark} ==="
  scripts/run-benchmark.sh "$benchmark"
done

echo
echo "All benchmarks completed."
