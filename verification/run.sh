#!/usr/bin/env bash
# Build and run the solver verification suite; nonzero exit on any failure.
set -euo pipefail
cd "$(dirname "$0")/.."
cmake --build build-release --target verify_fisher_kolmogorov
./build-release/verify_fisher_kolmogorov verification/results
