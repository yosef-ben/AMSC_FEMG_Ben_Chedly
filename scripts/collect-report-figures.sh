#!/usr/bin/env bash
# Copy the report-ready figures of the connectome chapter into report/images.
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p report/images

copy() { cp "benchmarks/$1" "report/images/$(basename "$1")"; }

copy 18_fisher_kolmogorov_1d_sensitivity/results/sensitivity.pdf
copy 18_fisher_kolmogorov_1d_sensitivity/results/time_step_study.pdf
copy 19_fisher_kolmogorov_fornari83/results/biomarker_comparison.pdf
copy 19_fisher_kolmogorov_fornari83/results/connectome_topology.pdf
copy 20_fisher_kolmogorov_alpha_sensitivity/results/alpha_sensitivity.pdf
copy 21_fisher_kolmogorov_corti83/results/regional_averages.pdf
copy 22_fisher_kolmogorov_sequential_performance/results/sequential_performance.pdf
copy 23_fisher_kolmogorov_diffusion_scaling/results/diffusion_scaling.pdf
copy 24_connectome_topology/results/connectome_regions.pdf
copy 24_connectome_topology/results/connectome_connectogram.pdf
copy 24_connectome_topology/results/connectome_views.pdf
copy 24_connectome_topology/results/lobe_connectivity.pdf
copy 25_connectome_seeding_vulnerability/results/seeding_vulnerability.pdf
copy 26_connectome_progression/results/activation_order.pdf
copy 27_connectome_seeding_patterns/results/seeding_patterns_expected.pdf

echo "Copied $(ls report/images | wc -l) figures into report/images"
