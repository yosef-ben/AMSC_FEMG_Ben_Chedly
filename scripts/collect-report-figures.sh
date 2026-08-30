#!/usr/bin/env bash
# Stage the report figures into report/images. The folder is rebuilt from
# scratch, so after a run it holds exactly the set the report uses, one
# file per figure, every one copied from the results of its benchmark.
set -euo pipefail
cd "$(dirname "$0")/.."
rm -rf report/images
mkdir -p report/images

copy() { cp "benchmarks/$1" "report/images/$(basename "$1")"; }

# Chapter 3: validation benchmarks.
copy 02_star_constant/results/star_constant.png
copy 03_star_linear/results/star_linear.png
copy 04_star_sine/results/star_sine.png
copy 04_star_sine/results/space_convergence.pdf
copy 05_star_radial_decay/results/star_radial_decay.png
copy 07_spectral_star/results/star_eigenmodes_00_05.png
copy 08_spectral_graphene/results/graphene_eigenmodes_00_05.png
copy 08_spectral_graphene/results/graphene_eigenmodes_06_09.png
copy 09_spectral_tree/results/tree_eigenmodes_00_07.png
copy 09_spectral_tree/results/tree_eigenmodes_08_15.png
copy 10_heat_graphene_eigenmode/results/graphene_eigenmode_time_snapshots.png
copy 10_heat_graphene_eigenmode/results/graphene_eigenmode_l2_decay.png
copy 11_heat_tree_eigenmode/results/tree_eigenmode_time_snapshots.png
copy 11_heat_tree_eigenmode/results/tree_eigenmode_l2_decay.png
copy 12_heat_graphene_eigenmode_time_convergence/results/time_convergence.png
copy 13_spectral_comparison/results/spectral_comparison.png

# Chapter 4: the connectome study.
copy 18_fisher_kolmogorov_1d_sensitivity/results/sensitivity.pdf
copy 18_fisher_kolmogorov_1d_sensitivity/results/time_step_study.pdf
copy 19_fisher_kolmogorov_fornari83/results/biomarker_comparison.pdf
copy 19_fisher_kolmogorov_fornari83/results/connectome_topology.pdf
copy 22_fisher_kolmogorov_sequential_performance/results/sequential_performance.pdf
copy 23_fisher_kolmogorov_diffusion_scaling/results/diffusion_scaling.pdf
copy 24_connectome_topology/results/connectome_regions.pdf
copy 24_connectome_topology/results/connectome_connectogram.pdf
copy 24_connectome_topology/results/connectome_views.pdf
copy 24_connectome_topology/results/lobe_connectivity.pdf
copy 27_connectome_seeding_patterns/results/seeding_patterns_expected.pdf
copy 27_connectome_seeding_patterns/results/seeding_patterns_regional_expected.pdf

# Chapter 5: the performance analysis.
copy 28_sequential_ordering_study/results/performance_ordering.pdf
copy 29_static_condensation/results/performance_condensation.pdf
copy 32_hybrid_condensation/results/performance_scaling.pdf

echo "Staged $(ls report/images | wc -l) figures into report/images"
