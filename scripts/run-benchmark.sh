#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
  cat <<'USAGE'
Usage:
  scripts/run-benchmark.sh <benchmark>

Available benchmarks:
  02_star_constant
  03_star_linear
  04_star_sine
  05_star_radial_decay
  06_star_localized
USAGE
  exit 1
fi

benchmark="$1"

configure_and_build() {
  local target="$1"
  cmake -S . -B build
  cmake --build build --target "$target"
}

copy_vtp_series() {
  local source_dir="$1"
  local destination_dir="$2"

  mkdir -p "$destination_dir"
  cp "$source_dir/solution.pvd" "$destination_dir/solution.pvd"

  grep -o 'solution_[0-9]\{4\}\.vtp' "$source_dir/solution.pvd" \
    | sort -u \
    | while read -r vtp_file; do
        cp "$source_dir/$vtp_file" "$destination_dir/$vtp_file"
      done
}

try_plot_space_convergence() {
  local csv_path="$1"
  if python3 scripts/plot-convergence.py "$csv_path"; then
    cp output/convergence/star_sine_space_convergence.pdf \
      benchmarks/04_star_sine/results/space_convergence.pdf
  else
    echo "Warning: unable to generate space convergence PDF. CSV was still copied." >&2
  fi
}

try_plot_time_convergence() {
  local csv_path="$1"
  local pdf_path="$2"
  if ! python3 scripts/plot-time-convergence.py "$csv_path" "$pdf_path"; then
    echo "Warning: unable to generate time convergence PDF. CSV was still copied." >&2
  fi
}

case "$benchmark" in
  02_star_constant)
    configure_and_build test_heat_star_constant
    ./build/test_heat_star_constant
    copy_vtp_series \
      output/visualization/star_constant \
      benchmarks/02_star_constant/results
    ;;

  03_star_linear)
    configure_and_build test_heat_star_linear
    ./build/test_heat_star_linear
    copy_vtp_series \
      output/visualization/star_linear \
      benchmarks/03_star_linear/results
    ;;

  04_star_sine)
    configure_and_build test_heat_star_sine
    cmake --build build --target test_heat_star_sine_convergence
    ./build/test_heat_star_sine
    ./build/test_heat_star_sine_convergence
    copy_vtp_series \
      output/visualization/star_sine \
      benchmarks/04_star_sine/results
    cp output/convergence/star_sine_space_convergence.csv \
      benchmarks/04_star_sine/results/space_convergence.csv
    try_plot_space_convergence output/convergence/star_sine_space_convergence.csv
    ;;

  05_star_radial_decay)
    configure_and_build test_heat_star_radial_decay
    cmake --build build --target test_heat_star_radial_decay_time_convergence
    ./build/test_heat_star_radial_decay
    ./build/test_heat_star_radial_decay_time_convergence
    copy_vtp_series \
      output/visualization/star_radial_decay \
      benchmarks/05_star_radial_decay/results
    cp output/convergence/star_radial_decay_time_convergence.csv \
      benchmarks/05_star_radial_decay/results/time_convergence.csv
    try_plot_time_convergence \
      output/convergence/star_radial_decay_time_convergence.csv \
      benchmarks/05_star_radial_decay/results/time_convergence.pdf
    ;;

  06_star_localized)
    configure_and_build test_heat_star_localized
    ./build/test_heat_star_localized
    copy_vtp_series \
      output/visualization/star_localized \
      benchmarks/06_star_localized/results
    ;;

  *)
    echo "Unknown benchmark: $benchmark" >&2
    echo "Run scripts/run-benchmark.sh without arguments to list available benchmarks." >&2
    exit 1
    ;;
esac

echo "Benchmark completed: $benchmark"
