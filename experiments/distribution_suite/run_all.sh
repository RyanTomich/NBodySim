#!/usr/bin/env bash
# Run the three production solvers through the shared C benchmark driver.
set -euo pipefail

suite_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
repo_dir=$(cd -- "$suite_dir/../.." && pwd -P)

usage() {
    cat <<'HELP'
Usage: run_all.sh [--list] [INPUT.txt ...]

With no inputs, run every inputs/*/*.txt file beside this script.
The same benchmark executable selects n2, barnes-hut and fmm internally.

Environment settings:
  BENCH_BIN   Path to benchmark binary (default: NBodySim/build/bench)
  THREADS     Space-separated thread budgets (default: "1")
  DT          Timestep (default: 1e-5)
  STEPS       Leapfrog steps per simulation (default: 10)
  SAMPLES     Timing samples per configuration (default: 30, minimum 2)
  WARMUPS     Untimed simulation warmups (default: 3)
  THETA       Barnes-Hut opening threshold (default: 0.3)
  ORDER       FMM expansion order (default: 4)
  OUTPUT_DIR  New output directory (default: results/distribution-suite/<UTC>-<pid>)

Examples, from NBodySim:
  experiments/distribution_suite/run_all.sh --list
  SAMPLES=3 STEPS=2 WARMUPS=1 experiments/distribution_suite/run_all.sh \
    experiments/distribution_suite/inputs/pbbs_cube/n256.txt
  THREADS="1 2 4 8" experiments/distribution_suite/run_all.sh

All three use unsoftened G=1 and the same saved input, dt and duration.
OUTPUT_DIR must not exist. No builds or input generation are performed here.
HELP
}

fail() { printf 'Error: %s\n' "$*" >&2; exit 2; }

list_only=0
case "${1:-}" in
    -h|--help) usage; exit 0 ;;
    --list) list_only=1; shift ;;
    --) shift ;;
esac
shopt -s nullglob
if (( $# )); then
    inputs=("$@")
else
    inputs=("$suite_dir"/inputs/*/*.txt)
fi
(( ${#inputs[@]} )) || fail "No inputs found. Run python3 $suite_dir/prepare.py first."
for input in "${inputs[@]}"; do
    [[ -f "$input" ]] || fail "Input does not exist: $input"
    [[ "$input" != *$'\t'* && "$input" != *$'\n'* ]] || fail "Input paths must not contain tabs or newlines."
    if (( list_only )); then printf '%s\n' "$input"; fi
done
if (( list_only )); then exit 0; fi

bench_bin=${BENCH_BIN:-"$repo_dir/build/bench"}
dt=${DT:-1e-5}
steps=${STEPS:-10}
samples=${SAMPLES:-30}
warmups=${WARMUPS:-3}
theta=${THETA:-0.3}
order=${ORDER:-4}
[[ -x "$bench_bin" ]] || fail "Benchmark binary is missing: $bench_bin. Build with make, or set BENCH_BIN."
[[ "$steps" =~ ^[1-9][0-9]*$ ]] || fail 'STEPS must be a positive integer.'
[[ "$samples" =~ ^[1-9][0-9]*$ ]] && (( samples >= 2 )) || fail 'SAMPLES must be an integer >= 2.'
[[ "$warmups" =~ ^(0|[1-9][0-9]*)$ ]] || fail 'WARMUPS must be a nonnegative integer.'
[[ "$order" =~ ^[1-6]$ ]] || fail 'ORDER must be an integer from 1 to 6.'
read -r -a thread_counts <<< "${THREADS:-1}"
(( ${#thread_counts[@]} )) || fail 'THREADS is empty.'
seen_threads=' '
for threads in "${thread_counts[@]}"; do
    [[ "$threads" =~ ^[1-9][0-9]*$ ]] || fail 'THREADS must contain positive integers.'
    [[ "$seen_threads" != *" $threads "* ]] || fail "Repeated thread count: $threads"
    seen_threads="$seen_threads$threads "
done
build_info=$("$bench_bin" --build-info)
[[ "$build_info" =~ \"instrumented\"[[:space:]]*:[[:space:]]*false ]] || fail 'Use an uninstrumented benchmark binary.'
for threads in "${thread_counts[@]}"; do
    if (( threads > 1 )); then
        [[ "$build_info" =~ \"openmp\"[[:space:]]*:[[:space:]]*true ]] || fail 'Multiple threads require an OPENMP=1 build.'
    fi
done

output_dir=${OUTPUT_DIR:-"$repo_dir/results/distribution-suite/$(date -u +%Y%m%dT%H%M%SZ)-$$"}
mkdir -p -- "$(dirname -- "$output_dir")"
mkdir -- "$output_dir" || fail 'OUTPUT_DIR must be a new directory.'
output_dir=$(cd -- "$output_dir" && pwd -P)
active_folder=''
trap 'result=$?; if (( result == 0 )); then printf "ok\n"; else printf "failed\n"; if [[ -n "$active_folder" ]]; then printf "failed\n" > "$active_folder/status.txt"; fi; fi > "$output_dir/status.txt"' EXIT
printf 'running\n' > "$output_dir/status.txt"
cp -- "$bench_bin" "$output_dir/bench"
bench_bin="$output_dir/bench"
printf '%s\n' "$build_info" > "$output_dir/build-info.json"
printf 'dt=%s\nsteps=%s\nsamples=%s\nwarmups=%s\ntheta=%s\norder=%s\nthreads=%s\nphysics=unsoftened G=1\n' \
    "$dt" "$steps" "$samples" "$warmups" "$theta" "$order" "${thread_counts[*]}" > "$output_dir/settings.txt"
printf 'input\tbodies\tsolver\tthreads\tdt\tsteps\tsamples\ttheta\torder\tmean_ms_per_step\tstddev_ms_per_step\tspeedup_vs_n2_same_threads\tposition_rms_vs_n2\tvelocity_rms_vs_n2\tstatus\toutput\n' > "$output_dir/summary.tsv"
printf 'Results: %s\n' "$output_dir"

input_number=0
for input in "${inputs[@]}"; do
    input_number=$((input_number + 1))
    input=$(cd -- "$(dirname -- "$input")" && printf '%s/%s' "$PWD" "$(basename -- "$input")")
    family=$(basename -- "$(dirname -- "$input")")
    stem=$(basename -- "$input" .txt)
    label=${family//[^a-zA-Z0-9_.-]/_}-${stem//[^a-zA-Z0-9_.-]/_}
    printf -v case_name '%04d-%s' "$input_number" "$label"
    case_dir="$output_dir/$case_name"
    mkdir -- "$case_dir"
    cp -- "$input" "$case_dir/input.txt"
    if [[ -f "${input%.txt}.source.json" ]]; then cp -- "${input%.txt}.source.json" "$case_dir/source.json"; fi
    read -r bodies < "$case_dir/input.txt"
    for threads in "${thread_counts[@]}"; do
        baseline_ms=''
        reference_file="$case_dir/n2-${threads}threads/final.txt"
        for solver in n2 barnes-hut fmm; do
            folder="$case_dir/$solver-${threads}threads"
            mkdir -- "$folder"
            active_folder=$folder
            args=("$bench_bin" "$case_dir/input.txt" "$dt" "$steps" "$samples" "$warmups" timing "$folder/final.txt" "$solver")
            selected_theta=''; selected_order=''
            case "$solver" in
                barnes-hut) args+=(theta "$theta"); selected_theta=$theta ;;
                fmm) args+=(order "$order"); selected_order=$order ;;
            esac
            printf '[%d/%d] %s: %s, %s thread(s)\n' "$input_number" "${#inputs[@]}" "$label" "$solver" "$threads"
            if OMP_NUM_THREADS="$threads" OMP_THREAD_LIMIT="$threads" OMP_DYNAMIC=FALSE OMP_MAX_ACTIVE_LEVELS=1 \
                    "${args[@]}" > "$folder/samples.csv" 2> "$folder/stderr.log"; then
                :
            else
                result=$?
                printf 'failed\n' > "$folder/status.txt"
                printf 'Solver failed (exit %s). See %s\n' "$result" "$folder/stderr.log" >&2
                exit "$result"
            fi
            # Stable sample variance of whole-run seconds; normalize by steps.
            timing=$(awk -F, -v expected="$samples" -v steps="$steps" '
                NR == 1 { if ($1 != "sample" || $2 != "seconds") bad=1; next }
                { if ($1 != n+1 || $2 !~ /^[0-9]+([.][0-9]*)?([eE][-+]?[0-9]+)?$/) bad=1;
                  n++; delta=$2-mean; mean+=delta/n; m2+=delta*($2-mean) }
                END { if (bad || n != expected || mean <= 0) exit 1;
                      printf "%.17g\t%.17g", 1000*mean/steps, 1000*sqrt(m2/(n-1))/steps }
            ' "$folder/samples.csv") || fail "Invalid timing samples: $folder/samples.csv"
            IFS=$'\t' read -r mean_ms stddev_ms <<< "$timing"
            if [[ "$solver" == n2 ]]; then baseline_ms=$mean_ms; fi
            speedup=$(awk -v b="$baseline_ms" -v t="$mean_ms" 'BEGIN { printf "%.17g", b/t }')
            errors=$(awk '
                NR == FNR { if (FNR == 1) { count=$1; next }
                            for (k=2; k<=7; k++) ref[FNR,k]=$k; next }
                FNR == 1 { if ($1 != count) bad=1; next }
                { n++; if (NF != 10) bad=1;
                  for (k=2; k<=4; k++) pos+=($k-ref[FNR,k])^2;
                  for (k=5; k<=7; k++) vel+=($k-ref[FNR,k])^2 }
                END { if (bad || n != count || n == 0) exit 1;
                      printf "%.17g\t%.17g", sqrt(pos/n), sqrt(vel/n) }
            ' "$reference_file" "$folder/final.txt") || fail "Invalid final state: $folder/final.txt"
            IFS=$'\t' read -r position_error velocity_error <<< "$errors"
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tok\t%s\n' \
                "$input" "$bodies" "$solver" "$threads" "$dt" "$steps" "$samples" "$selected_theta" "$selected_order" \
                "$mean_ms" "$stddev_ms" "$speedup" "$position_error" "$velocity_error" "$folder" >> "$output_dir/summary.tsv"
            printf 'ok\n' > "$folder/status.txt"
            active_folder=''
        done
    done
done
printf 'Completed. Summary: %s/summary.tsv\n' "$output_dir"
