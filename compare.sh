#!/usr/bin/env bash
set -euo pipefail

DATA_DIR=./test/datasets
BUILD_DIR=./build
SKIP_BUILD=false
SELECTED_ALGOS=()
ALGO_TIMEOUT_SECONDS=60
PARALLEL_JOBS=1
CURRENT_CHILD_PID=""
DEFAULT_THRESHOLDS=(0 1 2 3 4 5 6)
MISSING_EDGE_TEST_GROUPS=(
    "missing_edge_threshold_0:0"
    "missing_edge_threshold_1:1"
    "missing_edge_threshold_2:2"
    "missing_edge_threshold_3:3"
    "missing_edge_threshold_4:4"
    "missing_edge_threshold_5:5"
    "missing_edge_threshold_6:6"
)

usage() {
    cat <<'EOF'
Usage: ./compare.sh [options]

Options:
  -a, --algorithms a,b,c   Only compare the specified algorithm keys.
                           Example: --algorithms cde_match,treespan
  -b, --build-dir DIR      Build directory. Default: ./build
  -d, --data-dir DIR       Dataset root. Default: ./test/datasets
                           The script scans synthetic/ and real_graphs/ under this root.
  -p, --parallel N         Run up to N algorithm tasks concurrently. Default: 1.
      --skip-build         Reuse existing binaries in the build directory.
  -h, --help              Show this help message.
EOF
}

on_interrupt() {
    trap - INT TERM
    echo
    echo "Interrupted by user (Ctrl+C). Stopping compare.sh."
    stop_current_child
    exit 130
}

stop_current_child() {
    local pids=()
    local pid

    if [ -n "${CURRENT_CHILD_PID:-}" ]; then
        pids+=("$CURRENT_CHILD_PID")
    fi

    while IFS= read -r pid; do
        if [ -n "$pid" ]; then
            pids+=("$pid")
        fi
    done < <(jobs -pr 2>/dev/null || true)

    if [ ${#pids[@]} -eq 0 ]; then
        return 0
    fi

    for pid in "${pids[@]}"; do
        kill -TERM -- "-$pid" 2>/dev/null || true
        pkill -TERM -P "$pid" 2>/dev/null || true
        kill -TERM "$pid" 2>/dev/null || true
    done
    sleep 1
    for pid in "${pids[@]}"; do
        kill -KILL -- "-$pid" 2>/dev/null || true
        pkill -KILL -P "$pid" 2>/dev/null || true
        kill -KILL "$pid" 2>/dev/null || true
    done
    CURRENT_CHILD_PID=""
}

run_with_timeout() {
    local output_file="$1"
    shift

    timeout --kill-after=5s "${ALGO_TIMEOUT_SECONDS}s" "$@" > "$output_file" 2>&1 &
    CURRENT_CHILD_PID=$!
    wait "$CURRENT_CHILD_PID"
    local rc=$?
    CURRENT_CHILD_PID=""
    return "$rc"
}

trap on_interrupt INT TERM

while [[ $# -gt 0 ]]; do
    case "$1" in
        -a|--algorithms)
            IFS=',' read -r -a SELECTED_ALGOS <<< "$2"
            shift 2
            ;;
        -b|--build-dir)
            BUILD_DIR="$2"
            shift 2
            ;;
        -d|--data-dir)
            DATA_DIR="$2"
            shift 2
            ;;
        -p|--parallel)
            if [ $# -lt 2 ] || ! [[ "$2" =~ ^[1-9][0-9]*$ ]]; then
                echo "Error: --parallel requires a positive integer" >&2
                exit 1
            fi
            PARALLEL_JOBS="$2"
            shift 2
            ;;
        --skip-build)
            SKIP_BUILD=true
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown argument: $1" >&2
            usage
            exit 1
            ;;
    esac
done

timestamp=$(date +"%Y%m%d_%H%M%S")
RESULT_DIR="./result/${timestamp}"
mkdir -p "${RESULT_DIR}"

OUT_LOG="${RESULT_DIR}/out.log"
echo "out.log will be saved to: ${PWD}/${OUT_LOG#./}"
exec > "$OUT_LOG" 2>&1

echo "===== RESULTS WILL BE SAVED TO: ${RESULT_DIR} ====="

if [ "$SKIP_BUILD" = false ]; then
    echo "===== Building algorithm executables ====="
    ./build.sh "$BUILD_DIR"
fi

discover_algorithms() {
    local build_dir="$1"
    shopt -s nullglob
    local candidates=("${build_dir}"/ssm_*)
    local resolved=()

    for exe in "${candidates[@]}"; do
        if [ -x "$exe" ] && [ ! -d "$exe" ]; then
            resolved+=("$exe")
        fi
    done

    printf '%s\n' "${resolved[@]}" | sort
}

normalize_algo_key() {
    local name
    name="$(basename "$1")"
    name="${name#ssm_}"
    name="${name%.exe}"
    printf '%s\n' "$name"
}

resolve_selected_algorithms() {
    local build_dir="$1"
    local discovered=()
    mapfile -t discovered < <(discover_algorithms "$build_dir")

    if [ ${#discovered[@]} -eq 0 ]; then
        echo "Error: no algorithm executables found under ${build_dir}" >&2
        exit 1
    fi

    if [ ${#SELECTED_ALGOS[@]} -eq 0 ]; then
        printf '%s\n' "${discovered[@]}"
        return
    fi

    local resolved=()
    local raw key candidate found
    for raw in "${SELECTED_ALGOS[@]}"; do
        key="${raw#ssm_}"
        key="${key%.exe}"
        found=false
        for candidate in "${build_dir}/ssm_${key}" "${build_dir}/ssm_${key}.exe"; do
            if [ -x "$candidate" ] && [ ! -d "$candidate" ]; then
                resolved+=("$candidate")
                found=true
                break
            fi
        done

        if [ "$found" = false ]; then
            echo "Error: algorithm '${raw}' was not found in ${build_dir}" >&2
            exit 1
        fi
    done

    printf '%s\n' "${resolved[@]}"
}

discover_queries() {
    local query_dir="$1"
    shopt -s nullglob
    local query_files=("${query_dir}"/*.txt)

    if [ ${#query_files[@]} -eq 0 ]; then
        return 0
    fi

    printf '%s\n' "${query_files[@]}" | \
        awk '
            function basename(path, parts, n) {
                n = split(path, parts, "/")
                return parts[n]
            }
            {
                file = basename($0)
                name = file
                sub(/\.txt$/, "", name)

                size_key = 2147483647
                index_key = 2147483647

                if (match(name, /^[0-9]+$/)) {
                    index_key = name + 0
                } else if (match(name, /^query_[0-9]+_[0-9]+$/)) {
                    split(name, fields, "_")
                    size_key = fields[2] + 0
                    index_key = fields[3] + 0
                }

                printf "%d\t%d\t%s\t%s\n", size_key, index_key, name, $0
            }
        ' | sort -t $'\t' -k1,1n -k2,2n -k3,3 | cut -f4-
}

discover_query_groups() {
    local qdir="$1"
    shopt -s nullglob

    local groups=()
    local direct_queries=("${qdir}"/*.txt)
    if [ ${#direct_queries[@]} -gt 0 ]; then
        groups+=(".	${qdir}")
    fi

    local group_dir group_name group_queries
    for group_dir in "${qdir}"/*/; do
        group_queries=("${group_dir}"/*.txt)
        if [ ${#group_queries[@]} -eq 0 ]; then
            continue
        fi

        group_name="$(basename "$group_dir")"
        groups+=("${group_name}	${group_dir%/}")
    done

    if [ ${#groups[@]} -eq 0 ]; then
        return 0
    fi

    printf '%s\n' "${groups[@]}" | \
        awk -F '\t' '
            {
                name = $1
                rank = 100
                value = 2147483647

                if (name == "." || name == "baseline") {
                    rank = 0
                    value = 0
                } else if (match(name, /^vertices_num_[0-9]+$/)) {
                    rank = 10
                    value = substr(name, 14) + 0
                } else if (match(name, /^avg_degree_[0-9]+$/)) {
                    rank = 20
                    value = substr(name, 12) + 0
                } else if (match(name, /^missing_edge_threshold_[0-9]+$/)) {
                    rank = 30
                    value = substr(name, 24) + 0
                }

                printf "%d\t%d\t%s\t%s\n", rank, value, name, $0
            }
        ' | sort -t $'\t' -k1,1n -k2,2n -k3,3 | cut -f4-
}

discover_dataset_dirs() {
    local data_dir="$1"
    local roots=()

    if [ -d "${data_dir}/synthetic" ]; then
        roots+=("synthetic	${data_dir}/synthetic")
    fi
    if [ -d "${data_dir}/real_graphs" ]; then
        roots+=("real	${data_dir}/real_graphs")
    fi

    if [ ${#roots[@]} -eq 0 ]; then
        roots+=("custom	${data_dir}")
    fi

    local root_entry dataset_group search_root
    for root_entry in "${roots[@]}"; do
        dataset_group="${root_entry%%$'\t'*}"
        search_root="${root_entry#*$'\t'}"

        find "$search_root" -name graph_g.txt -not -path '*/query_graph/*' -exec dirname {} \; | \
            awk -v group="$dataset_group" '{ printf "%s\t%s\n", group, $0 }'
    done | sort -t $'\t' -k1,1 -k2,2
}

thresholds_for_query_group() {
    local query_group="$1"
    local entry group threshold

    for entry in "${MISSING_EDGE_TEST_GROUPS[@]}"; do
        group="${entry%%:*}"
        threshold="${entry#*:}"

        if [ "$query_group" = "$group" ]; then
            printf '%s\n' "$threshold"
            return 0
        fi
    done

    printf '%s\n' "${DEFAULT_THRESHOLDS[@]}"
}

parse_summary_line() {
    local line="$1"
    SUMMARY_ALGO=""
    SUMMARY_COUNT=""
    SUMMARY_LOAD_MS=""
    SUMMARY_RUN_MS=""
    SUMMARY_TOTAL_MS=""

    for token in $line; do
        case "$token" in
            algorithm=*) SUMMARY_ALGO="${token#algorithm=}" ;;
            count=*) SUMMARY_COUNT="${token#count=}" ;;
            load_ms=*) SUMMARY_LOAD_MS="${token#load_ms=}" ;;
            run_ms=*) SUMMARY_RUN_MS="${token#run_ms=}" ;;
            total_ms=*) SUMMARY_TOTAL_MS="${token#total_ms=}" ;;
        esac
    done

    if [ -z "$SUMMARY_ALGO" ] || [ -z "$SUMMARY_COUNT" ] || [ -z "$SUMMARY_RUN_MS" ] || [ -z "$SUMMARY_TOTAL_MS" ]; then
        return 1
    fi
}

parse_output_metadata() {
    local output_file="$1"
    local fallback_label="$2"

    OUTPUT_LABEL="$fallback_label"
    OUTPUT_RECURSION_CALLS="NA"

    if [ ! -f "$output_file" ]; then
        return 0
    fi

    local parsed_label parsed_recursion
    parsed_label=$(sed -n 's/^\(.*\) Results:$/\1/p' "$output_file" | tail -n 1 || true)
    if [ -n "$parsed_label" ]; then
        OUTPUT_LABEL="$parsed_label"
    fi

    parsed_recursion=$(sed -n 's/^Recursion Calls:[[:space:]]*//p' "$output_file" | tail -n 1 || true)
    parsed_recursion="${parsed_recursion//[[:space:]]/}"
    if [ -n "$parsed_recursion" ]; then
        OUTPUT_RECURSION_CALLS="$parsed_recursion"
    fi
}

sanitize_tsv_field() {
    local value="$1"
    value="${value//$'\t'/ }"
    value="${value//$'\r'/ }"
    value="${value//$'\n'/ }"
    printf '%s' "$value"
}

task_status_file() {
    local case_id="$1"
    local algo="$2"
    printf '%s/%s_%s.status\n' "$STATUS_DIR" "$case_id" "$algo"
}

write_task_status() {
    local status_file="$1"
    shift

    local fields=()
    local value
    for value in "$@"; do
        fields+=("$(sanitize_tsv_field "$value")")
    done

    printf '%s\n' "$(IFS=$'\t'; echo "${fields[*]}")" > "$status_file"
}

run_task() {
    local case_id="$1"
    local algo_index="$2"
    local exe="$3"
    local algo="$4"
    local gfile="$5"
    local qfile="$6"
    local threshold="$7"
    local out="$8"
    local abs_out="$9"

    local status_file rc status summary_line count run_value recursion label
    status_file="$(task_status_file "$case_id" "$algo")"

    set +e
    run_with_timeout "$out" "$exe" -d "$gfile" -q "$qfile" -t "$threshold"
    rc=$?
    set -e

    parse_output_metadata "$out" "$algo"
    label="$OUTPUT_LABEL"
    recursion="$OUTPUT_RECURSION_CALLS"
    count="NA"
    run_value="NA"

    if [ "$rc" -ne 0 ]; then
        if [ "$rc" -eq 124 ] || [ "$rc" -eq 137 ]; then
            status="Timeout"
        else
            status="RunError"
        fi
        write_task_status "$status_file" "$status" "$rc" "$count" "$run_value" "$recursion" "$abs_out" "$label"
        return 0
    fi

    summary_line=$(grep '^SSM_SUMMARY ' "$out" | tail -n 1 || true)
    if parse_summary_line "$summary_line"; then
        count="$SUMMARY_COUNT"
        run_value="$SUMMARY_RUN_MS"
        if grep -Eq '(Output|Intermediate) Limit:.*\(reached\)' "$out"; then
            status="ResourceLimit"
        else
            status="OK"
        fi
    else
        status="ParseError"
    fi

    write_task_status "$status_file" "$status" "$rc" "$count" "$run_value" "$recursion" "$abs_out" "$label"
    return 0
}

run_all_tasks() {
    local launched=0
    local case_id algo_index exe algo gfile qfile threshold out abs_out
    local sem_fifo="${RESULT_DIR}/parallel.sem"
    local token i

    if [ "$TOTAL_ALGO_RUNS" -eq 0 ]; then
        return 0
    fi

    echo "Running algorithm tasks with parallel=${PARALLEL_JOBS}"

    rm -f "$sem_fifo"
    mkfifo "$sem_fifo"
    exec 9<>"$sem_fifo"
    rm -f "$sem_fifo"

    for ((i = 0; i < PARALLEL_JOBS; i++)); do
        printf '.' >&9
    done

    while IFS=$'\t' read -r case_id algo_index exe algo gfile qfile threshold out abs_out; do
        launched=$((launched + 1))
        IFS= read -r -n 1 token <&9
        (
            trap 'printf "." >&9' EXIT
            run_task "$case_id" "$algo_index" "$exe" "$algo" "$gfile" "$qfile" "$threshold" "$out" "$abs_out"
        ) &

        if [ $((launched % PARALLEL_JOBS)) -eq 0 ] || [ "$launched" -eq "$TOTAL_ALGO_RUNS" ]; then
            echo "  Launched ${launched}/${TOTAL_ALGO_RUNS} algorithm tasks"
        fi
    done < "$TASK_MANIFEST"

    wait || true
    exec 9>&-
    exec 9<&-
}

record_failure() {
    local dataset_group="$1"
    local dirname="$2"
    local query_group_label="$3"
    local qname="$4"
    local threshold="$5"
    local algo="$6"
    local status="$7"
    local output="$8"

    case "$status" in
        Timeout)
            printf "%s\t%s\t%s\t%s\tt=%s\t%s\tTimeout(%ss)\t%s\n" \
                "$dataset_group" "$dirname" "$query_group_label" "$qname" "$threshold" "$algo" "$ALGO_TIMEOUT_SECONDS" "$output" >> "$FAIL_LIST"
            ;;
        RunError|ParseError|MissingStatus)
            printf "%s\t%s\t%s\t%s\tt=%s\t%s\t%s\t%s\n" \
                "$dataset_group" "$dirname" "$query_group_label" "$qname" "$threshold" "$algo" "$status" "$output" >> "$FAIL_LIST"
            ;;
        *)
            printf "%s\t%s\t%s\t%s\tt=%s\t%s\t%s\t%s\n" \
                "$dataset_group" "$dirname" "$query_group_label" "$qname" "$threshold" "$algo" "$status" "$output" >> "$FAIL_LIST"
            ;;
    esac
}

aggregate_results() {
    local case_id dataset_group dirname query_group_label qname t

    while IFS=$'\t' read -r case_id dataset_group dirname query_group_label qname t; do
        declare -A counts=()
        declare -A run_ms=()
        declare -A outputs=()
        declare -A statuses=()
        declare -A display_names=()
        declare -A recursion_calls=()

        local case_failed=false
        local case_reason=""
        local idx algo status_file status rc count run_value recursion output label

        for idx in "${!ALGO_KEYS[@]}"; do
            algo="${ALGO_KEYS[$idx]}"
            status_file="$(task_status_file "$case_id" "$algo")"

            if [ ! -f "$status_file" ]; then
                status="MissingStatus"
                rc="NA"
                count="NA"
                run_value="NA"
                recursion="NA"
                output="NA"
                label="$algo"
            else
                IFS=$'\t' read -r status rc count run_value recursion output label < "$status_file"
            fi

            counts["$algo"]="$count"
            run_ms["$algo"]="$run_value"
            recursion_calls["$algo"]="$recursion"
            outputs["$algo"]="$output"
            statuses["$algo"]="$status"
            display_names["$algo"]="${label:-$algo}"

            if [ "$status" != "OK" ]; then
                case_failed=true
                record_failure "$dataset_group" "$dirname" "$query_group_label" "$qname" "$t" "$algo" "$status" "$output"
            fi
        done

        local expected_count=""
        if [ "$case_failed" = false ]; then
            local baseline_algo="${ALGO_KEYS[0]}"
            expected_count="${counts[$baseline_algo]}"
            for algo in "${ALGO_KEYS[@]}"; do
                if [ "${counts[$algo]}" != "$expected_count" ]; then
                    case_failed=true
                    case_reason="CountMismatch"
                    record_failure "$dataset_group" "$dirname" "$query_group_label" \
                        "$qname" "$t" "$algo" \
                        "CountMismatch(expected=${expected_count},actual=${counts[$algo]})" \
                        "${outputs[$algo]}"
                fi
            done
        fi

        local row=("$dataset_group" "$dirname" "$query_group_label" "$qname" "$t")
        if [ "$case_failed" = true ]; then
            FAIL=$((FAIL + 1))
            if [ -z "$case_reason" ]; then
                case_reason="RuntimeOrParseError"
            fi
            echo "      [FAIL] ${dataset_group}/${dirname} ${query_group_label}/${qname} t=$t"
            echo "             Reason: $case_reason"
            row+=("FAIL" "$expected_count")
        else
            echo "      [DONE] ${dataset_group}/${dirname} ${query_group_label}/${qname} t=$t"
            row+=("OK" "$expected_count")
        fi

        for algo in "${ALGO_KEYS[@]}"; do
            local algo_label="${display_names[$algo]:-$algo}"
            local algo_status="${statuses[$algo]:-NA}"
            local algo_output="${outputs[$algo]:-NA}"
            local algo_run="${run_ms[$algo]:-NA}"
            local algo_recursion="${recursion_calls[$algo]:-NA}"

            if [ "$algo_status" = "OK" ]; then
                echo "             ${algo_label}: ${algo_run} ms, Recursion Calls: ${algo_recursion}, Output: ${algo_output}"
            elif [ "$algo_status" = "Timeout" ]; then
                echo "             ${algo_label}: TIMEOUT(${ALGO_TIMEOUT_SECONDS}s), Recursion Calls: ${algo_recursion}, Output: ${algo_output}"
            else
                echo "             ${algo_label}: ${algo_status}, Recursion Calls: ${algo_recursion}, Output: ${algo_output}"
            fi
        done

        for algo in "${ALGO_KEYS[@]}"; do
            row+=("${counts[$algo]:-NA}" "${run_ms[$algo]:-NA}" "${recursion_calls[$algo]:-NA}" "${outputs[$algo]:-NA}")
        done
        printf '%s\n' "$(IFS=$'\t'; echo "${row[*]}")" >> "$SUMMARY_FILE"
    done < "$CASE_MANIFEST"
}

mapfile -t ALGO_EXECUTABLES < <(resolve_selected_algorithms "$BUILD_DIR")
ALGO_KEYS=()
for exe in "${ALGO_EXECUTABLES[@]}"; do
    ALGO_KEYS+=("$(normalize_algo_key "$exe")")
done

echo "Algorithms under comparison: ${ALGO_KEYS[*]}"

SUMMARY_FILE="${RESULT_DIR}/summary.tsv"
FAIL_LIST="${RESULT_DIR}/failures.txt"
: > "$FAIL_LIST"

summary_header=(dataset_group dataset query_group query threshold status expected_count)
for algo in "${ALGO_KEYS[@]}"; do
    summary_header+=("${algo}_count" "${algo}_run_ms" "${algo}_recursion_calls" "${algo}_output")
done
printf '%s\n' "$(IFS=$'\t'; echo "${summary_header[*]}")" > "$SUMMARY_FILE"

TOTAL_CASES=0
TOTAL_ALGO_RUNS=0
FAIL=0

TASK_MANIFEST="${RESULT_DIR}/tasks.tsv"
CASE_MANIFEST="${RESULT_DIR}/cases.tsv"
STATUS_DIR="${RESULT_DIR}/task_status"
mkdir -p "$STATUS_DIR"
: > "$TASK_MANIFEST"
: > "$CASE_MANIFEST"

shopt -s nullglob

mapfile -t DATASET_DIRS < <(discover_dataset_dirs "$DATA_DIR")

if [ ${#DATASET_DIRS[@]} -eq 0 ]; then
    echo "Error: no graph_g.txt datasets found under ${DATA_DIR}/synthetic or ${DATA_DIR}/real_graphs" >&2
    exit 1
fi

for dataset_entry in "${DATASET_DIRS[@]}"; do
    dataset_group="${dataset_entry%%$'\t'*}"
    folder="${dataset_entry#*$'\t'}"
    dirname=$(basename "$folder")
    gfile="${folder}/graph_g.txt"
    qdir="${folder}/query_graph"

    if [ ! -f "$gfile" ]; then
        echo "[Skip] ${dataset_group}/${dirname}: graph_g.txt not found"
        continue
    fi

    if [ ! -d "$qdir" ]; then
        echo "[Skip] ${dataset_group}/${dirname}: query_graph/ not found"
        continue
    fi

    echo "Dataset: ${dataset_group}/${dirname}"
    mkdir -p "${RESULT_DIR}/${dataset_group}/${dirname}"

    mapfile -t QUERY_GROUPS < <(discover_query_groups "$qdir")
    if [ ${#QUERY_GROUPS[@]} -eq 0 ]; then
        echo "[Skip] ${dataset_group}/${dirname}: no query .txt files found under query_graph/"
        continue
    fi

    for query_group_entry in "${QUERY_GROUPS[@]}"; do
        query_group="${query_group_entry%%$'\t'*}"
        query_group_dir="${query_group_entry#*$'\t'}"
        query_group_label="$query_group"
        if [ "$query_group_label" = "." ]; then
            query_group_label="legacy"
        fi

        echo "  Query group: $query_group_label"

        mapfile -t QUERY_FILES < <(discover_queries "$query_group_dir")
        mapfile -t THRESHOLDS < <(thresholds_for_query_group "$query_group")
        for qfile in "${QUERY_FILES[@]}"; do
            qname=$(basename "$qfile" .txt)
            echo "    Query: $qname"

            for t in "${THRESHOLDS[@]}"; do
                TOTAL_CASES=$((TOTAL_CASES + 1))
                case_id=$(printf 'case_%06d' "$TOTAL_CASES")
                printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
                    "$case_id" "$dataset_group" "$dirname" "$query_group_label" "$qname" "$t" >> "$CASE_MANIFEST"

                for idx in "${!ALGO_EXECUTABLES[@]}"; do
                    exe="${ALGO_EXECUTABLES[$idx]}"
                    algo="${ALGO_KEYS[$idx]}"
                    algo_result_dir="${RESULT_DIR}/${dataset_group}/${dirname}/${query_group_label}/${algo}"
                    mkdir -p "${algo_result_dir}"
                    out="${algo_result_dir}/${qname}_t=${t}.txt"
                    abs_out="${PWD}/${out#./}"

                    TOTAL_ALGO_RUNS=$((TOTAL_ALGO_RUNS + 1))
                    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
                        "$case_id" "$idx" "$exe" "$algo" "$gfile" "$qfile" "$t" "$out" "$abs_out" >> "$TASK_MANIFEST"
                done
            done
        done
    done
done

echo "Queued comparison cases: $TOTAL_CASES"
echo "Queued algorithm tasks:  $TOTAL_ALGO_RUNS"

run_all_tasks

echo "Aggregating results..."
aggregate_results

echo "===== ALL DONE ====="
echo "Results saved at: $(realpath "${RESULT_DIR}")"
echo "Summary TSV: $(realpath "${SUMMARY_FILE}")"
echo "Algorithms: ${ALGO_KEYS[*]}"
echo "Comparison cases: $TOTAL_CASES"
echo "Algorithm runs:   $TOTAL_ALGO_RUNS"
echo "Parallel jobs:    $PARALLEL_JOBS"
echo "Failed cases:     $FAIL"

if [ "$FAIL" -ne 0 ]; then
    echo "Failure list: $(realpath "$FAIL_LIST")"
    exit 2
fi
