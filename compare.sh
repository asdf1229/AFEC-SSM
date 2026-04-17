#!/usr/bin/env bash
set -euo pipefail

DATA_DIR=../SSM-DatasetGen/datasets
BUILD_DIR=./build
SKIP_BUILD=false
SELECTED_ALGOS=()
ALGO_TIMEOUT_SECONDS=30
CURRENT_CHILD_PID=""

usage() {
    cat <<'EOF'
Usage: ./compare.sh [options]

Options:
  -a, --algorithms a,b,c   Only compare the specified algorithm keys.
                           Example: --algorithms cde_match,treespan
  -b, --build-dir DIR      Build directory. Default: ./build
  -d, --data-dir DIR       Dataset root. Default: ../SSM-DatasetGen/datasets
                           The script scans synthetic/ and raw/real_graphs/ under this root.
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
    local pid="${CURRENT_CHILD_PID:-}"
    if [ -z "$pid" ]; then
        return 0
    fi

    kill -TERM -- "-$pid" 2>/dev/null || true
    pkill -TERM -P "$pid" 2>/dev/null || true
    kill -TERM "$pid" 2>/dev/null || true
    sleep 1
    kill -KILL -- "-$pid" 2>/dev/null || true
    pkill -KILL -P "$pid" 2>/dev/null || true
    kill -KILL "$pid" 2>/dev/null || true
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
    local candidates=("${build_dir}"/ssm_ged_*)
    local resolved=()

    for exe in "${candidates[@]}"; do
        if [ -f "$exe" ]; then
            resolved+=("$exe")
        fi
    done

    printf '%s\n' "${resolved[@]}" | sort
}

normalize_algo_key() {
    local name
    name="$(basename "$1")"
    name="${name#ssm_ged_}"
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
        key="${raw#ssm_ged_}"
        key="${key%.exe}"
        found=false
        for candidate in "${build_dir}/ssm_ged_${key}" "${build_dir}/ssm_ged_${key}.exe"; do
            if [ -f "$candidate" ]; then
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
    if [ -d "${data_dir}/raw/real_graphs" ]; then
        roots+=("real	${data_dir}/raw/real_graphs")
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

shopt -s nullglob

mapfile -t DATASET_DIRS < <(discover_dataset_dirs "$DATA_DIR")

if [ ${#DATASET_DIRS[@]} -eq 0 ]; then
    echo "Error: no graph_g.txt datasets found under ${DATA_DIR}/synthetic or ${DATA_DIR}/raw/real_graphs" >&2
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
        for qfile in "${QUERY_FILES[@]}"; do
            qname=$(basename "$qfile" .txt)
            echo "    Query: $qname"

            for t in {0..3}; do
                TOTAL_CASES=$((TOTAL_CASES + 1))

                declare -A counts=()
                declare -A run_ms=()
                declare -A outputs=()
                declare -A statuses=()
                declare -A display_names=()
                declare -A recursion_calls=()

                case_failed=false
                case_reason=""

                for idx in "${!ALGO_EXECUTABLES[@]}"; do
                    exe="${ALGO_EXECUTABLES[$idx]}"
                    algo="${ALGO_KEYS[$idx]}"
                    algo_result_dir="${RESULT_DIR}/${dataset_group}/${dirname}/${query_group_label}/${algo}"
                    mkdir -p "${algo_result_dir}"
                    out="${algo_result_dir}/${qname}_t=${t}.txt"
                    abs_out="${PWD}/${out#./}"
                    outputs["$algo"]="$abs_out"
                    display_names["$algo"]="$algo"
                    recursion_calls["$algo"]="NA"

                    TOTAL_ALGO_RUNS=$((TOTAL_ALGO_RUNS + 1))

                    set +e
                    run_with_timeout "$out" "$exe" -d "$gfile" -q "$qfile" -t "$t"
                    rc=$?
                    set -e

                    parse_output_metadata "$out" "$algo"
                    display_names["$algo"]="$OUTPUT_LABEL"
                    recursion_calls["$algo"]="$OUTPUT_RECURSION_CALLS"

                    if [ $rc -ne 0 ]; then
                        if [ $rc -eq 130 ]; then
                            echo
                            echo "Interrupted by user (Ctrl+C). Stopping compare.sh."
                            exit 130
                        fi
                        case_failed=true
                        if [ $rc -eq 124 ] || [ $rc -eq 137 ]; then
                            statuses["$algo"]="Timeout"
                            printf "%s\t%s\t%s\t%s\tt=%s\t%s\tTimeout(%ss)\t%s\n" \
                                "$dataset_group" "$dirname" "$query_group_label" "$qname" "$t" "$algo" "$ALGO_TIMEOUT_SECONDS" "$abs_out" >> "$FAIL_LIST"
                        else
                            statuses["$algo"]="RunError"
                            printf "%s\t%s\t%s\t%s\tt=%s\t%s\tRunError\t%s\n" \
                                "$dataset_group" "$dirname" "$query_group_label" "$qname" "$t" "$algo" "$abs_out" >> "$FAIL_LIST"
                        fi
                        continue
                    fi

                    summary_line=$(grep '^SSM_GED_SUMMARY ' "$out" | tail -n 1 || true)
                    if ! parse_summary_line "$summary_line"; then
                        case_failed=true
                        statuses["$algo"]="ParseError"
                        printf "%s\t%s\t%s\t%s\tt=%s\t%s\tParseError\t%s\n" \
                            "$dataset_group" "$dirname" "$query_group_label" "$qname" "$t" "$algo" "$abs_out" >> "$FAIL_LIST"
                        continue
                    fi

                    counts["$algo"]="$SUMMARY_COUNT"
                    run_ms["$algo"]="$SUMMARY_RUN_MS"
                    statuses["$algo"]="OK"
                done

                expected_count=""
                if [ "$case_failed" = false ]; then
                    baseline_algo="${ALGO_KEYS[0]}"
                    expected_count="${counts[$baseline_algo]}"

                    for algo in "${ALGO_KEYS[@]}"; do
                        if [ "${counts[$algo]}" != "$expected_count" ]; then
                            case_failed=true
                            case_reason="CountMismatch ${baseline_algo}=${expected_count}, ${algo}=${counts[$algo]}"
                            printf "%s\t%s\t%s\t%s\tt=%s\t%s\tVerifyError\t%s\n" \
                                "$dataset_group" "$dirname" "$query_group_label" "$qname" "$t" "$algo" "${outputs[$algo]}" >> "$FAIL_LIST"
                        fi
                    done
                fi

                row=("$dataset_group" "$dirname" "$query_group_label" "$qname" "$t")
                if [ "$case_failed" = true ]; then
                    FAIL=$((FAIL + 1))
                    if [ -z "$case_reason" ]; then
                        case_reason="RuntimeOrParseError"
                    fi
                    echo "      [FAIL] t=$t"
                    echo "             Reason: $case_reason"
                    row+=("FAIL" "$expected_count")
                else
                    echo "      [DONE] t=$t"
                    row+=("OK" "$expected_count")
                fi

                for algo in "${ALGO_KEYS[@]}"; do
                    algo_label="${display_names[$algo]:-$algo}"
                    algo_status="${statuses[$algo]:-NA}"
                    algo_output="${outputs[$algo]:-NA}"
                    algo_run="${run_ms[$algo]:-NA}"
                    algo_recursion="${recursion_calls[$algo]:-NA}"

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
            done
        done
    done
done

echo "===== ALL DONE ====="
echo "Results saved at: $(realpath "${RESULT_DIR}")"
echo "Summary TSV: $(realpath "${SUMMARY_FILE}")"
echo "Algorithms: ${ALGO_KEYS[*]}"
echo "Comparison cases: $TOTAL_CASES"
echo "Algorithm runs:   $TOTAL_ALGO_RUNS"
echo "Failed cases:     $FAIL"

if [ "$FAIL" -ne 0 ]; then
    echo "Failure list: $(realpath "$FAIL_LIST")"
    exit 2
fi
