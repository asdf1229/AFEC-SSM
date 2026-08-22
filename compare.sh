#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESOURCE_RUNNER="${SCRIPT_DIR}/scripts/run_with_timeout.py"
PROGRESS_MONITOR="${SCRIPT_DIR}/scripts/monitor_progress.py"
DATA_DIR=./test/datasets
BUILD_DIR=./build
SKIP_BUILD=false
DRY_RUN=false
RESUME=false
RESULT_DIR_OPTION=""
SELECTED_ALGOS=()
SELECTED_THRESHOLDS=()
ALGO_TIMEOUT_SECONDS="${ALGO_TIMEOUT_SECONDS:-60}"
LOAD_TIMEOUT_SECONDS="${LOAD_TIMEOUT_SECONDS:-300}"
PARALLEL_JOBS=1
PROGRESS_INTERVAL_SECONDS="${PROGRESS_INTERVAL_SECONDS:-30}"
CURRENT_CHILD_PID=""
PROGRESS_PID=""
PROGRESS_STOP_FILE=""
PROGRESS_READY_FILE=""
PLAN_DIR=""
OUTPUT_LIMITS=(1000000 10000000 100000000 1000000000)
MAX_OUTPUT_LIMIT="${OUTPUT_LIMITS[${#OUTPUT_LIMITS[@]}-1]}"
DEFAULT_THRESHOLDS=(1 2 3 4 5 6)
MISSING_EDGE_TEST_GROUPS=(
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
                           Example: --algorithms afec,treespan
  -b, --build-dir DIR      Build directory. Default: ./build
  -d, --data-dir DIR       Dataset root. Default: ./test/datasets
                           The script scans synthetic/ and real_graphs/ under this root.
  -p, --parallel N         Run up to N algorithm tasks concurrently. Default: 1.
      --thresholds t1,t2   Override query-group defaults with these thresholds.
                           Example: --thresholds 2 or --thresholds 1,2,3,4,5,6
      --dry-run            Build and validate the task plan without running algorithms
                           or creating the target result directory.
      --result-dir DIR     Write results to DIR instead of a timestamped directory.
      --resume             Continue an existing --result-dir. Completed tasks,
                           including terminal failures, are not rerun.
      --skip-build         Reuse existing binaries in the build directory.
  -h, --help              Show this help message.

Environment:
  ALGO_TIMEOUT_SECONDS     Algorithm timeout after graph loading. Default: 60.
  LOAD_TIMEOUT_SECONDS     Graph-loading timeout in seconds. Default: 300.
  PROGRESS_INTERVAL_SECONDS  progress.tsv update interval. Default: 30.
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
        if [ -n "$pid" ] && [ "$pid" != "${PROGRESS_PID:-}" ]; then
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

cleanup_progress_monitor() {
    if [ -z "${PROGRESS_PID:-}" ]; then
        return 0
    fi
    if kill -0 "$PROGRESS_PID" 2>/dev/null; then
        kill -TERM "$PROGRESS_PID" 2>/dev/null || true
        wait "$PROGRESS_PID" 2>/dev/null || true
    fi
    PROGRESS_PID=""
    if [ -n "${PROGRESS_STOP_FILE:-}" ]; then
        rm -f "$PROGRESS_STOP_FILE"
    fi
    if [ -n "${PROGRESS_READY_FILE:-}" ]; then
        rm -f "$PROGRESS_READY_FILE"
    fi
}

cleanup_plan_dir() {
    if [ -z "${PLAN_DIR:-}" ] || [ ! -d "$PLAN_DIR" ]; then
        return 0
    fi

    rm -f "$PLAN_DIR/tasks.tsv" "$PLAN_DIR/cases.tsv" \
        "$PLAN_DIR/run_config.tsv" "$PLAN_DIR/pending_tasks.tsv"
    rmdir "$PLAN_DIR" 2>/dev/null || true
    PLAN_DIR=""
}

on_exit() {
    cleanup_progress_monitor
    cleanup_plan_dir
}

run_with_timeout() {
    local output_file="$1"
    local runner_status_file="$2"
    shift 2

    RUNNER_TIMED_OUT=""
    RUNNER_TIMEOUT_PHASE="NA"
    RUNNER_RETURN_CODE=""
    RUNNER_PEAK_RSS_KB="NA"
    RUNNER_READY_SEEN="0"
    RUNNER_LOAD_ELAPSED_MS="NA"
    RUNNER_ALGORITHM_ELAPSED_MS="NA"
    RUNNER_ERROR=""

    python3 -B "$RESOURCE_RUNNER" \
        --timeout "$ALGO_TIMEOUT_SECONDS" \
        --load-timeout "$LOAD_TIMEOUT_SECONDS" \
        --kill-after 5 \
        --output "$output_file" \
        --status "$runner_status_file" \
        -- "$@" &
    CURRENT_CHILD_PID=$!
    wait "$CURRENT_CHILD_PID"
    local rc=$?
    CURRENT_CHILD_PID=""

    if [ -f "$runner_status_file" ]; then
        local key value
        while IFS='=' read -r key value; do
            case "$key" in
                timed_out) RUNNER_TIMED_OUT="$value" ;;
                timeout_phase) RUNNER_TIMEOUT_PHASE="$value" ;;
                return_code) RUNNER_RETURN_CODE="$value" ;;
                peak_rss_kb) RUNNER_PEAK_RSS_KB="$value" ;;
                ready_seen) RUNNER_READY_SEEN="$value" ;;
                load_elapsed_ms) RUNNER_LOAD_ELAPSED_MS="$value" ;;
                algorithm_elapsed_ms) RUNNER_ALGORITHM_ELAPSED_MS="$value" ;;
                runner_error) RUNNER_ERROR="$value" ;;
            esac
        done < "$runner_status_file"
    fi
    return "$rc"
}

trap on_interrupt INT TERM
trap on_exit EXIT

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
        --thresholds)
            if [ $# -lt 2 ] ||
               ! [[ "$2" =~ ^(0|[1-9][0-9]*)(,(0|[1-9][0-9]*))*$ ]]; then
                echo "Error: --thresholds requires comma-separated nonnegative integers" >&2
                exit 1
            fi
            IFS=',' read -r -a SELECTED_THRESHOLDS <<< "$2"
            for ((threshold_index = 0; threshold_index < ${#SELECTED_THRESHOLDS[@]}; threshold_index++)); do
                for ((previous_index = 0; previous_index < threshold_index; previous_index++)); do
                    if [ "${SELECTED_THRESHOLDS[$threshold_index]}" = \
                         "${SELECTED_THRESHOLDS[$previous_index]}" ]; then
                        echo "Error: duplicate threshold: ${SELECTED_THRESHOLDS[$threshold_index]}" >&2
                        exit 1
                    fi
                done
            done
            shift 2
            ;;
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        --result-dir)
            if [ $# -lt 2 ] || [ -z "$2" ]; then
                echo "Error: --result-dir requires a directory" >&2
                exit 1
            fi
            RESULT_DIR_OPTION="$2"
            shift 2
            ;;
        --resume)
            RESUME=true
            shift
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

if ! command -v python3 >/dev/null 2>&1; then
    echo "Error: python3 is required for timeout and peak-memory accounting" >&2
    exit 1
fi
if ! [[ "$ALGO_TIMEOUT_SECONDS" =~ ^[0-9]+([.][0-9]+)?$ ]] ||
   [[ "$ALGO_TIMEOUT_SECONDS" =~ ^0+([.]0+)?$ ]]; then
    echo "Error: ALGO_TIMEOUT_SECONDS must be a positive number" >&2
    exit 1
fi
if ! [[ "$LOAD_TIMEOUT_SECONDS" =~ ^[0-9]+([.][0-9]+)?$ ]] ||
   [[ "$LOAD_TIMEOUT_SECONDS" =~ ^0+([.]0+)?$ ]]; then
    echo "Error: LOAD_TIMEOUT_SECONDS must be a positive number" >&2
    exit 1
fi
if ! [[ "$PROGRESS_INTERVAL_SECONDS" =~ ^[0-9]+([.][0-9]+)?$ ]] ||
   [[ "$PROGRESS_INTERVAL_SECONDS" =~ ^0+([.]0+)?$ ]]; then
    echo "Error: PROGRESS_INTERVAL_SECONDS must be a positive number" >&2
    exit 1
fi
if [ ! -f "$RESOURCE_RUNNER" ]; then
    echo "Error: resource runner not found: $RESOURCE_RUNNER" >&2
    exit 1
fi
if [ ! -f "$PROGRESS_MONITOR" ]; then
    echo "Error: progress monitor not found: $PROGRESS_MONITOR" >&2
    exit 1
fi

if [ "$RESUME" = true ] && [ -z "$RESULT_DIR_OPTION" ]; then
    echo "Error: --resume requires --result-dir DIR" >&2
    exit 1
fi

DATA_DIR=$(python3 -c 'import os, sys; print(os.path.realpath(sys.argv[1]))' "$DATA_DIR")
BUILD_DIR=$(python3 -c 'import os, sys; print(os.path.realpath(sys.argv[1]))' "$BUILD_DIR")
timestamp=$(date +"%Y%m%d_%H%M%S")
if [ -n "$RESULT_DIR_OPTION" ]; then
    RESULT_DIR=$(python3 -c 'import os, sys; print(os.path.realpath(sys.argv[1]))' "$RESULT_DIR_OPTION")
else
    RESULT_DIR=$(python3 -c 'import os, sys; print(os.path.realpath(sys.argv[1]))' "./result/${timestamp}")
fi

if [ "$RESUME" = true ]; then
    if [ ! -d "$RESULT_DIR" ]; then
        echo "Error: resume result directory does not exist: $RESULT_DIR" >&2
        exit 1
    fi
    for required in tasks.tsv cases.tsv run_config.tsv; do
        if [ ! -f "$RESULT_DIR/$required" ]; then
            echo "Error: resume metadata is missing: $RESULT_DIR/$required" >&2
            exit 1
        fi
    done
    if [ ! -d "$RESULT_DIR/task_status" ]; then
        echo "Error: resume task status directory is missing: $RESULT_DIR/task_status" >&2
        exit 1
    fi
    if [ ! -d "$RESULT_DIR/bin" ]; then
        echo "Error: resume binary snapshot is missing: $RESULT_DIR/bin" >&2
        exit 1
    fi
    # Rebuilding here could mix different binaries in one result set.  Resume
    # always reuses the binaries whose hashes are stored in run_config.tsv.
    SKIP_BUILD=true
else
    if [ -e "$RESULT_DIR" ] && [ ! -d "$RESULT_DIR" ]; then
        echo "Error: result path exists and is not a directory: $RESULT_DIR" >&2
        exit 1
    fi
    if [ -d "$RESULT_DIR" ] &&
       [ -n "$(find "$RESULT_DIR" -mindepth 1 -maxdepth 1 -print -quit)" ]; then
        echo "Error: result directory is not empty: $RESULT_DIR" >&2
        echo "Use --resume with the same --result-dir to continue it." >&2
        exit 1
    fi
fi

if [ "$DRY_RUN" = false ]; then
    mkdir -p "$RESULT_DIR"
    OUT_LOG="${RESULT_DIR}/out.log"
    if [ "$RESUME" = true ]; then
        echo "out.log will be appended to: $OUT_LOG"
        exec >> "$OUT_LOG" 2>&1
        echo
        echo "===== RESUMING RESULTS IN: ${RESULT_DIR} ====="
    else
        echo "out.log will be saved to: $OUT_LOG"
        exec > "$OUT_LOG" 2>&1
        echo "===== RESULTS WILL BE SAVED TO: ${RESULT_DIR} ====="
    fi
else
    echo "===== DRY RUN: NO ALGORITHM TASKS WILL BE EXECUTED ====="
    echo "Planned result directory: $RESULT_DIR"
fi

if [ "$SKIP_BUILD" = false ]; then
    echo "===== Building algorithm executables ====="
    "${SCRIPT_DIR}/build.sh" "$BUILD_DIR"
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

    # Put the complete AFEC configuration first so
    # it is also the correctness baseline and the default primary algorithm
    # for table/report tools that preserve summary-column order.
    printf '%s\n' "${resolved[@]}" | \
        awk '
            {
                executable = $0
                name = executable
                sub(/^.*\//, "", name)
                sub(/\.exe$/, "", name)
                rank = (name == "ssm_afec") ? 0 : 1
                printf "%d\t%s\n", rank, executable
            }
        ' | sort -t $'\t' -k1,1n -k2,2 | cut -f2-
}

normalize_algo_key() {
    local name
    name="$(basename "$1")"
    name="${name#ssm_}"
    name="${name%.exe}"
    printf '%s\n' "$name"
}

snapshot_selected_algorithms() {
    local idx source destination temporary

    mkdir -p "$SNAPSHOT_BIN_DIR"
    for idx in "${!SOURCE_ALGO_EXECUTABLES[@]}"; do
        source="${SOURCE_ALGO_EXECUTABLES[$idx]}"
        destination="${ALGO_EXECUTABLES[$idx]}"
        temporary="${destination}.tmp.$$.${RANDOM}"
        cp -p "$source" "$temporary"
        chmod u+x "$temporary"
        mv -f "$temporary" "$destination"
    done
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

file_sha256() {
    python3 - "$1" <<'PY'
import hashlib
import sys

digest = hashlib.sha256()
with open(sys.argv[1], "rb") as handle:
    for chunk in iter(lambda: handle.read(1024 * 1024), b""):
        digest.update(chunk)
print(digest.hexdigest())
PY
}

write_run_config() {
    local config_file="$1"
    local idx

    {
        printf 'format\tSSM_RUN_CONFIG_V3\n'
        printf 'data_dir\t%s\n' "$DATA_DIR"
        printf 'build_dir\t%s\n' "$BUILD_DIR"
        printf 'snapshot_bin_dir\t%s\n' "$SNAPSHOT_BIN_DIR"
        printf 'algorithm_timeout_seconds\t%s\n' "$ALGO_TIMEOUT_SECONDS"
        printf 'load_timeout_seconds\t%s\n' "$LOAD_TIMEOUT_SECONDS"
        printf 'parallel_jobs\t%s\n' "$PARALLEL_JOBS"
        if [ ${#SELECTED_THRESHOLDS[@]} -gt 0 ]; then
            printf 'threshold_mode\texplicit\n'
            printf 'requested_thresholds\t%s\n' "${SELECTED_THRESHOLDS[*]}"
        else
            printf 'threshold_mode\tquery_group_defaults\n'
            printf 'requested_thresholds\t%s\n' "${DEFAULT_THRESHOLDS[*]}"
            printf 'threshold_group_overrides\t%s\n' "${MISSING_EDGE_TEST_GROUPS[*]}"
        fi
        printf 'output_limits\t%s\n' "${OUTPUT_LIMITS[*]}"
        printf 'resource_runner\t%s\t%s\n' "$RESOURCE_RUNNER" \
            "$(file_sha256 "$RESOURCE_RUNNER")"
        for idx in "${!ALGO_EXECUTABLES[@]}"; do
            printf 'algorithm\t%s\t%s\t%s\n' "${ALGO_KEYS[$idx]}" \
                "${ALGO_EXECUTABLES[$idx]}" \
                "$(file_sha256 "${ALGO_HASH_EXECUTABLES[$idx]}")"
        done
    } > "$config_file"
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

    if [ ${#SELECTED_THRESHOLDS[@]} -gt 0 ]; then
        printf '%s\n' "${SELECTED_THRESHOLDS[@]}"
        return 0
    fi

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
    SUMMARY_PREPROCESSING_MS=""
    SUMMARY_SEARCH_MS=""
    SUMMARY_PEAK_RSS_KB=""
    SUMMARY_FILTER_CANDIDATES=""
    SUMMARY_RECURSION_CALLS=""
    SUMMARY_REQUESTED_THRESHOLD=""
    SUMMARY_EFFECTIVE_THRESHOLD=""
    SUMMARY_OUTPUT_LIMIT=""
    SUMMARY_OUTPUT_LIMIT_REACHED=""

    for token in $line; do
        case "$token" in
            algorithm=*) SUMMARY_ALGO="${token#algorithm=}" ;;
            count=*) SUMMARY_COUNT="${token#count=}" ;;
            load_ms=*) SUMMARY_LOAD_MS="${token#load_ms=}" ;;
            run_ms=*) SUMMARY_RUN_MS="${token#run_ms=}" ;;
            total_ms=*) SUMMARY_TOTAL_MS="${token#total_ms=}" ;;
            preprocessing_ms=*) SUMMARY_PREPROCESSING_MS="${token#preprocessing_ms=}" ;;
            search_ms=*) SUMMARY_SEARCH_MS="${token#search_ms=}" ;;
            peak_rss_kb=*) SUMMARY_PEAK_RSS_KB="${token#peak_rss_kb=}" ;;
            filter_candidates=*) SUMMARY_FILTER_CANDIDATES="${token#filter_candidates=}" ;;
            recursion_calls=*) SUMMARY_RECURSION_CALLS="${token#recursion_calls=}" ;;
            requested_threshold=*) SUMMARY_REQUESTED_THRESHOLD="${token#requested_threshold=}" ;;
            effective_threshold=*) SUMMARY_EFFECTIVE_THRESHOLD="${token#effective_threshold=}" ;;
            output_limit=*) SUMMARY_OUTPUT_LIMIT="${token#output_limit=}" ;;
            output_limit_reached=*) SUMMARY_OUTPUT_LIMIT_REACHED="${token#output_limit_reached=}" ;;
        esac
    done

    if [ -z "$SUMMARY_ALGO" ] || [ -z "$SUMMARY_COUNT" ] ||
       [ -z "$SUMMARY_LOAD_MS" ] ||
       [ -z "$SUMMARY_RUN_MS" ] || [ -z "$SUMMARY_TOTAL_MS" ] ||
       [ -z "$SUMMARY_PREPROCESSING_MS" ] || [ -z "$SUMMARY_SEARCH_MS" ] ||
       [ -z "$SUMMARY_PEAK_RSS_KB" ] ||
       [ -z "$SUMMARY_FILTER_CANDIDATES" ] ||
       [ -z "$SUMMARY_RECURSION_CALLS" ] ||
       [ -z "$SUMMARY_REQUESTED_THRESHOLD" ] ||
       [ -z "$SUMMARY_EFFECTIVE_THRESHOLD" ] ||
       [ -z "$SUMMARY_OUTPUT_LIMIT" ] ||
       [ -z "$SUMMARY_OUTPUT_LIMIT_REACHED" ]; then
        return 1
    fi
}

parse_ready_line() {
    local line="$1"
    READY_ALGO=""
    READY_LOAD_MS=""
    READY_REQUESTED_THRESHOLD=""
    READY_EFFECTIVE_THRESHOLD=""

    local token
    for token in $line; do
        case "$token" in
            algorithm=*) READY_ALGO="${token#algorithm=}" ;;
            load_ms=*) READY_LOAD_MS="${token#load_ms=}" ;;
            requested_threshold=*) READY_REQUESTED_THRESHOLD="${token#requested_threshold=}" ;;
            effective_threshold=*) READY_EFFECTIVE_THRESHOLD="${token#effective_threshold=}" ;;
        esac
    done

    [ -n "$READY_ALGO" ] && [ -n "$READY_LOAD_MS" ] &&
        [ -n "$READY_REQUESTED_THRESHOLD" ] &&
        [ -n "$READY_EFFECTIVE_THRESHOLD" ]
}

parse_preprocessed_line() {
    local line="$1"
    PREPROCESSED_ALGO=""
    PREPROCESSED_MS=""
    PREPROCESSED_FILTER_CANDIDATES=""

    local token
    for token in $line; do
        case "$token" in
            algorithm=*) PREPROCESSED_ALGO="${token#algorithm=}" ;;
            preprocessing_ms=*) PREPROCESSED_MS="${token#preprocessing_ms=}" ;;
            filter_candidates=*) PREPROCESSED_FILTER_CANDIDATES="${token#filter_candidates=}" ;;
        esac
    done

    [ -n "$PREPROCESSED_ALGO" ] && [ -n "$PREPROCESSED_MS" ] &&
        [ -n "$PREPROCESSED_FILTER_CANDIDATES" ]
}

parse_checkpoint_line() {
    local line="$1"
    CHECKPOINT_ALGO=""
    CHECKPOINT_LIMIT=""
    CHECKPOINT_COUNT=""
    CHECKPOINT_RUN_MS=""
    CHECKPOINT_PREPROCESSING_MS=""
    CHECKPOINT_SEARCH_MS=""
    CHECKPOINT_PEAK_RSS_KB=""
    CHECKPOINT_FILTER_CANDIDATES=""
    CHECKPOINT_RECURSION_CALLS=""

    local token
    for token in $line; do
        case "$token" in
            algorithm=*) CHECKPOINT_ALGO="${token#algorithm=}" ;;
            output_limit=*) CHECKPOINT_LIMIT="${token#output_limit=}" ;;
            count=*) CHECKPOINT_COUNT="${token#count=}" ;;
            run_ms=*) CHECKPOINT_RUN_MS="${token#run_ms=}" ;;
            preprocessing_ms=*) CHECKPOINT_PREPROCESSING_MS="${token#preprocessing_ms=}" ;;
            search_ms=*) CHECKPOINT_SEARCH_MS="${token#search_ms=}" ;;
            peak_rss_kb=*) CHECKPOINT_PEAK_RSS_KB="${token#peak_rss_kb=}" ;;
            filter_candidates=*) CHECKPOINT_FILTER_CANDIDATES="${token#filter_candidates=}" ;;
            recursion_calls=*) CHECKPOINT_RECURSION_CALLS="${token#recursion_calls=}" ;;
        esac
    done

    if [ -z "$CHECKPOINT_ALGO" ] || [ -z "$CHECKPOINT_LIMIT" ] ||
       [ -z "$CHECKPOINT_COUNT" ] || [ -z "$CHECKPOINT_RUN_MS" ] ||
       [ -z "$CHECKPOINT_PREPROCESSING_MS" ] ||
       [ -z "$CHECKPOINT_SEARCH_MS" ] ||
       [ -z "$CHECKPOINT_PEAK_RSS_KB" ] ||
       [ -z "$CHECKPOINT_FILTER_CANDIDATES" ] ||
       [ -z "$CHECKPOINT_RECURSION_CALLS" ]; then
        return 1
    fi
}

is_output_limit() {
    local candidate="$1"
    local output_limit
    for output_limit in "${OUTPUT_LIMITS[@]}"; do
        if [ "$candidate" = "$output_limit" ]; then
            return 0
        fi
    done
    return 1
}

is_unsigned_integer() {
    [[ "$1" =~ ^[0-9]+$ ]]
}

is_nonnegative_decimal() {
    [[ "$1" =~ ^[0-9]+([.][0-9]+)?$ ]]
}

is_optional_unsigned_integer() {
    [ "$1" = "NA" ] || is_unsigned_integer "$1"
}

is_optional_nonnegative_decimal() {
    [ "$1" = "NA" ] || is_nonnegative_decimal "$1"
}

is_integer() {
    [[ "$1" =~ ^-?[0-9]+$ ]]
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

    local temporary="${status_file}.tmp.$$.${RANDOM}"
    printf '%s\n' "$(IFS=$'\t'; echo "${fields[*]}")" > "$temporary"
    mv -f "$temporary" "$status_file"
}

read_task_status() {
    local status_file="$1"
    local line=""
    local fields=()

    if ! IFS= read -r line < "$status_file" && [ -z "$line" ]; then
        return 1
    fi
    IFS=$'\t' read -r -a fields <<< "$line"
    if [ "${#fields[@]}" -ne 19 ] || [ "${fields[0]}" != "SSM_TASK_V2" ]; then
        return 1
    fi

    TASK_STATUS="${fields[1]}"
    TASK_RETURN_CODE="${fields[2]}"
    TASK_FINAL_AVAILABLE="${fields[3]}"
    TASK_READY_SEEN="${fields[4]}"
    TASK_COUNT="${fields[5]}"
    TASK_LOAD_MS="${fields[6]}"
    TASK_RUN_MS="${fields[7]}"
    TASK_OBSERVED_RUN_MS="${fields[8]}"
    TASK_PEAK_RSS_KB="${fields[9]}"
    TASK_FILTER_CANDIDATES="${fields[10]}"
    TASK_RECURSION_CALLS="${fields[11]}"
    TASK_PREPROCESSING_MS="${fields[12]}"
    TASK_SEARCH_MS="${fields[13]}"
    TASK_EFFECTIVE_THRESHOLD="${fields[14]}"
    TASK_TIMEOUT_PHASE="${fields[15]}"
    TASK_OUTPUT="${fields[16]}"
    TASK_LABEL="${fields[17]}"
    TASK_CHECKPOINT_PROTOCOL="${fields[18]}"

    case "$TASK_STATUS" in
        OK|Timeout|LoadTimeout|Killed|RunError|RunnerError|ParseError|ConfigError) ;;
        *) return 1 ;;
    esac
    is_integer "$TASK_RETURN_CODE" || return 1
    [[ "$TASK_FINAL_AVAILABLE" =~ ^[01]$ ]] || return 1
    [[ "$TASK_READY_SEEN" =~ ^[01]$ ]] || return 1
    is_optional_unsigned_integer "$TASK_COUNT" || return 1
    is_optional_nonnegative_decimal "$TASK_LOAD_MS" || return 1
    is_optional_nonnegative_decimal "$TASK_RUN_MS" || return 1
    is_optional_nonnegative_decimal "$TASK_OBSERVED_RUN_MS" || return 1
    is_optional_unsigned_integer "$TASK_PEAK_RSS_KB" || return 1
    is_optional_unsigned_integer "$TASK_FILTER_CANDIDATES" || return 1
    is_optional_unsigned_integer "$TASK_RECURSION_CALLS" || return 1
    is_optional_nonnegative_decimal "$TASK_PREPROCESSING_MS" || return 1
    is_optional_nonnegative_decimal "$TASK_SEARCH_MS" || return 1
    is_optional_unsigned_integer "$TASK_EFFECTIVE_THRESHOLD" || return 1
    case "$TASK_TIMEOUT_PHASE" in NA|load|algorithm) ;; *) return 1 ;; esac
    case "$TASK_CHECKPOINT_PROTOCOL" in valid|invalid) ;; *) return 1 ;; esac
    [ -n "$TASK_OUTPUT" ] && [ -n "$TASK_LABEL" ] || return 1

    if [ "$TASK_STATUS" = "OK" ] &&
       { [ "$TASK_FINAL_AVAILABLE" != "1" ] || [ "$TASK_READY_SEEN" != "1" ] ||
         [ "$TASK_CHECKPOINT_PROTOCOL" != "valid" ]; }; then
        return 1
    fi
    if [ "$TASK_FINAL_AVAILABLE" = "1" ]; then
        is_unsigned_integer "$TASK_COUNT" || return 1
        is_nonnegative_decimal "$TASK_RUN_MS" || return 1
        is_nonnegative_decimal "$TASK_PREPROCESSING_MS" || return 1
        is_nonnegative_decimal "$TASK_SEARCH_MS" || return 1
    fi
    return 0
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

    local status_file runner_status_file rc return_code status summary_line ready_line
    local preprocessed_line
    local final_available count load_value run_value peak_rss filter_candidates recursion
    local observed_run_value preprocessing_value search_value effective_threshold
    local timeout_phase label summary_error ready_output_valid
    status_file="$(task_status_file "$case_id" "$algo")"
    runner_status_file="${status_file}.runner"

    set +e
    run_with_timeout "$out" "$runner_status_file" \
        "$exe" -d "$gfile" -q "$qfile" -t "$threshold"
    rc=$?
    set -e

    parse_output_metadata "$out" "$algo"
    label="$OUTPUT_LABEL"
    recursion="$OUTPUT_RECURSION_CALLS"
    filter_candidates="NA"
    final_available="0"
    count="NA"
    load_value="$RUNNER_LOAD_ELAPSED_MS"
    run_value="NA"
    observed_run_value="$RUNNER_ALGORITHM_ELAPSED_MS"
    peak_rss="$RUNNER_PEAK_RSS_KB"
    preprocessing_value="NA"
    search_value="NA"
    effective_threshold="NA"
    timeout_phase="${RUNNER_TIMEOUT_PHASE:-NA}"
    return_code="${RUNNER_RETURN_CODE:-$rc}"

    ready_output_valid="0"
    ready_line=$(grep '^SSM_READY ' "$out" | tail -n 1 || true)
    if parse_ready_line "$ready_line" &&
       [ "$READY_ALGO" = "$algo" ] &&
       is_nonnegative_decimal "$READY_LOAD_MS" &&
       [ "$READY_REQUESTED_THRESHOLD" = "$threshold" ] &&
       is_unsigned_integer "$READY_EFFECTIVE_THRESHOLD"; then
        load_value="$READY_LOAD_MS"
        effective_threshold="$READY_EFFECTIVE_THRESHOLD"
        ready_output_valid="1"
    fi

    preprocessed_line=$(grep '^SSM_PREPROCESSED ' "$out" | tail -n 1 || true)
    if parse_preprocessed_line "$preprocessed_line" &&
       [ "$PREPROCESSED_ALGO" = "$algo" ] &&
       is_nonnegative_decimal "$PREPROCESSED_MS" &&
       is_optional_unsigned_integer "$PREPROCESSED_FILTER_CANDIDATES"; then
        preprocessing_value="$PREPROCESSED_MS"
        filter_candidates="$PREPROCESSED_FILTER_CANDIDATES"
    fi

    summary_error="ParseError"
    summary_line=$(grep '^SSM_SUMMARY ' "$out" | tail -n 1 || true)
    if parse_summary_line "$summary_line"; then
        if [ "$SUMMARY_ALGO" = "$algo" ] &&
           is_unsigned_integer "$SUMMARY_COUNT" &&
           is_nonnegative_decimal "$SUMMARY_LOAD_MS" &&
           is_nonnegative_decimal "$SUMMARY_RUN_MS" &&
           is_nonnegative_decimal "$SUMMARY_TOTAL_MS" &&
           is_nonnegative_decimal "$SUMMARY_PREPROCESSING_MS" &&
           is_nonnegative_decimal "$SUMMARY_SEARCH_MS" &&
           is_unsigned_integer "$SUMMARY_PEAK_RSS_KB" &&
           is_optional_unsigned_integer "$SUMMARY_FILTER_CANDIDATES" &&
           is_optional_unsigned_integer "$SUMMARY_RECURSION_CALLS" &&
           [ "$SUMMARY_REQUESTED_THRESHOLD" = "$threshold" ] &&
           is_unsigned_integer "$SUMMARY_EFFECTIVE_THRESHOLD" &&
           [[ "$SUMMARY_OUTPUT_LIMIT_REACHED" =~ ^[01]$ ]]; then
            local expected_limit_reached="0"
            if (( SUMMARY_COUNT >= MAX_OUTPUT_LIMIT )); then
                expected_limit_reached="1"
            fi
            if [ "$SUMMARY_OUTPUT_LIMIT" = "$MAX_OUTPUT_LIMIT" ] &&
               [ "$SUMMARY_OUTPUT_LIMIT_REACHED" = "$expected_limit_reached" ] &&
               [ "$RUNNER_READY_SEEN" = "1" ]; then
                final_available="1"
                count="$SUMMARY_COUNT"
                load_value="$SUMMARY_LOAD_MS"
                run_value="$SUMMARY_RUN_MS"
                preprocessing_value="$SUMMARY_PREPROCESSING_MS"
                search_value="$SUMMARY_SEARCH_MS"
                filter_candidates="$SUMMARY_FILTER_CANDIDATES"
                recursion="$SUMMARY_RECURSION_CALLS"
                effective_threshold="$SUMMARY_EFFECTIVE_THRESHOLD"
                summary_error=""
                if ! is_unsigned_integer "$peak_rss" ||
                   (( SUMMARY_PEAK_RSS_KB > peak_rss )); then
                    peak_rss="$SUMMARY_PEAK_RSS_KB"
                fi
            else
                summary_error="ConfigError"
            fi
        fi
    fi

    if [ "$rc" -eq 0 ] && [ "$RUNNER_READY_SEEN" != "1" ]; then
        status="RunnerError"
        final_available="0"
    elif [ "$rc" -eq 0 ] && [ "$ready_output_valid" != "1" ]; then
        status="ParseError"
        final_available="0"
    elif [ "$rc" -eq 0 ]; then
        if [ "$final_available" = "1" ]; then
            status="OK"
        else
            status="$summary_error"
        fi
    elif [ -n "$RUNNER_ERROR" ]; then
        status="RunnerError"
    elif [ "$RUNNER_TIMED_OUT" = "1" ] && [ "$timeout_phase" = "load" ]; then
        status="LoadTimeout"
    elif [ "$RUNNER_TIMED_OUT" = "1" ]; then
        status="Timeout"
    elif [[ "$RUNNER_RETURN_CODE" =~ ^-[0-9]+$ ]]; then
        status="Killed"
    else
        status="RunError"
    fi

    local checkpoint_protocol="invalid"
    if [ "$RUNNER_READY_SEEN" = "1" ] && [ "$ready_output_valid" = "1" ]; then
        checkpoint_protocol="valid"
    fi
    write_task_status "$status_file" "SSM_TASK_V2" "$status" "$return_code" \
        "$final_available" "${RUNNER_READY_SEEN:-0}" "$count" "$load_value" \
        "$run_value" "$observed_run_value" "$peak_rss" "$filter_candidates" \
        "$recursion" "$preprocessing_value" "$search_value" \
        "$effective_threshold" "$timeout_phase" "$abs_out" "$label" \
        "$checkpoint_protocol"
    return 0
}

run_all_tasks() {
    local launched=0
    local case_id algo_index exe algo gfile qfile threshold out abs_out
    local sem_fifo="${RESULT_DIR}/parallel.sem"
    local token i pid
    local task_pids=()

    if [ "$TOTAL_PENDING_RUNS" -eq 0 ]; then
        echo "No pending algorithm tasks. Rebuilding summaries from saved statuses."
        return 0
    fi

    echo "Running ${TOTAL_PENDING_RUNS} pending algorithm tasks with parallel=${PARALLEL_JOBS}"

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
        task_pids+=("$!")

        if [ $((launched % PARALLEL_JOBS)) -eq 0 ] || [ "$launched" -eq "$TOTAL_PENDING_RUNS" ]; then
            echo "  Launched ${launched}/${TOTAL_PENDING_RUNS} pending algorithm tasks"
        fi
    done < "$RUN_MANIFEST"

    for pid in "${task_pids[@]}"; do
        wait "$pid" || true
    done
    exec 9>&-
    exec 9<&-
}

start_progress_monitor() {
    local readiness_checks=0
    PROGRESS_STOP_FILE="${RESULT_DIR}/.progress.stop.$$"
    PROGRESS_READY_FILE="${RESULT_DIR}/.progress.ready.$$"
    rm -f "$PROGRESS_STOP_FILE" "$PROGRESS_READY_FILE"
    python3 -B "$PROGRESS_MONITOR" \
        --tasks "$TASK_MANIFEST" \
        --status-dir "$STATUS_DIR" \
        --output "${RESULT_DIR}/progress.tsv" \
        --stop-file "$PROGRESS_STOP_FILE" \
        --ready-file "$PROGRESS_READY_FILE" \
        --baseline-completed "$RESUMED_TASKS_SKIPPED" \
        --interval "$PROGRESS_INTERVAL_SECONDS" &
    PROGRESS_PID=$!

    while [ ! -f "$PROGRESS_READY_FILE" ]; do
        readiness_checks=$((readiness_checks + 1))
        if [ "$readiness_checks" -ge 200 ] ||
           ! kill -0 "$PROGRESS_PID" 2>/dev/null; then
            wait "$PROGRESS_PID" 2>/dev/null || true
            PROGRESS_PID=""
            echo "Error: progress monitor failed to initialize" >&2
            return 1
        fi
        sleep 0.05
    done
    rm -f "$PROGRESS_READY_FILE"
}

finish_progress_monitor() {
    local monitor_rc=0
    if [ -z "${PROGRESS_PID:-}" ]; then
        return 0
    fi
    if kill -0 "$PROGRESS_PID" 2>/dev/null; then
        : > "$PROGRESS_STOP_FILE"
        kill -TERM "$PROGRESS_PID" 2>/dev/null || true
    fi
    wait "$PROGRESS_PID" || monitor_rc=$?
    PROGRESS_PID=""
    rm -f "$PROGRESS_STOP_FILE" "$PROGRESS_READY_FILE"
    if [ "$monitor_rc" -ne 0 ]; then
        echo "Error: progress monitor exited with status $monitor_rc" >&2
        return 1
    fi
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
        LoadTimeout)
            printf "%s\t%s\t%s\t%s\tt=%s\t%s\tLoadTimeout(%ss)\t%s\n" \
                "$dataset_group" "$dirname" "$query_group_label" "$qname" "$threshold" "$algo" "$LOAD_TIMEOUT_SECONDS" "$output" >> "$FAIL_LIST"
            ;;
        RunError|RunnerError|ParseError|MissingStatus)
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
        declare -A load_ms=()
        declare -A run_ms=()
        declare -A observed_run_ms=()
        declare -A preprocessing_ms=()
        declare -A search_ms=()
        declare -A peak_rss_kb=()
        declare -A filter_candidates=()
        declare -A outputs=()
        declare -A statuses=()
        declare -A return_codes=()
        declare -A final_availabilities=()
        declare -A timeout_phases=()
        declare -A effective_thresholds=()
        declare -A checkpoint_protocols=()
        declare -A checkpoint_parse_statuses=()
        declare -A display_names=()
        declare -A recursion_calls=()
        declare -A checkpoint_counts=()
        declare -A checkpoint_run_ms=()
        declare -A checkpoint_preprocessing_ms=()
        declare -A checkpoint_search_ms=()
        declare -A checkpoint_peak_rss_kb=()
        declare -A checkpoint_filter_candidates=()
        declare -A checkpoint_recursion_calls=()
        declare -A checkpoint_sources=()

        local case_failed=false
        local case_reason=""
        local idx algo status_file status return_code final_available ready_seen count load_value
        local run_value observed_run_value preprocessing_value search_value
        local peak_rss filter_value recursion effective_threshold checkpoint_protocol
        local output timeout_phase label
        local checkpoint_line checkpoint_key output_limit checkpoint_stream_valid seen_gap

        for idx in "${!ALGO_KEYS[@]}"; do
            algo="${ALGO_KEYS[$idx]}"
            status_file="$(task_status_file "$case_id" "$algo")"

            if [ ! -f "$status_file" ]; then
                status="MissingStatus"
                return_code="NA"
                final_available="0"
                ready_seen="0"
                count="NA"
                load_value="NA"
                run_value="NA"
                observed_run_value="NA"
                preprocessing_value="NA"
                search_value="NA"
                peak_rss="NA"
                filter_value="NA"
                recursion="NA"
                effective_threshold="NA"
                checkpoint_protocol="invalid"
                timeout_phase="NA"
                output="NA"
                label="$algo"
            elif read_task_status "$status_file"; then
                status="$TASK_STATUS"
                return_code="$TASK_RETURN_CODE"
                final_available="$TASK_FINAL_AVAILABLE"
                ready_seen="$TASK_READY_SEEN"
                count="$TASK_COUNT"
                load_value="$TASK_LOAD_MS"
                run_value="$TASK_RUN_MS"
                observed_run_value="$TASK_OBSERVED_RUN_MS"
                peak_rss="$TASK_PEAK_RSS_KB"
                filter_value="$TASK_FILTER_CANDIDATES"
                recursion="$TASK_RECURSION_CALLS"
                preprocessing_value="$TASK_PREPROCESSING_MS"
                search_value="$TASK_SEARCH_MS"
                effective_threshold="$TASK_EFFECTIVE_THRESHOLD"
                timeout_phase="$TASK_TIMEOUT_PHASE"
                output="$TASK_OUTPUT"
                label="$TASK_LABEL"
                checkpoint_protocol="$TASK_CHECKPOINT_PROTOCOL"
            else
                status="StatusParseError"
                return_code="NA"
                final_available="0"
                ready_seen="0"
                count="NA"
                load_value="NA"
                run_value="NA"
                observed_run_value="NA"
                preprocessing_value="NA"
                search_value="NA"
                peak_rss="NA"
                filter_value="NA"
                recursion="NA"
                effective_threshold="NA"
                checkpoint_protocol="invalid"
                timeout_phase="NA"
                output="NA"
                label="$algo"
            fi

            checkpoint_stream_valid=true
            if [ "$checkpoint_protocol" = "valid" ] &&
               [ "$output" != "NA" ] && [ -f "$output" ]; then
                while IFS= read -r checkpoint_line; do
                    if ! parse_checkpoint_line "$checkpoint_line"; then
                        continue
                    fi
                    if [ "$CHECKPOINT_ALGO" != "$algo" ] ||
                       ! is_output_limit "$CHECKPOINT_LIMIT"; then
                        continue
                    fi
                    if ! is_unsigned_integer "$CHECKPOINT_COUNT" ||
                       [ "$CHECKPOINT_COUNT" != "$CHECKPOINT_LIMIT" ] ||
                       ! is_nonnegative_decimal "$CHECKPOINT_RUN_MS" ||
                       ! is_nonnegative_decimal "$CHECKPOINT_PREPROCESSING_MS" ||
                       ! is_nonnegative_decimal "$CHECKPOINT_SEARCH_MS" ||
                       ! is_unsigned_integer "$CHECKPOINT_PEAK_RSS_KB" ||
                       ! is_optional_unsigned_integer "$CHECKPOINT_FILTER_CANDIDATES" ||
                       ! is_optional_unsigned_integer "$CHECKPOINT_RECURSION_CALLS"; then
                        continue
                    fi

                    checkpoint_key="${algo}|${CHECKPOINT_LIMIT}"
                    if [ -n "${checkpoint_counts[$checkpoint_key]+present}" ]; then
                        checkpoint_stream_valid=false
                        continue
                    fi
                    checkpoint_counts["$checkpoint_key"]="$CHECKPOINT_COUNT"
                    checkpoint_run_ms["$checkpoint_key"]="$CHECKPOINT_RUN_MS"
                    checkpoint_preprocessing_ms["$checkpoint_key"]="$CHECKPOINT_PREPROCESSING_MS"
                    checkpoint_search_ms["$checkpoint_key"]="$CHECKPOINT_SEARCH_MS"
                    checkpoint_peak_rss_kb["$checkpoint_key"]="$CHECKPOINT_PEAK_RSS_KB"
                    checkpoint_filter_candidates["$checkpoint_key"]="$CHECKPOINT_FILTER_CANDIDATES"
                    checkpoint_recursion_calls["$checkpoint_key"]="$CHECKPOINT_RECURSION_CALLS"
                    checkpoint_sources["$checkpoint_key"]="checkpoint"
                done < <(grep '^SSM_CHECKPOINT ' "$output" 2>/dev/null || true)
            fi

            # Checkpoints form a strict prefix. A later checkpoint without all
            # earlier flushed milestones indicates truncated/corrupt output.
            seen_gap=false
            for output_limit in "${OUTPUT_LIMITS[@]}"; do
                checkpoint_key="${algo}|${output_limit}"
                if [ -n "${checkpoint_counts[$checkpoint_key]+present}" ]; then
                    if [ "$seen_gap" = true ]; then
                        checkpoint_stream_valid=false
                    fi
                else
                    seen_gap=true
                fi
            done

            if [ "$checkpoint_stream_valid" = false ]; then
                for output_limit in "${OUTPUT_LIMITS[@]}"; do
                    checkpoint_key="${algo}|${output_limit}"
                    unset "checkpoint_counts[$checkpoint_key]"
                    unset "checkpoint_run_ms[$checkpoint_key]"
                    unset "checkpoint_preprocessing_ms[$checkpoint_key]"
                    unset "checkpoint_search_ms[$checkpoint_key]"
                    unset "checkpoint_peak_rss_kb[$checkpoint_key]"
                    unset "checkpoint_filter_candidates[$checkpoint_key]"
                    unset "checkpoint_recursion_calls[$checkpoint_key]"
                    unset "checkpoint_sources[$checkpoint_key]"
                done
                if [ "$status" = "OK" ]; then
                    status="ParseError"
                fi
            fi
            if [ "$checkpoint_protocol" != "valid" ]; then
                if [ "$status" = "LoadTimeout" ]; then
                    checkpoint_parse_statuses["$algo"]="NotStarted"
                else
                    checkpoint_parse_statuses["$algo"]="ProtocolError"
                fi
            elif [ "$checkpoint_stream_valid" = false ]; then
                checkpoint_parse_statuses["$algo"]="ParseError"
            else
                checkpoint_parse_statuses["$algo"]="OK"
            fi

            if [ "$final_available" = "1" ]; then
                for output_limit in "${OUTPUT_LIMITS[@]}"; do
                    checkpoint_key="${algo}|${output_limit}"
                    if [ -n "${checkpoint_counts[$checkpoint_key]+present}" ]; then
                        continue
                    fi

                    # If enumeration naturally finishes below a requested
                    # ceiling, its final sample is exactly that ceiling's
                    # result.  A missing checkpoint at or below final count is
                    # instead malformed output and must not be hidden.
                    if [[ "$count" =~ ^[0-9]+$ ]] && (( count < output_limit )); then
                        checkpoint_counts["$checkpoint_key"]="$count"
                        checkpoint_run_ms["$checkpoint_key"]="$run_value"
                        checkpoint_preprocessing_ms["$checkpoint_key"]="$preprocessing_value"
                        checkpoint_search_ms["$checkpoint_key"]="$search_value"
                        checkpoint_peak_rss_kb["$checkpoint_key"]="$peak_rss"
                        checkpoint_filter_candidates["$checkpoint_key"]="$filter_value"
                        checkpoint_recursion_calls["$checkpoint_key"]="$recursion"
                        checkpoint_sources["$checkpoint_key"]="final"
                    else
                        # A completed/final count at or above this milestone
                        # proves that its flushed checkpoint should exist.
                        # Keep earlier valid checkpoints usable, but expose the
                        # missing milestone as a checkpoint-stream parse error.
                        checkpoint_parse_statuses["$algo"]="ParseError"
                        if [ "$status" = "OK" ]; then
                            status="ParseError"
                        fi
                    fi
                done
            fi

            counts["$algo"]="$count"
            load_ms["$algo"]="$load_value"
            run_ms["$algo"]="$run_value"
            observed_run_ms["$algo"]="$observed_run_value"
            preprocessing_ms["$algo"]="$preprocessing_value"
            search_ms["$algo"]="$search_value"
            peak_rss_kb["$algo"]="$peak_rss"
            filter_candidates["$algo"]="$filter_value"
            recursion_calls["$algo"]="$recursion"
            outputs["$algo"]="$output"
            statuses["$algo"]="$status"
            return_codes["$algo"]="$return_code"
            final_availabilities["$algo"]="$final_available"
            timeout_phases["$algo"]="$timeout_phase"
            effective_thresholds["$algo"]="$effective_threshold"
            checkpoint_protocols["$algo"]="$checkpoint_protocol"
            display_names["$algo"]="${label:-$algo}"

            if [ "$status" != "OK" ]; then
                case_failed=true
                record_failure "$dataset_group" "$dirname" "$query_group_label" "$qname" "$t" "$algo" "$status" "$output"
            fi
        done

        local expected_count="NA"
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
            local algo_load="${load_ms[$algo]:-NA}"
            local algo_peak_rss="${peak_rss_kb[$algo]:-NA}"
            local algo_filter="${filter_candidates[$algo]:-NA}"
            local algo_recursion="${recursion_calls[$algo]:-NA}"

            if [ "$algo_status" = "OK" ]; then
                echo "             ${algo_label}: load ${algo_load} ms, run ${algo_run} ms, Peak RSS: ${algo_peak_rss} KiB, Filter Candidates: ${algo_filter}, Recursion Calls: ${algo_recursion}, Output: ${algo_output}"
            elif [ "$algo_status" = "Timeout" ]; then
                echo "             ${algo_label}: TIMEOUT(${ALGO_TIMEOUT_SECONDS}s), Peak RSS: ${algo_peak_rss} KiB, Recursion Calls: ${algo_recursion}, Output: ${algo_output}"
            elif [ "$algo_status" = "LoadTimeout" ]; then
                echo "             ${algo_label}: LOAD TIMEOUT(${LOAD_TIMEOUT_SECONDS}s), Peak RSS: ${algo_peak_rss} KiB, Output: ${algo_output}"
            else
                echo "             ${algo_label}: ${algo_status}, Peak RSS: ${algo_peak_rss} KiB, Recursion Calls: ${algo_recursion}, Output: ${algo_output}"
            fi
        done

        for algo in "${ALGO_KEYS[@]}"; do
            row+=("${counts[$algo]:-NA}" "${run_ms[$algo]:-NA}" "${peak_rss_kb[$algo]:-NA}" \
                "${recursion_calls[$algo]:-NA}" "${outputs[$algo]:-NA}")
            row+=("${load_ms[$algo]:-NA}" "${filter_candidates[$algo]:-NA}" \
                "${statuses[$algo]:-NA}" "${return_codes[$algo]:-NA}" \
                "${timeout_phases[$algo]:-NA}" "${final_availabilities[$algo]:-0}")
            row+=("${observed_run_ms[$algo]:-NA}" "${preprocessing_ms[$algo]:-NA}" \
                "${search_ms[$algo]:-NA}" "${effective_thresholds[$algo]:-NA}" \
                "${checkpoint_protocols[$algo]:-invalid}" \
                "${checkpoint_parse_statuses[$algo]:-ProtocolError}")
            for output_limit in "${OUTPUT_LIMITS[@]}"; do
                checkpoint_key="${algo}|${output_limit}"
                row+=("${checkpoint_counts[$checkpoint_key]:-NA}" \
                    "${checkpoint_run_ms[$checkpoint_key]:-NA}" \
                    "${checkpoint_preprocessing_ms[$checkpoint_key]:-NA}" \
                    "${checkpoint_search_ms[$checkpoint_key]:-NA}" \
                    "${checkpoint_peak_rss_kb[$checkpoint_key]:-NA}" \
                    "${checkpoint_filter_candidates[$checkpoint_key]:-NA}" \
                    "${checkpoint_recursion_calls[$checkpoint_key]:-NA}")
            done
        done
        printf '%s\n' "$(IFS=$'\t'; echo "${row[*]}")" >> "$SUMMARY_FILE"

        # Each per-limit table is self-contained and retains every case. A
        # valid target checkpoint remains usable even if the 1e9 run later
        # times out; terminal_status records that later outcome separately.
        for output_limit in "${OUTPUT_LIMITS[@]}"; do
            local limit_all_available=true
            local limit_available_count=0
            local limit_count_mismatch=false
            local limit_expected_count="NA"
            local limit_fields=()

            for algo in "${ALGO_KEYS[@]}"; do
                checkpoint_key="${algo}|${output_limit}"
                local result_available result_source result_status
                local result_count result_run result_preprocessing result_search
                local result_peak result_filter result_recursion
                if [ -n "${checkpoint_counts[$checkpoint_key]+present}" ]; then
                    result_available="1"
                    result_source="${checkpoint_sources[$checkpoint_key]:-unknown}"
                    result_status="OK"
                    result_count="${checkpoint_counts[$checkpoint_key]:-NA}"
                    result_run="${checkpoint_run_ms[$checkpoint_key]:-NA}"
                    result_preprocessing="${checkpoint_preprocessing_ms[$checkpoint_key]:-NA}"
                    result_search="${checkpoint_search_ms[$checkpoint_key]:-NA}"
                    result_peak="${checkpoint_peak_rss_kb[$checkpoint_key]:-NA}"
                    result_filter="${checkpoint_filter_candidates[$checkpoint_key]:-NA}"
                    result_recursion="${checkpoint_recursion_calls[$checkpoint_key]:-NA}"
                    limit_available_count=$((limit_available_count + 1))
                    if [ "$limit_expected_count" = "NA" ]; then
                        limit_expected_count="$result_count"
                    elif [ "$result_count" != "$limit_expected_count" ]; then
                        limit_count_mismatch=true
                    fi
                else
                    result_available="0"
                    result_source="none"
                    result_status="${statuses[$algo]:-MissingStatus}"
                    result_count="NA"
                    result_run="NA"
                    result_preprocessing="NA"
                    result_search="NA"
                    result_peak="NA"
                    result_filter="NA"
                    result_recursion="NA"
                    limit_all_available=false
                fi

                local last_limit="NA"
                local last_count="NA"
                local last_run="NA"
                local last_preprocessing="NA"
                local last_search="NA"
                local last_peak="NA"
                local last_filter="NA"
                local last_recursion="NA"
                local candidate_limit
                for candidate_limit in "${OUTPUT_LIMITS[@]}"; do
                    if (( candidate_limit > output_limit )); then
                        break
                    fi
                    local candidate_key="${algo}|${candidate_limit}"
                    if [ "${checkpoint_sources[$candidate_key]:-}" = "checkpoint" ]; then
                        last_limit="$candidate_limit"
                        last_count="${checkpoint_counts[$candidate_key]:-NA}"
                        last_run="${checkpoint_run_ms[$candidate_key]:-NA}"
                        last_preprocessing="${checkpoint_preprocessing_ms[$candidate_key]:-NA}"
                        last_search="${checkpoint_search_ms[$candidate_key]:-NA}"
                        last_peak="${checkpoint_peak_rss_kb[$candidate_key]:-NA}"
                        last_filter="${checkpoint_filter_candidates[$candidate_key]:-NA}"
                        last_recursion="${checkpoint_recursion_calls[$candidate_key]:-NA}"
                    fi
                done

                limit_fields+=("$result_available" "$result_source" "$result_status" \
                    "${statuses[$algo]:-MissingStatus}" "${return_codes[$algo]:-NA}" \
                    "${timeout_phases[$algo]:-NA}" \
                    "${checkpoint_parse_statuses[$algo]:-ProtocolError}" \
                    "${effective_thresholds[$algo]:-NA}" "$result_count" \
                    "${load_ms[$algo]:-NA}" "$result_run" "$result_preprocessing" \
                    "$result_search" "$result_peak" "$result_filter" \
                    "$result_recursion" "${observed_run_ms[$algo]:-NA}" \
                    "${preprocessing_ms[$algo]:-NA}" "${search_ms[$algo]:-NA}" \
                    "${peak_rss_kb[$algo]:-NA}" "$last_limit" "$last_count" \
                    "$last_run" "$last_preprocessing" "$last_search" \
                    "$last_peak" "$last_filter" "$last_recursion" \
                    "${outputs[$algo]:-NA}")
            done

            local limit_case_status
            if [ "$limit_all_available" = true ]; then
                if [ "$limit_count_mismatch" = true ]; then
                    limit_case_status="CountMismatch"
                else
                    limit_case_status="OK"
                fi
            elif [ "$limit_available_count" -gt 0 ]; then
                limit_case_status="PARTIAL"
            else
                limit_case_status="FAIL"
            fi

            local limit_row=("$dataset_group" "$dirname" "$query_group_label" "$qname" \
                "$t" "$output_limit" "$limit_case_status" "$limit_expected_count")
            limit_row+=("${limit_fields[@]}")
            printf '%s\n' "$(IFS=$'\t'; echo "${limit_row[*]}")" \
                >> "${LIMIT_SUMMARY_FILES[$output_limit]}"
        done
    done < "$CASE_MANIFEST"
}

declare -A LIMIT_SUMMARY_FILES=()

initialize_summary_files() {
    local algo output_limit limit_summary_file
    local -a summary_header limit_header

    SUMMARY_FILE="${RESULT_DIR}/summary.tsv"
    FAIL_LIST="${RESULT_DIR}/failures.txt"
    : > "$FAIL_LIST"

    summary_header=(dataset_group dataset query_group query threshold status expected_count)
    for algo in "${ALGO_KEYS[@]}"; do
        summary_header+=("${algo}_count" "${algo}_run_ms" "${algo}_peak_rss_kb" \
            "${algo}_recursion_calls" "${algo}_output")
        summary_header+=("${algo}_load_ms" "${algo}_filter_candidates" \
            "${algo}_terminal_status" "${algo}_return_code" \
            "${algo}_timeout_phase" "${algo}_final_available")
        summary_header+=("${algo}_observed_run_ms_value" "${algo}_preprocessing_ms" \
            "${algo}_search_ms" "${algo}_effective_threshold" \
            "${algo}_checkpoint_protocol" "${algo}_checkpoint_parse_status")
        for output_limit in "${OUTPUT_LIMITS[@]}"; do
            summary_header+=("${algo}_count_at_${output_limit}" \
                "${algo}_run_ms_at_${output_limit}" \
                "${algo}_preprocessing_ms_at_${output_limit}" \
                "${algo}_search_ms_at_${output_limit}" \
                "${algo}_peak_rss_kb_at_${output_limit}" \
                "${algo}_filter_candidates_at_${output_limit}" \
                "${algo}_recursion_calls_at_${output_limit}")
        done
    done
    printf '%s\n' "$(IFS=$'\t'; echo "${summary_header[*]}")" > "$SUMMARY_FILE"

    LIMIT_SUMMARY_FILES=()
    for output_limit in "${OUTPUT_LIMITS[@]}"; do
        limit_summary_file="${RESULT_DIR}/summary_at_${output_limit}.tsv"
        LIMIT_SUMMARY_FILES["$output_limit"]="$limit_summary_file"
        limit_header=(dataset_group dataset query_group query threshold output_limit \
            case_status expected_capped_count)
        for algo in "${ALGO_KEYS[@]}"; do
            limit_header+=("${algo}_result_available" "${algo}_result_source" \
                "${algo}_result_status" "${algo}_terminal_status" \
                "${algo}_return_code" "${algo}_timeout_phase" \
                "${algo}_checkpoint_parse_status" "${algo}_effective_threshold" \
                "${algo}_count" "${algo}_load_ms" "${algo}_run_ms" \
                "${algo}_preprocessing_ms" "${algo}_search_ms" \
                "${algo}_peak_rss_kb" "${algo}_filter_candidates" \
                "${algo}_recursion_calls" "${algo}_terminal_run_ms_observed" \
                "${algo}_terminal_preprocessing_ms" \
                "${algo}_terminal_search_ms" \
                "${algo}_terminal_peak_rss_kb" \
                "${algo}_last_checkpoint_limit_value" \
                "${algo}_last_checkpoint_count_value" \
                "${algo}_last_checkpoint_run_ms_value" \
                "${algo}_last_checkpoint_preprocessing_ms_value" \
                "${algo}_last_checkpoint_search_ms_value" \
                "${algo}_last_checkpoint_peak_rss_kb_value" \
                "${algo}_last_checkpoint_filter_candidates_value" \
                "${algo}_last_checkpoint_recursion_calls_value" "${algo}_output")
        done
        printf '%s\n' "$(IFS=$'\t'; echo "${limit_header[*]}")" > "$limit_summary_file"
    done
}

SNAPSHOT_BIN_DIR="${RESULT_DIR}/bin"
SOURCE_ALGO_EXECUTABLES=()
ALGO_EXECUTABLES=()
ALGO_HASH_EXECUTABLES=()

if [ "$RESUME" = true ]; then
    mapfile -t ALGO_EXECUTABLES < <(resolve_selected_algorithms "$SNAPSHOT_BIN_DIR")
else
    mapfile -t SOURCE_ALGO_EXECUTABLES < <(resolve_selected_algorithms "$BUILD_DIR")
fi

ALGO_KEYS=()
if [ "$RESUME" = true ]; then
    for exe in "${ALGO_EXECUTABLES[@]}"; do
        ALGO_KEYS+=("$(normalize_algo_key "$exe")")
    done
    ALGO_HASH_EXECUTABLES=("${ALGO_EXECUTABLES[@]}")
else
    for exe in "${SOURCE_ALGO_EXECUTABLES[@]}"; do
        algo="$(normalize_algo_key "$exe")"
        ALGO_KEYS+=("$algo")
        ALGO_EXECUTABLES+=("${SNAPSHOT_BIN_DIR}/ssm_${algo}")
    done
    if [ "$DRY_RUN" = true ]; then
        # A dry run must leave RESULT_DIR untouched, so hash the selected build
        # products while planning the paths that a real run would snapshot.
        ALGO_HASH_EXECUTABLES=("${SOURCE_ALGO_EXECUTABLES[@]}")
    else
        echo "===== Saving selected executables to ${SNAPSHOT_BIN_DIR} ====="
        snapshot_selected_algorithms
        ALGO_HASH_EXECUTABLES=("${ALGO_EXECUTABLES[@]}")
    fi
fi

echo "Algorithms under comparison: ${ALGO_KEYS[*]}"
if [ ${#SELECTED_THRESHOLDS[@]} -gt 0 ]; then
    echo "Explicit thresholds: ${SELECTED_THRESHOLDS[*]}"
else
    echo "Thresholds: query-group defaults"
fi

TOTAL_CASES=0
TOTAL_ALGO_RUNS=0
FAIL=0

SAVED_TASK_MANIFEST="${RESULT_DIR}/tasks.tsv"
SAVED_CASE_MANIFEST="${RESULT_DIR}/cases.tsv"
SAVED_RUN_CONFIG="${RESULT_DIR}/run_config.tsv"
STATUS_DIR="${RESULT_DIR}/task_status"

if [ "$DRY_RUN" = true ] || [ "$RESUME" = true ]; then
    PLAN_DIR=$(mktemp -d "${TMPDIR:-/tmp}/ssm_compare_plan.XXXXXX")
    TASK_MANIFEST="${PLAN_DIR}/tasks.tsv"
    CASE_MANIFEST="${PLAN_DIR}/cases.tsv"
    RUN_CONFIG_FILE="${PLAN_DIR}/run_config.tsv"
else
    TASK_MANIFEST="$SAVED_TASK_MANIFEST"
    CASE_MANIFEST="$SAVED_CASE_MANIFEST"
    RUN_CONFIG_FILE="$SAVED_RUN_CONFIG"
fi

if [ "$DRY_RUN" = false ]; then
    mkdir -p "$STATUS_DIR"
fi
: > "$TASK_MANIFEST"
: > "$CASE_MANIFEST"
write_run_config "$RUN_CONFIG_FILE"

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
    if [ "$DRY_RUN" = false ]; then
        mkdir -p "${RESULT_DIR}/${dataset_group}/${dirname}"
    fi

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
            if [ "$DRY_RUN" = false ]; then
                echo "    Query: $qname"
            fi

            for t in "${THRESHOLDS[@]}"; do
                TOTAL_CASES=$((TOTAL_CASES + 1))
                case_id=$(printf 'case_%06d' "$TOTAL_CASES")
                printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
                    "$case_id" "$dataset_group" "$dirname" "$query_group_label" "$qname" "$t" >> "$CASE_MANIFEST"

                for idx in "${!ALGO_EXECUTABLES[@]}"; do
                    exe="${ALGO_EXECUTABLES[$idx]}"
                    algo="${ALGO_KEYS[$idx]}"
                    algo_result_dir="${RESULT_DIR}/${dataset_group}/${dirname}/${query_group_label}/${algo}"
                    if [ "$DRY_RUN" = false ]; then
                        mkdir -p "${algo_result_dir}"
                    fi
                    out="${algo_result_dir}/${qname}_t=${t}.txt"
                    abs_out="$out"

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

if [ "$RESUME" = true ]; then
    if ! cmp -s "$RUN_CONFIG_FILE" "$SAVED_RUN_CONFIG"; then
        echo "Error: resume configuration does not match $SAVED_RUN_CONFIG" >&2
        echo "Algorithms, binary hashes, thresholds, timeouts, parallelism, or paths changed." >&2
        exit 1
    fi
    if ! cmp -s "$CASE_MANIFEST" "$SAVED_CASE_MANIFEST"; then
        echo "Error: discovered comparison cases do not match $SAVED_CASE_MANIFEST" >&2
        exit 1
    fi
    if ! cmp -s "$TASK_MANIFEST" "$SAVED_TASK_MANIFEST"; then
        echo "Error: discovered algorithm tasks do not match $SAVED_TASK_MANIFEST" >&2
        exit 1
    fi
    TASK_MANIFEST="$SAVED_TASK_MANIFEST"
    CASE_MANIFEST="$SAVED_CASE_MANIFEST"
fi

RUN_MANIFEST="$TASK_MANIFEST"
TOTAL_PENDING_RUNS="$TOTAL_ALGO_RUNS"
RESUMED_TASKS_SKIPPED=0
if [ "$RESUME" = true ]; then
    RUN_MANIFEST="${PLAN_DIR}/pending_tasks.tsv"
    : > "$RUN_MANIFEST"
    TOTAL_PENDING_RUNS=0
    while IFS=$'\t' read -r case_id idx exe algo gfile qfile t out abs_out; do
        status_file="$(task_status_file "$case_id" "$algo")"
        if [ -f "$status_file" ] && read_task_status "$status_file" &&
           [ "$TASK_OUTPUT" = "$abs_out" ] && [ -f "$TASK_OUTPUT" ]; then
            RESUMED_TASKS_SKIPPED=$((RESUMED_TASKS_SKIPPED + 1))
            continue
        fi
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "$case_id" "$idx" "$exe" "$algo" "$gfile" "$qfile" "$t" \
            "$out" "$abs_out" >> "$RUN_MANIFEST"
        TOTAL_PENDING_RUNS=$((TOTAL_PENDING_RUNS + 1))
    done < "$TASK_MANIFEST"
    echo "Resume-complete tasks skipped: $RESUMED_TASKS_SKIPPED"
    echo "Resume-pending tasks:          $TOTAL_PENDING_RUNS"
fi

if [ "$DRY_RUN" = true ]; then
    if [ "$RESUME" = true ]; then
        echo "Dry-run resume tasks to execute: $TOTAL_PENDING_RUNS"
    fi
    echo "Dry run complete. No algorithm tasks were executed."
    echo "Target result directory was not modified: $RESULT_DIR"
    exit 0
fi

start_progress_monitor
run_all_tasks
finish_progress_monitor

echo "Aggregating results..."
initialize_summary_files
aggregate_results

echo "===== ALL DONE ====="
echo "Results saved at: $(realpath "${RESULT_DIR}")"
echo "Summary TSV: $(realpath "${SUMMARY_FILE}")"
for output_limit in "${OUTPUT_LIMITS[@]}"; do
    echo "Summary at ${output_limit}: $(realpath "${LIMIT_SUMMARY_FILES[$output_limit]}")"
done
echo "Algorithms: ${ALGO_KEYS[*]}"
echo "Comparison cases: $TOTAL_CASES"
echo "Algorithm runs:   $TOTAL_ALGO_RUNS"
echo "Parallel jobs:    $PARALLEL_JOBS"
echo "Failed cases:     $FAIL"

if [ "$FAIL" -ne 0 ]; then
    echo "Failure list: $(realpath "$FAIL_LIST")"
    exit 2
fi
