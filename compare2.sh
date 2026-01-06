#!/bin/bash
set -e

DATA_DIR=./test/dataset
EXEC_NAME=ssm_ged

#############################################
# CREATE TIMESTAMPED RESULT DIRECTORY
#############################################
timestamp=$(date +"%Y%m%d_%H%M%S")
RESULT_DIR="./result/${timestamp}"
mkdir -p "${RESULT_DIR}"

echo "===== RESULTS WILL BE SAVED TO: ${RESULT_DIR} ====="

#############################################
# BUILD VERSION 1 & 2 (via build.sh)
#############################################
echo "===== Building Version 1 (Approximate_Matching) via build.sh ====="
./build.sh build_v1
cp "build_v1/${EXEC_NAME}" ./exec_v1

echo "===== Building Version 2 (Approximate_Matching_V2) via build.sh ====="
./build.sh -v2 build_v2
cp "build_v2/${EXEC_NAME}" ./exec_v2

#############################################
# [MOD] Helper: extract "count: <num>" from output file
#############################################
extract_count() {
    local file="$1"
    # 取最后一次出现的 count，并提取数字；兼容 "count: 655" / "count : 655"
    awk '
        {
            for (i=1; i<=NF; i++) {
                if ($i ~ /^count:$/ || ($i ~ /^count:.*$/ && $i != "count:")) {
                    # 处理两种： "count:" "655" 或 "count:655"
                    if ($i == "count:") {
                        if ((i+1) <= NF && $(i+1) ~ /^[0-9]+$/) last=$(i+1)
                    } else if ($i ~ /^count:[0-9]+$/) {
                        sub(/^count:/,"",$i); last=$i
                    }
                } else if ($i == "count" && (i+1)<=NF && $(i+1) ~ /^:$/ && (i+2)<=NF && $(i+2) ~ /^[0-9]+$/) {
                    # 兼容 "count : 655"
                    last=$(i+2)
                }
            }
        }
        END { if (last!="") print last; else exit 1 }
    ' "$file"
}

extract_total_time() {
    local file="$1"
    # 提取 "Total Time: 29.1100 ms" 行中的数字部分，确保只获取最后的 total time
    awk '
        /Total Time:/ {
            split($3, a, "ms"); print a[1];
        }
    ' "$file"
}

#############################################
# RUN BOTH VERSIONS FOR ALL QUERY GRAPHS
#############################################
echo "===== Running all query graphs ====="

shopt -s nullglob

# [MOD] 统计信息
TOTAL=0
OK=0
DIFF=0
MISSING=0

for folder in ${DATA_DIR}/*/; do
    dirname=$(basename "$folder")
    gfile="${folder}/graph_g.txt"
    qdir="${folder}/query_graph"

    if [ ! -f "$gfile" ]; then
        echo "[Skip] $dirname: graph_g.txt not found"
        continue
    fi

    if [ ! -d "$qdir" ]; then
        echo "[Skip] $dirname: query_graph/ not found"
        continue
    fi

    echo "Dataset: $dirname"

    mkdir -p "${RESULT_DIR}/${dirname}/v1"
    mkdir -p "${RESULT_DIR}/${dirname}/v2"

    for qfile in ${qdir}/*.txt; do
        qname=$(basename "$qfile" .txt)
        echo "  Query: $qname"

        for t in {0..2}; do
            out1="${RESULT_DIR}/${dirname}/v1/${qname}_t=${t}.txt"
            out2="${RESULT_DIR}/${dirname}/v2/${qname}_t=${t}.txt"

            # 先运行并保存输出
            ./exec_v1 -d "$gfile" -q "$qfile" -t $t > "$out1"
            ./exec_v2 -d "$gfile" -q "$qfile" -t $t > "$out2"

            # [MOD] 比较 count 字段（通过 extract_count）
            TOTAL=$((TOTAL+1))

            c1=""
            c2=""
            if c1=$(extract_count "$out1" 2>/dev/null) && c2=$(extract_count "$out2" 2>/dev/null); then
                if [ "$c1" = "$c2" ]; then
                    OK=$((OK+1))
                    echo "    [OK] t=$t count: $c1"
                else
                    DIFF=$((DIFF+1))
                    echo "    [DIFF] t=$t v1_count=$c1 v2_count=$c2"
                    echo "           v1_file=$out1"
                    echo "           v2_file=$out2"
                fi
            else
                MISSING=$((MISSING+1))
                echo "    [MISSING] t=$t count field not found in one/both outputs"
                echo "              v1_file=$out1"
                echo "              v2_file=$out2"
            fi
        done
    done
done

echo "===== ALL DONE ====="
echo "Results saved at: ${RESULT_DIR}"

# [MOD] 汇总
echo "===== SUMMARY ====="
echo "Total runs:   ${TOTAL}"
echo "Matched:      ${OK}"
echo "Different:    ${DIFF}"
echo "Count missing:${MISSING}"

# [MOD] 如果存在不一致，脚本以非 0 退出，方便 CI/批处理发现问题
if [ "${DIFF}" -ne 0 ] || [ "${MISSING}" -ne 0 ]; then
    exit 2
fi
