#!/bin/bash
set -e

# ===================== 配置 =====================
DATA_DIR=./test/dataset
EXEC_NAME=ssm_ged

# 选择要测试的版本（v1 或 v2），默认为 v1
VERSION=${1:-v2}

# 如果没有传入有效的版本，则退出
if [[ "$VERSION" != "v1" && "$VERSION" != "v2" ]]; then
    echo "Error: Please specify a valid version (v1 or v2)"
    exit 1
fi
# ===============================================

#############################################
# CREATE TIMESTAMPED RESULT DIRECTORY
#############################################
timestamp=$(date +"%Y%m%d_%H%M%S")
RESULT_DIR="./result/${timestamp}"
mkdir -p "${RESULT_DIR}"

echo "===== RESULTS WILL BE SAVED TO: ${RESULT_DIR} ====="

#############################################
# BUILD SELECTED VERSION (v1 or v2)
#############################################
if [ "$VERSION" == "v1" ]; then
    echo "===== Building Version 1 (Approximate_Matching) ====="
    ./build.sh build_v1
    EXEC_PATH="build_v1/${EXEC_NAME}"
else
    echo "===== Building Version 2 (Approximate_Matching_V2) ====="
    ./build.sh -v2 build_v2
    EXEC_PATH="build_v2/${EXEC_NAME}"
fi

#############################################
# [MOD] Helper: extract "Total Time" from output file
#############################################
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
# RUN SELECTED VERSION FOR ALL QUERY GRAPHS
#############################################
echo "===== Running version $VERSION for all query graphs ====="

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

    mkdir -p "${RESULT_DIR}/${dirname}/$VERSION"

    for qfile in ${qdir}/*.txt; do
        qname=$(basename "$qfile" .txt)
        echo "  Query: $qname"

        for t in {0..3}; do
            out="${RESULT_DIR}/${dirname}/${VERSION}/${qname}_t=${t}.txt"

            # 先运行并保存输出
            $EXEC_PATH -d "$gfile" -q "$qfile" -t $t > "$out"

            # [MOD] 比较 total time 字段
            TOTAL=$((TOTAL+1))

            t1=""
            if t1=$(extract_total_time "$out" 2>/dev/null); then
                echo "    [OK] t=$t Total Time: $t1 ms"
            else
                MISSING=$((MISSING+1))
                echo "    [MISSING] t=$t Total Time field not found"
                echo "              output_file=$out"
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
echo "Missing Total Time:${MISSING}"

# [MOD] 如果存在不一致，脚本以非 0 退出，方便 CI/批处理发现问题
if [ "${MISSING}" -ne 0 ]; then
    exit 2
fi
