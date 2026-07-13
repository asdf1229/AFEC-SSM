#!/usr/bin/env bash

BUILD_DIR="build"
GENERATOR=""
SOURCE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

for arg in "$@"; do
    case "$arg" in
        -h|--help)
            echo "Usage: ./build.sh [build_directory]"
            echo
            echo "Description:"
            echo "  Configure and build all algorithm executables."
            echo "  Default build directory is './build'."
            exit 0
            ;;
        *)
            if [[ ! "$arg" =~ ^- ]]; then
                BUILD_DIR="$arg"
            fi
            ;;
    esac
done

if command -v ninja >/dev/null 2>&1; then
    GENERATOR="Ninja"
elif command -v mingw32-make >/dev/null 2>&1; then
    GENERATOR="MinGW Makefiles"
elif command -v make >/dev/null 2>&1; then
    GENERATOR="Unix Makefiles"
fi

echo "========================================"
echo "Build Configuration:"
echo "  Directory : $BUILD_DIR"
echo "  Build type: Release"
if [ -n "$GENERATOR" ]; then
    echo "  Generator : $GENERATOR"
else
    echo "  Generator : default"
fi
echo "========================================"

if [[ "$BUILD_DIR" = /* ]]; then
    BUILD_PATH="$BUILD_DIR"
else
    BUILD_PATH="${SOURCE_DIR}/${BUILD_DIR#./}"
fi

if [ -d "$BUILD_PATH" ]; then
    echo "Cleaning build directory '$BUILD_DIR'..."
    rm -rf "$BUILD_PATH"
fi

mkdir -p "$BUILD_PATH"
cd "$BUILD_PATH" || exit 1

echo "Running CMake..."
if [ -n "$GENERATOR" ]; then
    cmake -G "$GENERATOR" -DCMAKE_BUILD_TYPE=Release "$SOURCE_DIR" || exit 1
else
    cmake -DCMAKE_BUILD_TYPE=Release "$SOURCE_DIR" || exit 1
fi

echo "Compiling..."
cmake --build . --parallel || exit 1

shopt -s nullglob
executables=()
for candidate in ssm_ged_*; do
    if [ -f "$candidate" ]; then
        executables+=("$candidate")
    fi
done

if [ ${#executables[@]} -eq 0 ]; then
    echo "Error: no algorithm executables were produced."
    exit 1
fi

echo "========================================"
echo "Build successful!"
echo "Executables:"
for exe in "${executables[@]}"; do
    echo "  $BUILD_DIR/$exe"
done
echo "========================================"
