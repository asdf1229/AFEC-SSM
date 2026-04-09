#!/bin/bash

BUILD_DIR="build"
GENERATOR=""

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
if [ -n "$GENERATOR" ]; then
    echo "  Generator : $GENERATOR"
else
    echo "  Generator : default"
fi
echo "========================================"

if [ -d "$BUILD_DIR" ]; then
    echo "Cleaning build directory '$BUILD_DIR'..."
    rm -rf "$BUILD_DIR"
fi

mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR" || exit 1

echo "Running CMake..."
if [ -n "$GENERATOR" ]; then
    cmake -G "$GENERATOR" ..
else
    cmake ..
fi

echo "Compiling..."
cmake --build . --parallel

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
