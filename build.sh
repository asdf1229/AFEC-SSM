#!/bin/bash

BUILD_DIR="build"
ENABLE_V2="OFF"
CLEAN_BUILD=true
BUILD_DIR_SET_BY_USER=false

# Parse command-line options and arguments
for arg in "$@"
do
    case $arg in
        -v2|--version2)
        ENABLE_V2="ON"
        ;;
        -h|--help)
        echo "Usage: ./build.sh [options] [build_directory]"
        echo "Options:"
        echo "  -v2, --version2   Enable APPROXIMATE_MATCHING_V2"
        echo "  -h,  --help       Show this help message"
        exit 0
        ;;
        *)
        # If the argument does not begin with '-', assume it is the build directory name
        if [[ ! $arg =~ ^- ]]; then
            BUILD_DIR=$arg
            BUILD_DIR_SET_BY_USER=true
        fi
        ;;
    esac
done

# If v2 enabled and user didn't specify build dir, use build_v2
if [[ "$ENABLE_V2" == "ON" && "$BUILD_DIR_SET_BY_USER" == false ]]; then
    BUILD_DIR="build_v2"
fi

echo "========================================"
echo "Build Configuration:"
echo "  Directory : $BUILD_DIR"
echo "  Version 2 : $ENABLE_V2"
echo "========================================"

# clean
if [ -d "$BUILD_DIR" ]; then
    echo "Cleaning build directory '$BUILD_DIR'..."
    rm -rf "$BUILD_DIR"
fi

mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR" || exit

# run CMake
echo "Running CMake..."
cmake -DAPPROXIMATE_MATCHING_V2=$ENABLE_V2 ..
# Check if CMake configuration was successful
if [ $? -ne 0 ]; then
    echo "Error: CMake configuration failed."
    exit 1
fi

# run Make
echo "Compiling..."
make -j$(nproc)

if [ $? -eq 0 ]; then
    echo "========================================"
    echo "Build successful!"
    echo "Executable: $BUILD_DIR/ssm_ged"
    echo "========================================"
else
    echo "Error: Compilation failed."
    exit 1
fi
