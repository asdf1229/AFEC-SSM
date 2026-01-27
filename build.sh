#!/bin/bash

# build/
#   ├── CMakeCache.txt
#   ├── Makefile
#   └── ssm_ged

BUILD_DIR="build"
ENABLE_V2="OFF"
CLEAN_BUILD=true
BUILD_DIR_SET_BY_USER=false

# Parse command-line options and arguments
for arg in "$@"
do
    case $arg in
        -h|--help)
            echo "Usage: ./build.sh [build_directory]"
            echo
            echo "Description:"
            echo "  Build the project using CMake and Make."
            echo "  Default build directory is './build'."
            exit 0
            ;;
        *)
            # If the argument does not begin with '-', assume it is the build directory name
            if [[ ! $arg =~ ^- ]]; then
                BUILD_DIR="$arg"
            fi
            ;;
    esac
done

echo "========================================"
echo "Build Configuration:"
echo "  Directory : $BUILD_DIR"
echo "========================================"

# clean
if [ -d "$BUILD_DIR" ]; then
    echo "Cleaning build directory '$BUILD_DIR'..."
    rm -rf "$BUILD_DIR"
fi

mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR" || exit 1

# Run CMake
echo "Running CMake..."
cmake ..

# Run Make
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
