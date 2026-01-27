#!/bin/bash
# Simple compilation script for AdaptiveBinning tool
# Usage: ./compile_adaptive_binning.sh

echo "========================================"
echo "Compiling Adaptive Binning Tool"
echo "========================================"

# Check if ROOT is available
if ! command -v root-config &> /dev/null; then
    echo "ERROR: root-config not found in PATH"
    echo "Please make sure ROOT is installed and sourced"
    echo "Example: source /path/to/root/bin/thisroot.sh"
    exit 1
fi

echo "ROOT found: $(which root)"
echo "ROOT version: $(root-config --version)"

# Get ROOT compilation flags
ROOT_CFLAGS=$(root-config --cflags)
ROOT_LIBS=$(root-config --glibs)

# Compile AdaptiveBinning
echo ""
echo "Compiling AdaptiveBinning.cpp..."
g++ AdaptiveBinning.cpp -o AdaptiveBinning \
    -I../include \
    $ROOT_CFLAGS \
    $ROOT_LIBS \
    -std=c++20

if [ $? -eq 0 ]; then
    echo "✓ Successfully compiled AdaptiveBinning"
    echo "  Executable: ./AdaptiveBinning"
else
    echo "✗ Compilation failed for AdaptiveBinning"
    exit 1
fi

# Compile Example_UseBinning
echo ""
echo "Compiling Example_UseBinning.cpp..."
g++ Example_UseBinning.cpp -o Example_UseBinning \
    -I../include \
    $ROOT_CFLAGS \
    $ROOT_LIBS \
    -std=c++20

if [ $? -eq 0 ]; then
    echo "✓ Successfully compiled Example_UseBinning"
    echo "  Executable: ./Example_UseBinning"
else
    echo "✗ Compilation failed for Example_UseBinning"
    exit 1
fi

echo ""
echo "========================================"
echo "Compilation Complete!"
echo "========================================"
echo ""
echo "Run the adaptive binning tool:"
echo "  ./AdaptiveBinning"
echo ""
echo "Run the example (after generating binning):"
echo "  ./Example_UseBinning binning_scheme.txt"
echo ""
