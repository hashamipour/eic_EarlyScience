#!/bin/bash

# Interactive Binning Scheme Plotter - Build and Run Script
# Usage: ./run_interactive_binning.sh [input_file.root]

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}Interactive Binning Scheme Plotter${NC}"
echo -e "${GREEN}========================================${NC}"

# Check if input file is provided
if [ $# -eq 0 ]; then
    echo -e "${YELLOW}Usage: $0 <input_file.root>${NC}"
    echo -e "${YELLOW}Example: $0 ../DDIS_Combined_output.root${NC}"
    exit 1
fi

INPUT_FILE="$1"

# Check if input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo -e "${RED}Error: Input file '$INPUT_FILE' not found${NC}"
    exit 1
fi

echo -e "${GREEN}Input file: $INPUT_FILE${NC}"

# Get the directory of this script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_ROOT="$SCRIPT_DIR/.."

# Check if we should use CMake build or direct compilation
EXECUTABLE=""

if [ -f "$PROJECT_ROOT/build/ddis_plot_binning_interactive" ]; then
    echo -e "${GREEN}Using CMake-built executable${NC}"
    EXECUTABLE="$PROJECT_ROOT/build/ddis_plot_binning_interactive"
elif [ -f "$SCRIPT_DIR/Plot_BinningScheme_Interactive" ]; then
    echo -e "${GREEN}Using locally compiled executable${NC}"
    EXECUTABLE="$SCRIPT_DIR/Plot_BinningScheme_Interactive"
else
    echo -e "${YELLOW}Executable not found. Attempting to compile...${NC}"

    # Try to compile directly
    cd "$SCRIPT_DIR"

    if command -v root-config &> /dev/null; then
        echo -e "${GREEN}Compiling Plot_BinningScheme_Interactive.cpp...${NC}"
        g++ -std=c++17 -o Plot_BinningScheme_Interactive \
            Plot_BinningScheme_Interactive.cpp \
            `root-config --cflags --libs`

        if [ $? -eq 0 ]; then
            echo -e "${GREEN}Compilation successful!${NC}"
            EXECUTABLE="$SCRIPT_DIR/Plot_BinningScheme_Interactive"
        else
            echo -e "${RED}Compilation failed${NC}"
            exit 1
        fi
    else
        echo -e "${RED}Error: ROOT is not installed or root-config not in PATH${NC}"
        echo -e "${YELLOW}Please install ROOT or use CMake to build:${NC}"
        echo -e "${YELLOW}  cd $PROJECT_ROOT/build${NC}"
        echo -e "${YELLOW}  cmake ..${NC}"
        echo -e "${YELLOW}  make ddis_plot_binning_interactive${NC}"
        exit 1
    fi
fi

# Create figs directory if it doesn't exist
mkdir -p "$PROJECT_ROOT/figs"

echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}Starting interactive session...${NC}"
echo -e "${GREEN}========================================${NC}"
echo ""

# Run the interactive tool
cd "$PROJECT_ROOT"
"$EXECUTABLE" "$INPUT_FILE"

EXIT_CODE=$?

if [ $EXIT_CODE -eq 0 ]; then
    echo ""
    echo -e "${GREEN}========================================${NC}"
    echo -e "${GREEN}Session completed successfully!${NC}"
    echo -e "${GREEN}========================================${NC}"
    echo -e "${GREEN}Check the figs/ directory for output plots${NC}"
else
    echo -e "${RED}Error: Program exited with code $EXIT_CODE${NC}"
fi

exit $EXIT_CODE
