#!/bin/bash

# Classification Workflow Testing Script
# This script provides easy access to common testing scenarios

set -e  # Exit on any error

# Color output for better readability
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}Classification Workflow Testing${NC}"
echo "================================="

# Check if Python is available
if ! command -v python3 &> /dev/null; then
    echo -e "${RED}Error: Python3 is required but not found${NC}"
    exit 1
fi

# Function to run a test with parameters
run_test() {
    local test_name="$1"
    local output_dir="$2"
    local max_seqs="$3"
    local curated_flag="$4"
    
    echo -e "${YELLOW}Starting $test_name...${NC}"
    
    cmd="python3 test_classification_workflow.py --output-dir $output_dir"
    
    if [ ! -z "$max_seqs" ]; then
        cmd="$cmd --max-sequences $max_seqs"
    fi
    
    if [ "$curated_flag" == "true" ]; then
        cmd="$cmd --curated-only"
    fi
    
    echo "Running: $cmd"
    
    if $cmd; then
        echo -e "${GREEN}✓ $test_name completed successfully${NC}"
        echo -e "Results saved to: ${BLUE}$output_dir${NC}"
    else
        echo -e "${RED}✗ $test_name failed${NC}"
        return 1
    fi
}

# Parse command line arguments
case "${1:-menu}" in
    "quick"|"q")
        echo "Running quick test (10 sequences)..."
        run_test "Quick Test" "test_results_quick" "10" "false"
        ;;
    
    "curated"|"c")
        echo "Running curated test (50 sequences)..."
        run_test "Curated Test" "test_results_curated" "50" "true"
        ;;
    
    "small"|"s")
        echo "Running small test (100 sequences)..."
        run_test "Small Test" "test_results_small" "100" "false"
        ;;
        
    "medium"|"m")
        echo "Running medium test (500 sequences)..."
        run_test "Medium Test" "test_results_medium" "500" "false"
        ;;
    
    "full"|"f")
        echo -e "${RED}WARNING: This will test ALL sequences in the database!${NC}"
        echo "This could take many hours or days depending on your database size."
        read -p "Are you sure you want to continue? (yes/no): " -r
        if [[ $REPLY =~ ^[Yy][Ee][Ss]$ ]]; then
            run_test "Full Test" "test_results_full" "" "false"
        else
            echo "Test cancelled."
        fi
        ;;
    
    "interactive"|"i")
        echo "Starting interactive mode..."
        python3 run_classification_tests.py
        ;;
    
    "analyze"|"a")
        if [ -z "$2" ]; then
            echo "Usage: $0 analyze <results_directory>"
            echo "Example: $0 analyze test_results_quick"
            exit 1
        fi
        
        if [ ! -d "$2" ]; then
            echo -e "${RED}Error: Directory '$2' does not exist${NC}"
            exit 1
        fi
        
        echo "Analyzing results in $2..."
        python3 -c "
import sys
sys.path.append('.')
from run_classification_tests import analyze_previous_results
analyze_previous_results('$2')
"
        ;;
    
    "help"|"h"|"--help")
        echo "Usage: $0 [command]"
        echo ""
        echo "Commands:"
        echo "  quick, q          Run quick test (10 sequences)"
        echo "  curated, c        Run curated test (50 sequences)"
        echo "  small, s          Run small test (100 sequences)"
        echo "  medium, m         Run medium test (500 sequences)"
        echo "  full, f           Run full database test (all sequences)"
        echo "  interactive, i    Start interactive mode"
        echo "  analyze, a <dir>  Analyze previous results"
        echo "  help, h           Show this help message"
        echo ""
        echo "Examples:"
        echo "  $0 quick              # Quick test"
        echo "  $0 curated            # Test curated sequences"
        echo "  $0 analyze test_results_quick  # Analyze previous results"
        echo "  $0                    # Show menu (default)"
        ;;
    
    "menu"|*)
        echo ""
        echo "Available test options:"
        echo "1. Quick test (10 sequences) - 5-10 minutes"
        echo "2. Curated test (50 curated sequences) - 30-60 minutes"
        echo "3. Small test (100 sequences) - 1-2 hours"
        echo "4. Medium test (500 sequences) - 4-8 hours"
        echo "5. Full database test - Many hours/days"
        echo "6. Interactive mode"
        echo "7. Analyze previous results"
        echo ""
        read -p "Select an option (1-7) or 'h' for help: " -r choice
        
        case $choice in
            1) ./run_tests.sh quick ;;
            2) ./run_tests.sh curated ;;
            3) ./run_tests.sh small ;;
            4) ./run_tests.sh medium ;;
            5) ./run_tests.sh full ;;
            6) ./run_tests.sh interactive ;;
            7) 
                read -p "Enter results directory: " -r results_dir
                ./run_tests.sh analyze "$results_dir"
                ;;
            h|H) ./run_tests.sh help ;;
            *) 
                echo -e "${RED}Invalid option${NC}"
                exit 1
                ;;
        esac
        ;;
esac

echo -e "${GREEN}Done!${NC}" 