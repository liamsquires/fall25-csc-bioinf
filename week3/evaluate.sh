#!/bin/bash

echo "Language    Runtime"
echo "-------------------"

# Change to week3 directory if we're not already there
if [ ! -f "test_phylo.py" ]; then
    cd week3
fi

# Run Python tests and capture timing
python_output=$(python test_phylo.py 2>&1)
python_time=$(echo "$python_output" | grep -o '[0-9]\+ms' | head -1)
if [ -z "$python_time" ]; then
    python_time="ERROR"
    echo "Python error output: $python_output" >&2
fi
echo "python      $python_time"

# Run Codon tests and capture timing
# Try to find the correct libpython path
if [ -n "$Python3_ROOT_DIR" ]; then
    # GitHub Actions environment
    export CODON_PYTHON="$Python3_ROOT_DIR/lib/libpython3.13.so"
elif [ -f "/usr/lib/x86_64-linux-gnu/libpython3.12.so" ]; then
    # Local environment
    export CODON_PYTHON="/usr/lib/x86_64-linux-gnu/libpython3.12.so"
else
    # Try to find any libpython
    export CODON_PYTHON=$(find /usr/lib* /opt/* -name "libpython*.so*" 2>/dev/null | head -1)
fi

# echo "Using CODON_PYTHON: $CODON_PYTHON" >&2
codon_output=$(codon run test_phylo.py 2>&1)
codon_time=$(echo "$codon_output" | grep -o '[0-9]\+ms' | head -1)
if [ -z "$codon_time" ]; then
    codon_time="ERROR"
    echo "Codon error output: $codon_output" >&2
fi
echo "codon       $codon_time"