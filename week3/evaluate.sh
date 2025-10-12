#!/bin/bash

echo "Language    Runtime"
echo "-------------------"

# Run Python tests and capture timing
python_output=$(python test_phylo.py 2>&1)
python_time=$(echo "$python_output" | grep -o '[0-9]\+ms' | head -1)
if [ -z "$python_time" ]; then
    python_time="ERROR"
fi
echo "python      $python_time"

# Run Codon tests and capture timing
export CODON_PYTHON=/usr/lib/x86_64-linux-gnu/libpython3.12.so
codon_output=$(codon run test_phylo.py 2>&1)
codon_time=$(echo "$codon_output" | grep -o '[0-9]\+ms' | head -1)
if [ -z "$codon_time" ]; then
    codon_time="ERROR"
fi
echo "codon       $codon_time"