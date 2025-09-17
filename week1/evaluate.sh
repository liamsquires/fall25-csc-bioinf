#!/bin/bash
set -euo pipefail

# Helper: compute N50 given an assembly_output.txt file
compute_n50() {
  awk '{print $2}' "$1" > contig_lengths.txt
  sort -nr contig_lengths.txt -o contig_lengths.txt
  total=$(awk '{s+=$1} END{print s}' contig_lengths.txt)
  if [ "$total" -eq 0 ]; then
    echo 0
    return
  fi
  half=$(( total / 2 ))
  awk -v half=$half '
  {
    sum+=$1
    if(sum >= half) {
      print $1
      exit
    }
  }' contig_lengths.txt
}

# Helper: run one program (python or codon), record runtime and N50
run_eval() {
  dataset=$1
  lang=$2
  exe=$3

  # Special stack size for data4
  if [ "$dataset" = "data4" ]; then
    ulimit -s 8192000
  fi

  start=$(date +%s)

  if [ "$lang" = "python" ]; then
    python3 code/genome-assembly-Zhongyu-Chen/main.py \
      code/genome-assembly-Zhongyu-Chen/"$dataset" > assembly_output.txt
  else
    "$exe" code/genome-assembly-Zhongyu-Chen/"$dataset" > assembly_output.txt
  fi

  end=$(date +%s)
  runtime=$(( end - start ))

  # Format runtime as H:MM:SS
  hours=$(( runtime / 3600 ))
  mins=$(( (runtime % 3600) / 60 ))
  secs=$(( runtime % 60 ))
  formatted=$(printf "%d:%02d:%02d" $hours $mins $secs)

  n50=$(compute_n50 assembly_output.txt)

  echo -e "${dataset}\t${lang}\t\t${formatted}\t${n50}"
}

echo -e "Dataset\tLanguage\tRuntime\tN50"
echo "-------------------------------------------------------------------------------------------------------"

for d in data1 data2 data3 data4; do
  run_eval "$d" "python" "python"
  run_eval "$d" "codon" "code/genome-assembly-codon/assembler"
done
