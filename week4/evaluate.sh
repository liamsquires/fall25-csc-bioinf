#!/bin/bash
set -e

# Change to week4 directory if we're not already there
if [ ! -f "global.py" ]; then
    cd week4
fi

measure_time() {
    start=$(date +%s%3N)
    "$@"
    end=$(date +%s%3N)
    echo $((end - start))
}

parse_fasta() {
    awk '/^>/{if(seq){print seq; seq=""} next} {seq=seq$0} END{print seq}' "$1"
}

queries=($(parse_fasta q1.fa))
targets=($(parse_fasta t1.fa))

declare -A results

# q1–q5 vs t1–t5
for i in "${!queries[@]}"; do
    q="${queries[$i]}"
    t="${targets[$i]}"
    qi="q$((i+1))"

    results["global-$qi"]=$(measure_time python3 global.py "$t" "$q")
    results["local-$qi"]=$(measure_time python3 local.py "$t" "$q")
    results["fitting-$qi"]=$(measure_time python3 fitting.py "$t" "$q")
    results["affine-$qi"]=$(measure_time python3 affine.py "$t" "$q")
done

# Human vs Orangutan
human_seq=$(awk '!/^>/{printf "%s", $0}' MT-human.fa)
orang_seq=$(awk '!/^>/{printf "%s", $0}' MT-orang.fa)

results["global-mt_human"]=$(measure_time python3 global.py "$human_seq" "$orang_seq")
results["local-mt_human"]=$(measure_time python3 local.py "$human_seq" "$orang_seq")
results["fitting-mt_human"]=$(measure_time python3 fitting.py "$human_seq" "$orang_seq")
results["affine-mt_human"]=$(measure_time python3 affine.py "$human_seq" "$orang_seq")

printf "%-20s %-10s %-10s\n" "Method" "Language" "Runtime"
printf "%-20s\n" "--------------------------------------"

for name in global local fitting affine; do
    for qi in q1 q2 q3 q4 q5 mt_human; do
        key="${name}-${qi}"
        if [[ -n "${results[$key]}" ]]; then
            printf "%-20s %-10s %sms\n" "$key" "python" "${results[$key]}"
        fi
    done
done
