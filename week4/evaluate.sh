#!/bin/bash
set -e

# Change to week4 directory if we're not already there
if [ ! -f "global.py" ]; then
    cd week4
fi

# Optional: compile Codon versions
for script in global local fitting; do
    echo "Compiling ${script}.py..."
    codon build -release -o ${script}_codon ${script}.py
done

measure_time() {
    start=$(date +%s%3N)
    "$@" >/dev/null 2>&1
    end=$(date +%s%3N)
    echo $((end - start))
}

parse_fasta() {
    awk '/^>/{if(seq){print seq; seq=""} next} {seq=seq$0} END{print seq}' "$1"
}

queries=($(parse_fasta q1.fa))
targets=($(parse_fasta t1.fa))

declare -A results

printf "%-20s %-10s %-10s\n" "Method" "Language" "Runtime"
printf "%-20s\n" "--------------------------------------"

for i in "${!queries[@]}"; do
    q="${queries[$i]}"
    t="${targets[$i]}"
    qi="q$((i+1))"

    for name in global local fitting; do
        
        results["${name}-${qi}"]=$(measure_time python3 ${name}.py "$t" "$q")
        printf "%-20s %-10s %sms\n" "${name}-${qi}" "python" "${results[${name}-${qi}]}"
        sleep 0.1
    done
done

human_seq=$(awk '!/^>/{printf "%s", $0}' MT-human.fa)
orang_seq=$(awk '!/^>/{printf "%s", $0}' MT-orang.fa)

for name in global local fitting; do
    
    results["${name}-mt_human"]=$(measure_time python3 ${name}.py "$human_seq" "$orang_seq")
    printf "%-20s %-10s %sms\n" "${name}-mt_human" "python" "${results[${name}-mt_human]}"
    sleep 0.1
done

results=()

for i in "${!queries[@]}"; do
    q="${queries[$i]}"
    t="${targets[$i]}"
    qi="q$((i+1))"

    for name in global local fitting; do
        
        results["${name}-${qi}"]=$(measure_time ./${name}_codon "$t" "$q")
        printf "%-20s %-10s %sms\n" "${name}-${qi}" "codon" "${results[${name}-${qi}]}"
        sleep 0.1
    done
done

for name in global local fitting; do
    
    results["${name}-mt_human"]=$(measure_time ./${name}_codon "$human_seq" "$orang_seq")
    printf "%-20s %-10s %sms\n" "${name}-mt_human" "codon" "${results[${name}-mt_human]}"
    sleep 0.1
done
