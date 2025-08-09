#!/bin/bash

set -ueo pipefail

ALG="joint"
MODE="${1:-""}"

if [[ "$MODE" == "-h" ]] || [[ "$MODE" == "--help" ]] || [[ -z $MODE ]]; then
    echo "Usage: mode (runs / zeros), input file, k, test mode parameters (C - complements, N - counting check (slow))"
    exit 0
fi

INPUT="$2"
K="$3"
TEST_MODE="${4:-""}"

# TEST_MODE:
# "N" for counting kmers check - with python script, slow
# "C" for also running the computation with complements
# "F" for output in one line without explanatory texts
# "L" for using local directory instead of temp one
# "S" for only checking previously computed files (only when L was used before)

if [[ "$TEST_MODE" == *"L"* ]] || [[ "$TEST_MODE" == *"S"* ]]; then
    TEMP_DIR="testing_outputs"
    mkdir -p "$TEMP_DIR"
else
    TEMP_DIR="$(mktemp -d)"
    trap '{ rm -rf -- "$TEMP_DIR"; }' EXIT
fi

TIME_FORMAT_STRING="time:\t%U\nmemory:\t%M"

if [[ "$TEST_MODE" == *"F"* ]]; then
    TIME_FORMAT_STRING="%U %M"
fi

if [[ "$MODE" == "runs" ]]; then
    MASK_OPT="minrun"
else
    MASK_OPT="maxone"
fi

if [[ "$TEST_MODE" != *"S"* ]]; then
    if [[ "$TEST_MODE" == *"C"* ]]; then
        /usr/bin/time -f "$TIME_FORMAT_STRING" -o "$TEMP_DIR"/"$ALG"_time.txt \
            ./kmercamel compute -k "$K" -a "$ALG" -O "$MODE" "$INPUT" > "$TEMP_DIR"/"$ALG".txt
        COMMAND="./kmercamel compute -k "$K" "$INPUT" > "$TEMP_DIR"/gg_raw.txt && ./kmercamel maskopt -t "$MASK_OPT" -k "$K" "$TEMP_DIR"/gg_raw.txt > "$TEMP_DIR"/gg.txt"
        /usr/bin/time -f "$TIME_FORMAT_STRING" -o "$TEMP_DIR"/gg_time.txt /bin/sh -c "$COMMAND"
    else
        /usr/bin/time -f "$TIME_FORMAT_STRING" -o "$TEMP_DIR"/"$ALG"_time.txt \
            ./kmercamel compute -u -k "$K" -a "$ALG" -O "$MODE" "$INPUT" > "$TEMP_DIR"/"$ALG".txt
        COMMAND="./kmercamel compute -k "$K" -u "$INPUT" > "$TEMP_DIR"/gg_raw.txt && ./kmercamel maskopt -u -t "$MASK_OPT" -k "$K" "$TEMP_DIR"/gg_raw.txt > "$TEMP_DIR"/gg.txt"
        /usr/bin/time -f "$TIME_FORMAT_STRING" -o "$TEMP_DIR"/gg_time.txt /bin/sh -c "$COMMAND"
    fi
fi

if [[ "$TEST_MODE" != *"F"* ]]; then
    echo =================
fi

if [[ "$TEST_MODE" == *"N"* ]]; then
    if [[ "$TEST_MODE" == *"C"* ]]; then
        O1="$(./scripts/count_kmers.py "$TEMP_DIR"/"$ALG".txt -k "$K")"
        O2="$(./scripts/count_kmers.py "$TEMP_DIR"/gg.txt -k "$K")"
    else
        O1="$(./scripts/count_noncomplement_kmers.py "$TEMP_DIR"/"$ALG".txt -k "$K")"
        O2="$(./scripts/count_noncomplement_kmers.py "$TEMP_DIR"/gg.txt -k "$K")"
    fi
fi

if [[ "$TEST_MODE" == *"C"* ]]; then
    L1="$(cat "$TEMP_DIR"/"$ALG".txt | tail -n 1 | wc -m)"
    L2="$(cat "$TEMP_DIR"/gg.txt | tail -n 1 | wc -m)"
else
    L1="$(cat "$TEMP_DIR"/"$ALG".txt | tail -n 1 | wc -m)"
    L2="$(cat "$TEMP_DIR"/gg.txt | tail -n 1 | wc -m)"
fi

if [[ "$TEST_MODE" == *"C"* ]]; then      
    R1="$(cat "$TEMP_DIR"/"$ALG".txt | tail -n 1 | tr [a-z] '0' | tr -s '0' | tr -d [A-Z] | wc -m)"
    R2="$(cat "$TEMP_DIR"/gg.txt | tail -n 1 | tr [a-z] '0' | tr -s '0' | tr -d [A-Z] | wc -m)"
else
    R1="$(cat "$TEMP_DIR"/"$ALG".txt | tail -n 1 | tr [a-z] '0' | tr -s '0' | tr -d [A-Z] | wc -m)"
    R2="$(cat "$TEMP_DIR"/gg.txt | tail -n 1 | tr [a-z] '0' | tr -s '0' | tr -d [A-Z] | wc -m)"
fi

if [[ "$TEST_MODE" == *"C"* ]]; then      
    Z1="$(cat "$TEMP_DIR"/"$ALG".txt | tail -n 1 | tr -d [A-Z] | wc -m)"
    Z2="$(cat "$TEMP_DIR"/gg.txt | tail -n 1 | tr -d [A-Z] | wc -m)"
else
    Z1="$(cat "$TEMP_DIR"/"$ALG".txt | tail -n 1 | tr -d [A-Z] | wc -m)"
    Z2="$(cat "$TEMP_DIR"/gg.txt | tail -n 1 | tr -d [A-Z] | wc -m)"
fi

if [[ "$TEST_MODE" == *"F"* ]]; then
    printf "%s %s" "$INPUT" "$K"
    if [[ "$TEST_MODE" == *"C"* ]]; then
        printf " C"
    else
        printf " -"
    fi
    printf " : "

    printf "%d %d %s;" "$L1" "$R1" "$Z1" "$(cat "$TEMP_DIR"/"$ALG"_time.txt)"
    printf "%d %d %s;" "$L2" "$R2" "$Z2" "$(cat "$TEMP_DIR"/gg_time.txt)"
    printf "\n"
else
    if [[ "$TEST_MODE" == *"N"* ]]; then
        echo kmers:
        printf "% 10d\t"$ALG"\n"    "$O1"
        printf "% 10d\tgg\n"        "$O2"
    fi

    echo lengths:
    printf "% 10d\t"$ALG"\n"        "$L1"
    printf "% 10d\tgg\n"            "$L2"

    echo runs:
    printf "% 10d\t"$ALG"\n"        "$R1"
    printf "% 10d\tgg\n"            "$R2"

    echo zeros:
    printf "% 10d\t"$ALG"\n"        "$Z1"
    printf "% 10d\tgg\n"            "$Z2"

    echo resources - "$ALG":
    echo "$(cat "$TEMP_DIR"/"$ALG"_time.txt)"
    echo resources - gg:
    echo "$(cat "$TEMP_DIR"/gg_time.txt)"
fi
