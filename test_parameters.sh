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

if [[ "$TEST_MODE" != *"S"* ]]; then
    /usr/bin/time -f "$TIME_FORMAT_STRING" -o "$TEMP_DIR"/"$ALG"_time.txt \
        ./kmercamel compute -u -k "$K" -a "$ALG" -O "$MODE" "$INPUT" > "$TEMP_DIR"/"$ALG".txt
    COMMAND="./kmercamel compute -k "$K" -u "$INPUT" > "$TEMP_DIR"/gg_raw.txt && ./kmercamel maskopt -u -t minrun -k "$K" "$TEMP_DIR"/gg_raw.txt > "$TEMP_DIR"/gg.txt"
    /usr/bin/time -f "$TIME_FORMAT_STRING" -o "$TEMP_DIR"/gg_time.txt /bin/sh -c "$COMMAND"

    if [[ "$TEST_MODE" == *"C"* ]]; then
        /usr/bin/time -f "$TIME_FORMAT_STRING" -o "$TEMP_DIR"/"$ALG"_c_time.txt \
            ./kmercamel compute -k "$K" -a "$ALG" -O "$MODE" "$INPUT" > "$TEMP_DIR"/"$ALG"_c.txt
        COMMAND="./kmercamel compute -k "$K" "$INPUT" > "$TEMP_DIR"/gg_c_raw.txt && ./kmercamel maskopt -t minrun -k "$K" "$TEMP_DIR"/gg_c_raw.txt > "$TEMP_DIR"/gg_c.txt"
        /usr/bin/time -f "$TIME_FORMAT_STRING" -o "$TEMP_DIR"/gg_c_time.txt /bin/sh -c "$COMMAND"
    fi
fi

if [[ "$TEST_MODE" != *"F"* ]]; then
    echo =================
fi

if [[ "$TEST_MODE" == *"N"* ]]; then
    O1="$(./scripts/count_noncomplement_kmers.py "$TEMP_DIR"/"$ALG".txt -k "$K")"
    O2="$(./scripts/count_noncomplement_kmers.py "$TEMP_DIR"/gg.txt -k "$K")"
    O3=""
    O4=""
    if [[ "$TEST_MODE" == *"C"* ]]; then
        O3="$(./scripts/count_kmers.py "$TEMP_DIR"/"$ALG"_c.txt -k "$K")"
        O4="$(./scripts/count_kmers.py "$TEMP_DIR"/gg_c.txt -k "$K")"
    fi
fi

L1="$(cat "$TEMP_DIR"/"$ALG".txt | tail -n 1 | wc -m)"
L2="$(cat "$TEMP_DIR"/gg.txt | tail -n 1 | wc -m)"
L3=""
L4=""
if [[ "$TEST_MODE" == *"C"* ]]; then
    L3="$(cat "$TEMP_DIR"/"$ALG"_c.txt | tail -n 1 | wc -m)"
    L4="$(cat "$TEMP_DIR"/gg_c.txt | tail -n 1 | wc -m)"
fi

R1="$(cat "$TEMP_DIR"/"$ALG".txt | tail -n 1 | tr [a-z] '0' | tr -s '0' | tr -d [A-Z] | wc -m)"
R2="$(cat "$TEMP_DIR"/gg.txt | tail -n 1 | tr [a-z] '0' | tr -s '0' | tr -d [A-Z] | wc -m)"
R3=""
R4=""
if [[ "$TEST_MODE" == *"C"* ]]; then      
    R3="$(cat "$TEMP_DIR"/"$ALG"_c.txt | tail -n 1 | tr [a-z] '0' | tr -s '0' | tr -d [A-Z] | wc -m)"
    R4="$(cat "$TEMP_DIR"/gg_c.txt | tail -n 1 | tr [a-z] '0' | tr -s '0' | tr -d [A-Z] | wc -m)"
fi

if [[ "$TEST_MODE" == *"F"* ]]; then
    printf "%s %s" "$INPUT" "$K"
    if [[ "$TEST_MODE" == *"C"* ]]; then
        printf " C"
    else
        printf " -"
    fi
    printf " : "

    printf "%d %d %s;" "$L1" "$R1" "$(cat "$TEMP_DIR"/"$ALG"_time.txt)"
    printf "%d %d %s;" "$L2" "$R2" "$(cat "$TEMP_DIR"/gg_time.txt)"
    if [[ "$TEST_MODE" == *"C"* ]]; then
        printf "%d %d %s;" "$L3" "$R3" "$(cat "$TEMP_DIR"/"$ALG"_c_time.txt)"
        printf "%d %d %s;" "$L4" "$R4" "$(cat "$TEMP_DIR"/gg_c_time.txt)"
    fi
    printf "\n"
else
    if [[ "$TEST_MODE" == *"N"* ]]; then
        echo kmers:
        printf "% 10d\t"$ALG"\t% 10d (c)\n"   "$O1" "$O3"
        printf "% 10d\tgg\t% 10d (c)\n"     "$O2" "$O4"
    fi

    echo lengths:
    printf "% 10d\t"$ALG"\t% 10d (c)\n"   "$L1" "$L3"
    printf "% 10d\tgg\t% 10d (c)\n"     "$L2" "$L4"

    echo runs:
    printf "% 10d\t"$ALG"\t% 10d (c)\n"   "$R1" "$R3"
    printf "% 10d\tgg\t% 10d (c)\n"     "$R2" "$R4"

    echo resources - "$ALG":
    echo "$(cat "$TEMP_DIR"/"$ALG"_time.txt)"
    echo resources - gg:
    echo "$(cat "$TEMP_DIR"/gg_time.txt)"

    if [[ "$TEST_MODE" == *"C"* ]]; then
        echo "resources - "$ALG" (c)":
        echo "$(cat "$TEMP_DIR"/"$ALG"_c_time.txt)"
        echo "resources - gg (c)":
        echo "$(cat "$TEMP_DIR"/gg_c_time.txt)"
    fi
fi
