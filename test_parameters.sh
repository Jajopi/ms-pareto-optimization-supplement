#!/bin/bash

set -ueo pipefail

ALG="joint"
MODE="${1:-""}"

if [[ "$MODE" == "-h" ]] || [[ "$MODE" == "--help" ]] || [[ -z $MODE ]]; then
    echo "Usage: mode (runs / zeros), input file, k, test mode parameters (C, N, F, L, S, P - see comments in the source code)"
    exit 0
fi

INPUT="$2"
K="$3"
TEST_MODE="${4:-""}"

# TEST_MODE:
# "C" to run the computation with complements
# "N" ! currently disabled - to perform check by counting kmers - with python script, very slow
# "F" to output stats for one measurement in one line without explanatory texts (useful for collecting data)
# "L" to use local directory instead of temp one for computation (is not deleted after)
# "S" to only check previously computed files (when parameter L was used before)
# "P" to print the computed masked superstrings (can be huge)

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

if [[ "$MODE" == "runs"  ]]; then MASK_OPT="minrun"; fi
if [[ "$MODE" == "zeros" ]]; then MASK_OPT="maxone"; fi
if [[ -z "${MASK_OPT:-""}" ]]; then
    echo Wrong mode "$MODE", must be \"runs\" or \"zeros\".
    exit 1
fi

if [[ "$TEST_MODE" != *"S"* ]]; then
    if [[ "$TEST_MODE" == *"C"* ]]; then COMPL=""; else COMPL="-u"; fi

    COMMAND="./kmercamel compute -k $K $COMPL -a $ALG -O $MODE $INPUT > $TEMP_DIR/$ALG.txt"
    /usr/bin/time -f "$TIME_FORMAT_STRING" -o "$TEMP_DIR"/"$ALG"_time.txt /bin/sh -c "$COMMAND"

    COMMAND="./kmercamel compute -k $K $COMPL $INPUT > $TEMP_DIR/gg_raw.txt && \
            ./kmercamel maskopt -t $MASK_OPT -k $K $COMPL $TEMP_DIR/gg_raw.txt > $TEMP_DIR/gg.txt"
    /usr/bin/time -f "$TIME_FORMAT_STRING" -o "$TEMP_DIR"/gg_time.txt /bin/sh -c "$COMMAND"
fi

# if [[ "$TEST_MODE" == *"N"* ]]; then
#     if [[ "$TEST_MODE" == *"C"* ]]; then
#         O1="$(./scripts/count_kmers.py "$TEMP_DIR"/"$ALG".txt -k "$K")"
#         O2="$(./scripts/count_kmers.py "$TEMP_DIR"/gg.txt -k "$K")"
#     else
#         O1="$(./scripts/count_noncomplement_kmers.py "$TEMP_DIR"/"$ALG".txt -k "$K")"
#         O2="$(./scripts/count_noncomplement_kmers.py "$TEMP_DIR"/gg.txt -k "$K")"
#     fi
# fi

L1="$(echo "$(cat "$TEMP_DIR"/"$ALG".txt    | tail -n 1 | wc -m)"-1 | bc)"
L2="$(echo "$(cat "$TEMP_DIR"/gg.txt        | tail -n 1 | wc -m)"-1 | bc)"

R1="$(echo "$(cat "$TEMP_DIR"/"$ALG".txt    | tail -n 1 | tr [a-z] '0' | tr -s '0' | tr -d [A-Z] | wc -m)"-1 | bc)"
R2="$(echo "$(cat "$TEMP_DIR"/gg.txt        | tail -n 1 | tr [a-z] '0' | tr -s '0' | tr -d [A-Z] | wc -m)"-1 | bc)"

Z1="$(echo "$(cat "$TEMP_DIR"/"$ALG".txt    | tail -n 1 | tr -d [A-Z] | wc -m)"-1 | bc)"
Z2="$(echo "$(cat "$TEMP_DIR"/gg.txt        | tail -n 1 | tr -d [A-Z] | wc -m)"-1 | bc)"

if [[ "$TEST_MODE" != *"F"* ]]; then
    echo ======== RESULTS: =========

    # if [[ "$TEST_MODE" == *"N"* ]]; then
    #     echo kmers:
    #     printf "% 10d\t"$ALG"\n"    "$O1"
    #     printf "% 10d\tgg\n"        "$O2"
    # fi

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
    cat "$TEMP_DIR"/"$ALG"_time.txt
    echo resources - gg:
    cat "$TEMP_DIR"/gg_time.txt

    if [[ "$TEST_MODE" == *"P"* ]]; then
        echo ""
        echo masked superstring - "$ALG":
        tail -n 1 "$TEMP_DIR"/"$ALG".txt
        echo ""
        echo masked superstring - gg:
        tail -n 1 "$TEMP_DIR"/gg.txt
    fi
else
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
fi
