#!/usr/bin/env bash
# Tests of the pareto optimization (run from the kmercamel-pareto directory, e.g. by `make pareto-test`).
#   BIN      binary to test (default ./kmercamel; use a sanitizer build to also check memory errors)
#   THREADS  number of threads for the multithreaded tests (default 4)
#   REF_BIN  optional reference binary; single-threaded outputs have to be byte-identical to it
set -u
BIN=${BIN:-./kmercamel}
THREADS=${THREADS:-4}
REF_BIN=${REF_BIN:-}
DIR=$(cd "$(dirname "$0")" && pwd)
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
VERIFY="python3 -I $DIR/verify_ms.py"
FAILED=0; PASSED=0

pass(){ PASSED=$((PASSED + 1)); }
fail(){ FAILED=$((FAILED + 1)); echo "FAIL: $*"; }

python3 -I "$DIR/gen_inputs.py" data/spneumoniae.fa "$WORK" || { echo "Cannot generate inputs"; exit 1; }

# Compute a masked superstring, verify it represents the k-mers and that the reported objective matches the recount
check(){ # fasta, options...
    local fasta=$1; shift
    local out="$WORK/out.fa" log="$WORK/log.txt"
    "$BIN" compute -a pareto "$@" -o "$out" "$fasta" > "$log" 2>&1
    local rc=$?
    local k=$(echo "$*" | grep -oE -- '-k [0-9]+' | cut -d' ' -f2)
    local unidirectional=$(echo "$*" | grep -o -- '-u')
    local objective=$(echo "$*" | grep -oE 'runs|zeros')
    local penalty=$(echo "$*" | grep -oE -- '-p [0-9]+' | cut -d' ' -f2)
    [ -z "$penalty" ] && { [ "$objective" = runs ] && penalty=7 || penalty=2; }
    local result
    result=$($VERIFY "$fasta" "$out" "$k" $unidirectional 2>&1)
    if [ $rc -ne 0 ] || ! echo "$result" | grep -q "ok=True" || grep -qE "ERROR|runtime error|Assertion" "$log"; then
        fail "[$(basename "$fasta") $*] rc=$rc $result $(grep -m1 -E 'ERROR|runtime error|Assertion|Exception' "$log")"; return
    fi
    local length=$(echo "$result" | grep -oP '(?<=\|S\|=)\d+') runs=$(echo "$result" | grep -oP '(?<=runs=)\d+') zeros=$(echo "$result" | grep -oP '(?<=zeros=)\d+')
    local expected=$(( length + penalty * zeros ))
    [ "$objective" = runs ] && expected=$(( length + penalty * runs ))
    local reported=$(grep -o "objective value: [0-9]*" "$log" | awk '{print $3}')
    [ "$reported" = "$expected" ] && pass || fail "[$(basename "$fasta") $*] reported objective $reported, recounted $expected"
}

OPTIONS=("-k 31 -O runs" "-k 31 -O zeros" "-k 30 -O runs" "-k 30 -O runs -u" "-k 15 -O runs" "-k 31 -O runs -p 40" "-k 5 -O runs" "-k 4 -O zeros" "-k 63 -O runs")

echo "== Validity in all modes"
for mode in "" "-t $THREADS" "-m" "-m -t $THREADS"; do
    for fasta in small mid random prefix_a; do
        for options in "${OPTIONS[@]}"; do
            [[ $fasta != small && $fasta != mid && $options == *"-k 63"* ]] && continue # Sequences of length 31 have no 63-mers
            check "$WORK/$fasta.fa" $options $mode
        done
    done
done

echo "== Input formats (N's, lowercase, gzip, stdin)"
check "$WORK/mixed.fa" -k 21 -O runs
check "$WORK/mixed.fa.gz" -k 21 -O runs
"$BIN" compute -a pareto -k 21 -O runs -o "$WORK/stdin.fa" - < "$WORK/mixed.fa" > /dev/null 2>&1
$VERIFY "$WORK/mixed.fa" "$WORK/stdin.fa" 21 > /dev/null && pass || fail "reading from stdin"

echo "== Low memory mode does not change the result"
for fasta in small mid random; do
    for options in "-k 31 -O runs" "-k 31 -O zeros -u" "-k 63 -O runs" "-k 31 -O runs -p 20"; do
        "$BIN" compute -a pareto $options -o "$WORK/default.fa" "$WORK/$fasta.fa" > /dev/null 2>&1
        "$BIN" compute -a pareto $options -m -o "$WORK/low.fa" "$WORK/$fasta.fa" > /dev/null 2>&1
        cmp -s "$WORK/default.fa" "$WORK/low.fa" && pass || fail "-m changes the result [$fasta $options]"
    done
done

echo "== Large penalties (16-bit priorities, phase 3)"
for options in "-k 31 -O zeros -p 300" "-k 31 -O runs -p 248" "-k 31 -O runs -p 260" "-k 31 -O runs -p 400 -m" "-k 5 -O zeros -p 2000" "-k 31 -O zeros -p 300 -t $THREADS"; do
    check "$WORK/small.fa" $options
done
check "$WORK/mid.fa" -k 31 -O zeros -p 50

echo "== Invalid options are rejected"
for options in "-a global -O runs" "-a global -p 5" "-a global -t 2" "-a global -m" "-a pareto" "-a pareto -O length" \
               "-a pareto -O runs -p -1" "-a pareto -O runs -p 70000" "-a pareto -O runs -t 0" "-a pareto -O zeros -p 3000"; do
    "$BIN" compute -k 31 $options -o "$WORK/invalid.fa" "$WORK/small.fa" > /dev/null 2>&1 && fail "accepted: $options" || pass
done
"$BIN" compute -k 63 -a pareto -O runs -o "$WORK/invalid.fa" "$WORK/random.fa" > /dev/null 2>&1 && fail "accepted input without k-mers" || pass
for options in "-a local" "-a pareto" "-O runs"; do
    "$BIN" lowerbound -k 31 $options "$WORK/small.fa" > /dev/null 2>&1 && fail "lowerbound accepted: $options" || pass
done

echo "== Lower bounds do not exceed computed objectives"
for objective in runs zeros; do
    bound=$("$BIN" lowerbound -k 31 -a pareto -O $objective "$WORK/mid.fa" 2>/dev/null)
    "$BIN" compute -k 31 -a pareto -O $objective -o "$WORK/lb.fa" "$WORK/mid.fa" > "$WORK/lb.log" 2>&1
    value=$(grep -o "objective value: [0-9]*" "$WORK/lb.log" | awk '{print $3}')
    [ -n "$bound" ] && [ -n "$value" ] && [ "$bound" -le "$value" ] && pass || fail "lower bound $bound vs objective $value ($objective)"
done
for options in "" "-O matchtig-count"; do
    [ -n "$("$BIN" lowerbound -k 31 $options "$WORK/mid.fa" 2>/dev/null)" ] && pass || fail "lowerbound $options"
done

if [ -n "$REF_BIN" ]; then
    echo "== Same output as $REF_BIN"
    for fasta in small mid random prefix_a; do
        for options in "${OPTIONS[@]}"; do
            [[ $fasta != small && $fasta != mid && $options == *"-k 63"* ]] && continue
            "$BIN" compute -a pareto $options -o "$WORK/a.fa" "$WORK/$fasta.fa" > /dev/null 2>&1
            "$REF_BIN" compute -a pareto $options -o "$WORK/b.fa" "$WORK/$fasta.fa" > /dev/null 2>&1
            cmp -s "$WORK/a.fa" "$WORK/b.fa" && pass || fail "different output [$fasta $options]"
        done
    done
fi

echo "Passed: $PASSED, failed: $FAILED"
[ $FAILED -eq 0 ]
