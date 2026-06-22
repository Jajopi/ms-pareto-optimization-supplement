# Test k-mer sets in masked superstring representations and compressed representations

include: "project_structure.smk"

JELLYFISH_FLAGS = "-s 500M -t 1 -C"

rule compute_kmer_count_in_ms_jellyfish:
    input:
        f"{SUPERSTRINGS}/{{dataset}}/{{k}}/{{method}}.fa.{{decompression_type}}{{merged_or_not}}"
    output:
        f"{SUPERSTRINGS}/{{dataset}}/{{k}}/{{method}}.fa.{{decompression_type}}{{merged_or_not}}.kmer_count.txt"
    wildcard_constraints:
        decompression_type="(de-co|de-co.max|no-comp)",
        merged_or_not=".*"
    shell:
        f"""
            mkdir -p $(dirname {{output}})
            jellyfish count {JELLYFISH_FLAGS} -m {{wildcards.k}} -o {{output}}.tmp {{input}}
            jellyfish stats {{output}}.tmp | head -n 2 | tail -n 1 | tr -s ' ' | cut -d ' ' -f 2 > {{output}}
            rm {{output}}.tmp
        """

rule merge_kmer_sets_for_count_testing:
    input:
        f"{DATA_DIR}/{{dataset}}.fa",
        f"{SUPERSTRINGS}/{{dataset}}/{{k}}/{{method}}.fa.{{decompression_type}}"
    output:
        f"{SUPERSTRINGS}/{{dataset}}/{{k}}/{{method}}.fa.{{decompression_type}}.merged"
    wildcard_constraints:
        decompression_type="(de-co|de-co.max|no-comp)"
    shell:
        f"""
            mkdir -p $(dirname {{output}})
            cat {{input[0]}} {{input[1]}} > {{output}}
        """

rule prepare_not_compressed_ms_for_jellyfish:
    input:
        f"{SUPERSTRINGS}/{{dataset}}/{{k}}/{{method}}.fa"
    output:
        f"{SUPERSTRINGS}/{{dataset}}/{{k}}/{{method}}.fa.no-comp"
    shell:
        f"""
            kmercamel ms2spss -k {{wildcards.k}} -o {{output}} {{input}}
        """

rule check_kmer_set_equality_in_regular_and_compressed_files:
    input:
        right_count_file=ancient(f"{SUPERSTRINGS}/{{dataset}}/{{k}}/kmer_count.txt"),
        tested_count_file=ancient(f"{SUPERSTRINGS}/{{dataset}}/{{k}}/{{method}}.fa.{{decompression_type}}.kmer_count.txt"),
        merged_count_file=ancient(f"{SUPERSTRINGS}/{{dataset}}/{{k}}/{{method}}.fa.{{decompression_type}}.merged.kmer_count.txt")
    output:
        f"{SUPERSTRINGS}/{{dataset}}/{{k}}/{{method}}.fa.{{decompression_type}}.check"
    wildcard_constraints:
        decompression_type="(de-co|de-co.max|no-comp)"
    shell:
        f"""
            mkdir -p $(dirname {{output}})
            RIGHT_COUNT=$(cat {{input.right_count_file}})
            ACTUAL_COUNT=$(cat {{input.tested_count_file}})
            MERGED_COUNT=$(cat {{input.merged_count_file}})
            if [[ "$RIGHT_COUNT" -ne "$ACTUAL_COUNT" ]] || [[ "$RIGHT_COUNT" -ne "$MERGED_COUNT" ]]; then
                echo "k-mer count mismatch for {SUPERSTRINGS}/{{wildcards.dataset}}/{{wildcards.k}}/{{wildcards.method}}.fa{{wildcards.decompression_type}}: right count = $RIGHT_COUNT, actual count = $ACTUAL_COUNT, merged count = $MERGED_COUNT" >&2
                exit 42
            fi
            touch {{output}}
        """
