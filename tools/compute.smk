# Compute various (masked) superstring representations, benchmark some of them

include: "project_structure.smk"

GGCAT_FLAGS = "-j 1 -s 1 -m 2 -p -t .temp"
JELLYFISH_FLAGS = "-s 500M -t 1 -C"

LOCAL_KMERCAMEL = f"{BASE_DIR}/kmercamel/kmercamel-pareto/kmercamel"
KMER_COUNTER_MS_SCRIPT = f"{BASE_DIR}/tools/kmer_counter_ms.py"
SUPERSTRING_SCRAMBLE_SCRIPT = f"{BASE_DIR}/tools/scramble_superstring.py"
BLOSSOM5 = f"{BASE_DIR}/tools/blossom5/blossom5"

rule compute_kmer_count_jellyfish:
    input:
        ancient(f"{DATA_DIR}/{{dataset}}.fa")
    output:
        f"{SUPERSTRINGS}/{{dataset}}/{{k}}/kmer_count.txt"
    shell:
        f"""
            mkdir -p $(dirname {{output}})
            jellyfish count {JELLYFISH_FLAGS} -m {{wildcards.k}} -o {{output}}.tmp {{input}}
            jellyfish stats {{output}}.tmp | head -n 2 | tail -n 1 | tr -s ' ' | cut -d ' ' -f 2 > {{output}}
            rm {{output}}.tmp
        """

rule compute_lowerbound_matchtig_count:
    input:
        ancient(f"{DATA_DIR}/{{dataset}}.fa")
    output:
        f"{SUPERSTRINGS}/{{dataset}}/{{k}}/lowerbound_matchtig_count.txt"
    shell:
        f"""
            mkdir -p $(dirname {{output}})
            {LOCAL_KMERCAMEL} lowerbound -k {{wildcards.k}} -O matchtig-count {{input}} > {{output}}
        """

rule compute_lowerbound_length_kmercamel:
    input:
        ancient(f"{DATA_DIR}/{{dataset}}.fa")
    output:
        f"{SUPERSTRINGS}/{{dataset}}/{{k}}/lowerbound_length.txt"
    shell:
        f"""
            mkdir -p $(dirname {{output}})
            kmercamel lowerbound -k {{wildcards.k}} {{input}} > {{output}}
        """

rule compute_greedy_simplitigs_kmercamel:
    input:
        ancient(f"{DATA_DIR}/{{dataset}}.fa")
    output:
        f"{SUPERSTRINGS}/{{dataset}}/{{k}}/simplitigs.fa"
    benchmark:
        f"{BENCHMARK_DIR}/kmercamel_simplitigs/{{dataset}}_{{k}}.tsv"
    shell:
        f"""
            mkdir -p $(dirname {{output}})
            kmercamel compute -a local -d 1 -k {{wildcards.k}} {{input}} > {{output}}
        """

rule compute_optimal_matchtigs_blossom5:
    input:
        unitigs = ancient(f"{SUPERSTRINGS}/{{dataset}}/{{k}}/unitigs.no-ms.fa"),
        blossom5 = ancient(BLOSSOM5)
    output:
        f"{SUPERSTRINGS}/{{dataset}}/{{k}}/optimal_matchtigs.no-ms.fa"
    benchmark:
        f"{BENCHMARK_DIR}/optimal_matchtigs/{{dataset}}_{{k}}.tsv"
    shell:
        f"""
            mkdir -p $(dirname {{output}})
            matchtigs -k {{wildcards.k}} -t 1 --fa-in {{input.unitigs}} --matchtigs-fa-out {{output}} --blossom5-command {{input.blossom5}} 2> /dev/null
        """

def get_ggcat_flag_for_input(input_type):
    if input_type == "matchtigs": return "--greedy-matchtigs"
    if input_type == "eulertigs": return "--eulertigs"
    if input_type == "unitigs": return ""
    raise Exception(f"Invalid input type specified: {input_type}.")

rule compute_spss_ggcat:
    input:
        ancient(f"{DATA_DIR}/{{dataset}}.fa")
    output:
        f"{SUPERSTRINGS}/{{dataset}}/{{k}}/{{representation}}.no-ms.fa" # matchtigs or eulertigs
    wildcard_constraints:
        representation="(unitigs|eulertigs|matchtigs)"
    params:
        type_flag = lambda wildcards: get_ggcat_flag_for_input(wildcards.representation)
    benchmark:
        f"{BENCHMARK_DIR}/ggcat_{{representation}}/{{dataset}}_{{k}}.tsv"
    shell:
        f"""
            mkdir -p $(dirname {{output}})
            ggcat build {{params.type_flag}} -k {{wildcards.k}} {{input}} {GGCAT_FLAGS} -o {{output}}
        """

rule convert_ggcat_or_matchtigs_output_to_ms:
    input:
        ancient(f"{SUPERSTRINGS}/{{dataset}}/{{k}}/{{representation}}.no-ms.fa")
    output:
        f"{SUPERSTRINGS}/{{dataset}}/{{k}}/{{representation}}.fa"
    wildcard_constraints:
        representation="(unitigs|eulertigs|matchtigs|optimal_matchtigs)"
    shell:
        f"""
            kmercamel spss2ms -k {{wildcards.k}} {{input}} > {{output}}
        """

rule compute_ms_pareto:
    input:
        ancient(f"{DATA_DIR}/{{dataset}}.fa")
    output:
        f"{SUPERSTRINGS}/{{dataset}}/{{k}}/pareto-{{run_penalty}}.fa"
    wildcard_constraints:
        run_penalty="\\d+"
    benchmark:
        f"{BENCHMARK_DIR}/pareto_optimization/{{dataset}}_{{k}}_{{run_penalty}}.tsv"
    shell:
        f"""
            mkdir -p $(dirname {{output}})
            {LOCAL_KMERCAMEL} compute -a pareto -O runs -p {{wildcards.run_penalty}} -k {{wildcards.k}} {{input}} > {{output}}
        """

rule compute_ms_minone_kmercamel:
    input:
        ancient(f"{DATA_DIR}/{{dataset}}.fa")
    output:
        f"{SUPERSTRINGS}/{{dataset}}/{{k}}/greedy-minone.fa"
    benchmark:
        f"{BENCHMARK_DIR}/kmercamel_greedy_minone/{{dataset}}_{{k}}.tsv"
    shell:
        f"""
            mkdir -p $(dirname {{output}})
            kmercamel compute -k {{wildcards.k}} {{input}} > {{output}}
        """

rule optimize_mask_kmercamel:
    input:
        ancient(f"{SUPERSTRINGS}/{{dataset}}/{{k}}/greedy-minone.fa")
    output:
        f"{SUPERSTRINGS}/{{dataset}}/{{k}}/greedy-{{mask_opt_criterion}}.fa"
    wildcard_constraints:
        mask_opt_criterion="(?!minone).*" # Supress cyclic dependency error
    benchmark:
        f"{BENCHMARK_DIR}/kmercamel_mask_opt_{{mask_opt_criterion}}/{{dataset}}_{{k}}.tsv"
    shell:
        f"""
            kmercamel maskopt -k {{wildcards.k}} -t {{wildcards.mask_opt_criterion}} {{input}} > {{output}}
        """

rule split_ms_into_mask_and_superstring:
    input:
        ancient(f"{SUPERSTRINGS}/{{dataset}}/{{k}}/{{method}}.fa")
    output:
        superstring = f"{SUPERSTRINGS}/{{dataset}}/{{k}}/{{method}}.superstring",
        mask = f"{SUPERSTRINGS}/{{dataset}}/{{k}}/{{method}}.mask"
    shell:
        f"""
            kmercamel ms2mssep -m {{output.mask}} -s {{output.superstring}} {{input}}
        """
