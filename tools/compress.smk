# Compress superstring and mask files

include: "project_structure.smk"

MASK_TO_POSITION_ENCODER = f"{BASE_DIR}/tools/mask_to_run_encoding.py"
ELIAS_FANO_COMPRESSOR = f"{BASE_DIR}/tools/EliasFano/EliasFano"

BZIP2_FLAGS = "--best"
GECO_FLAGS = "-l5"
GECO_MAX_FLAGS = "-lr 0.005 -hs 160 \
    -tm 1:1:1:0:0.6/0:0:0 -tm 1:1:0:0:0.6/0:0:0 -tm 2:1:2:0:0.90/0:0:0 \
    -tm	2:1:1:0:0.8/0:0:0 -tm 3:1:0:0:0.8/0:0:0 -tm 4:1:0:0:0.8/0:0:0 -tm 5:1:0:0:0.8/0:0:0 -tm 6:1:0:0:0.8/0:0:0 \
	-tm 7:1:1:0:0.7/0:0:0 -tm 8:1:0:0:0.85/0:0:0 -tm 9:1:1:0:0.88/0:0:0 -tm 11:10:2:0:0.9/0:0:0 \
	-tm	11:10:0:0:0.88/0:0:0 -tm 12:20:1:0:0.88/0:0:0 -tm 14:50:1:1:0.89/1:10:0.89 -tm 17:2000:1:10:0.88/2:50:0.88 \
	-tm 20:1200:1:160:0.88/3:15:0.88"
	# Source: https://github.com/cobilab/HumanGenome/blob/786d26090e5d804641dc424f433978c7475468f5/scripts/Run39.sh#L16
XZ_FLAGS = "-9"
XZ_MAX_FLAGS = "-e --lzma2=preset=9,dict=1500MiB,nice=250"

rule encode_spss_mask_into_run_lengths_and_k:
    input:
        ancient(f"{SUPERSTRINGS}/{{dataset}}/{{k}}/{{method}}.mask")
    output:
        f"{SUPERSTRINGS}/{{dataset}}/{{k}}/{{method}}.mask.lengths"
    wildcard_constraints:
        method="(unitigs|simplitigs|eulertigs|matchtigs|optimal_matchtigs)"
    shell:
        f"""
            {MASK_TO_POSITION_ENCODER} < {{input}} > {{output}}
        """

rule compress_elias_fano:
    input:
        ancient(f"{SUPERSTRINGS}/{{dataset}}/{{k}}/{{method}}.mask")
    output:
        f"{SUPERSTRINGS}/{{dataset}}/{{k}}/{{method}}.mask.efe"
    shell:
        f"""
            mkdir -p $(dirname {{output}})
            {ELIAS_FANO_COMPRESSOR} < {{input}} > {{output}}
        """

rule compress_bzip2:
    input:
        ancient(f"{SUPERSTRINGS}/{{dataset}}/{{k}}/{{method}}.{{file_type}}")
    output:
        f"{SUPERSTRINGS}/{{dataset}}/{{k}}/{{method}}.{{file_type}}.bz2"
    shell:
        f"""
            bzip2 -z -k {BZIP2_FLAGS} {input}
        """

rule compress_geco3:
    input:
        ancient(f"{SUPERSTRINGS}/{{dataset}}/{{k}}/{{method}}.{{file_type}}")
    output:
        f"{SUPERSTRINGS}/{{dataset}}/{{k}}/{{method}}.{{file_type}}.co"
    shell:
        f"""
            cat {{input}} | tr '01' 'AC' > {{input}}.geco.tmp
            GeCo3 {GECO_FLAGS} {{input}}.geco.tmp
            mv {{input}}.geco.tmp.co {{output}}
            rm {{input}}.geco.tmp
        """

rule compress_geco3_max:
    input:
        ancient(f"{SUPERSTRINGS}/{{dataset}}/{{k}}/{{method}}.{{file_type}}")
    output:
        f"{SUPERSTRINGS}/{{dataset}}/{{k}}/{{method}}.{{file_type}}.co.max"
    resources:
        mem_mb=20000
    shell:
        f"""
            cat {{input}} | tr '01' 'AC' > {{input}}.geco-max.tmp
            GeCo3 {GECO_MAX_FLAGS} {{input}}.geco-max.tmp
            mv {{input}}.geco-max.tmp.co {{output}}
            rm {{input}}.geco-max.tmp
        """

rule compress_xz:
    input:
        ancient(f"{SUPERSTRINGS}/{{dataset}}/{{k}}/{{method}}.{{file_type}}")
    output:
        f"{SUPERSTRINGS}/{{dataset}}/{{k}}/{{method}}.{{file_type}}.xz",
    shell:
        f"""
            xz -z -T1 {XZ_FLAGS} -c {{input}} > {{output}}
        """

rule compress_xz_max:
    input:
        ancient(f"{SUPERSTRINGS}/{{dataset}}/{{k}}/{{method}}.{{file_type}}")
    output:
        f"{SUPERSTRINGS}/{{dataset}}/{{k}}/{{method}}.{{file_type}}.xz.max",
    resources:
        mem_mb=15000
    shell:
        f"""
            xz -z -T1 {XZ_MAX_FLAGS} -c {{input}} > {{output}}
        """
