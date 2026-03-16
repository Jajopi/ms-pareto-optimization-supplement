# Download data and local version of kmercamel

include: "project_structure.smk"

DATASET_FILE = f"{BASE_DIR}/datasets.txt"
DATASET_URLS = dict()

def get_url_for_dataset(dataset_name: str):
    if not dataset_name in DATASET_URLS.keys():
        with open(DATASET_FILE) as file:
            for line in file.readlines():
                if len(line.strip()) == 0 or line.startswith("#"): continue
                name, url = map(lambda x: x.strip(), line.split("@"))
                name = name.replace(" ", "_")
                DATASET_URLS[name] = url
    return DATASET_URLS[dataset_name]

rule download_data:
    priority: 10
    output:
        f"{BASE_DIR}/data/{{dataset}}.fa"
    params:
        url = lambda wildcards: get_url_for_dataset(wildcards.dataset)
    shell:
        f"""
            mkdir -p $(dirname {{output}})
            if [[ ! -z "{{params.url}}" ]]; then
                curl -o {{output}}.xz {{params.url}}
                xz -d {{output}}.xz \
                    || mv {{output}}.xz {{output}} # In case the file is not compressed
            fi
        """
