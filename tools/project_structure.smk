# Define global path variables

from os import path, getcwd
BASE_DIR = path.abspath(getcwd()) + "/.."

DATA_DIR = f"{BASE_DIR}/data"
SUPERSTRINGS = f"{BASE_DIR}/computed"
BENCHMARK_DIR = f"{BASE_DIR}/benchmarks"

RESULTS = "results"
