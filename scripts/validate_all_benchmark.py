from openmc_fusion_benchmarks.validate_spec import validate_benchmark
from pathlib import Path

benchmarks_dir = Path("src/openmc_fusion_benchmarks/benchmarks")

# Loop through each benchmark directory
for benchmark in benchmarks_dir.iterdir():
    if benchmark.is_file():
        continue  # Skip files, process only directories

    validate_benchmark(benchmark.name)  # Pass only the directory name
