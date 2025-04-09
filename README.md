# flu-syn-rates

## Notebooks
* `compute_counts.ipynb`: reads in tree and computes counts
* `analyze_counts.ipybn`: analyzes counts from above notebook
* `analyze_metadata.ipynb`: analyzes the metadata file

## Snakemake pipeline

### Prerequisites
Install the two conda/mamba environments from file:
```
mamba env create -f environment.yml
mamba env create -f environment_usher.yml
```

Run the pipeline:
```
snakemake -j 1 --use-conda
```
The current version does not make use of multiple cores, but it's not hard.
