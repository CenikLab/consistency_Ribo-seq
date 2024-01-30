# consistency_Ribo-seq
This repository demonstrates the consistency and reliability of ribosome profiling data with an analysis of 132 HeLa samples gathered from 20 distinct studies. Included are 132  ribosome profiling files (ribo files), facilitating the reproduction of Spearman correlation analyses both within and across these studies. Detailed information about the ribo file format can be found [here](https://github.com/ribosomeprofiling/riboflow). The necessary scripts for these analyses are also provided within this repository.

To generate the ribo files, duplicated reads based on the same length and start position were removed. Then we implemented a read selection module to extract ribosome protected footprints by automatically determining the lower and upper fragment lengths based on the positions with the highest number of CDS-aligned reads. This process, involving the extension and adjustment of boundaries to include at least 85% of total CDS reads, ensured the selection of the most informative ribosome protected footprints. Using this approach we quantified ribosome occupancy for more than 12,000 genes (counts per million reads exceeding one in more than 70% of the samples). We used the CPM values to calculate the Spearman correlation between samples.

## Overview 
There are two directories in this repository:

- script - Code to run the pipeline and information about 132 HeLa samples (HELA\_list.csv)

- processed - Example intermediate file and the final Spearman correlation plot. You can find the following files under this folder:

1. ribo\_hela\_cpm.csv CPM values for 12045 genes that have CPM > 1 across 70% of samples.

2. ribo\_HeLa\_Spearman.csv spearman correlation between samples.

3. dedup\_ribo\_HeLa\_spearman\_cor.pdf visualization for the result.

## Getting started

### Input files
ribo files required to reproduce results need to be downloaded from this [Zenodo link](https://zenodo.org/uploads/10565283).

## Workflow
Let's walk through the process:

1. Download this repository and create a temporary output directory.
```
git clone https://github.com/CenikLab/consistency_Ribo-seq && cd consistency_Ribo-seq
WORK_DIR="/tmp/consistency_Ribo-seq_test"
```

2. Install [conda](https://docs.conda.io/en/latest/miniconda.html) for required software dependencies.

3. Run the following to generate the conda environment and activate it:
```
conda env create -f environment.yaml && conda activate ribo_consistency
```

3. Download the required ribo files and create a new directory to house all the ribo files.
```
cd "$WORK_DIR"
curl https://zenodo.org/uploads/10565283 && cd -
```

4. Run the following:
```
python $PWD/script/HeLa_correlation.py --riboinput "$WORK_DIR" --GSMinput $PWD/script/HELA_list.csv --outdir "$WORK_DIR"
```
You will get ribo\_hela\_cpm.csv in your output directory.

5. Finally run the following to generate a plot:
```
Rscript $PWD/script/HeLa_plot.R --CPMinput "$WORK_DIR" --GSMinput $PWD/script/HELA_list.csv --outdir "$WORK_DIR"
```

You will get ribo\_HeLa\_Spearman.csv and dedup\_ribo\_HeLa\_spearman\_cor.pdf in your output directory.

The final plot of the Spearman correlation between HeLa samples should be:
![consistency_Ribo-seq](https://github.com/CenikLab/consistency_Ribo-seq/blob/main/processed/dedup_ribo_HeLa_spearman_cor.jpg "compare")
