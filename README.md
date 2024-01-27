# consistency_Ribo-seq
This repository demonstrates the consistency and reliability of ribosome profiling data with an analysis of 132 HeLa samples gathered from 20 distinct studies. Included are 132 aligned HeLa ribosome profiling files (ribo files), facilitating the reproduction of Spearman correlation analyses both within and across these studies. The necessary scripts for these analyses are also provided within this repository.

## Overview 
There are 2 main folders in the repo:
- script/ - houses all the codes to run the pipeline and provides information about 132 HeLa samples (HELA_list.csv)
- processed/ - houses example intermediate file and the final Spearman correlation plot, you can find the following files under this folder:
1.ribo_hela_cpm.csv CPM values for 12045 genes that have CPM > 1 across 70% of samples.
2.ribo_HeLa_Spearman.csv spearman correlation between samples.
3.dedup_ribo_HeLa_spearman_cor.pdf visualization for the result.

## Getting started
### Files you need
ribo files required to reproduce results need to be downloaded from: [zenodo link](https://zenodo.org/uploads/10565283)

### Dependencies (please install before running the scripts)
- Python + libraries (pandas, ribopy,  numpy, bioinfokit)
- R + libraries (ggplot2, tidyverse, optparse)

## Workflow
Let's go through the process: 
1. Download this repo to your computer. 
2. Download the required ribo files and create a new directory to house all the ribo files.
3. Run `python script/HeLa_correlation.py --ribonput "input folder of ribo files" --GSMinput "directory for HELA_list.csv" --outdir "you output folder"`. And you will get ribo_hela_cpm.csv in your output directory.
4. Run `Rscript script/HeLa_plot.R --CPMinput "directory of step3 output" --GSMinput "directory for HELA_list.csv" --outdir "path/to/output"`. And you will get ribo_HeLa_Spearman.csv and dedup_ribo_HeLa_spearman_cor.pdf in your output directory.
The final plot of the Spearman correlation between HeLa samples is:
https://github.com/CenikLab/consistency_Ribo-seq/blob/main/processed/dedup_ribo_HeLa_spearman_cor.pdf
