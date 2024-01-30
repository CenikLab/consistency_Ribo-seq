# consistency_Ribo-seq
This repository demonstrates the consistency and reliability of ribosome profiling data with an analysis of 132 HeLa samples gathered from 20 distinct studies. Included are 132 aligned HeLa ribosome profiling files (ribo files, detailed information about ribo file can be found: https://github.com/ribosomeprofiling/riboflow), facilitating the reproduction of Spearman correlation analyses both within and across these studies. The necessary scripts for these analyses are also provided within this repository. The ribo files removed duplicated reads based on the same length and position. Then we implemented a read selection module to extract ribosome-protected fragments in this script, automatically determining the lower and upper fragment lengths based on the positions with the highest number of CDS-aligned reads. This process, involving the extension and adjustment of boundaries to include at least 85% of total CDS reads, ensured the selection of the most informative ribosome-protected mRNA footprint positions. Using this approach we quantified ribosome occupancy for more than 12,000 genes (counts per million reads exceeding one in more than 70% of the samples). We used the CPM values to calculate the Spearman correlation between samples.

## Overview 
There are 2 main folders in the repo:
- script/ - houses all the codes to run the pipeline and provides information about 132 HeLa samples (HELA_list.csv)
- processed/ - houses example intermediate file and the final Spearman correlation plot, you can find the following files under this folder:
1.ribo_hela_cpm.csv CPM values for 12045 genes that have CPM > 1 across 70% of samples.
2.ribo_HeLa_Spearman.csv spearman correlation between samples.
3.dedup_ribo_HeLa_spearman_cor.pdf visualization for the result.

## Getting started
### Files you need
ribo files (removed duplicated readds based on same position and length) required to reproduce results need to be downloaded from: [zenodo link](https://zenodo.org/uploads/10565283)

### Dependencies (please install before running the scripts)
- Python + libraries (pandas, ribopy,  numpy, bioinfokit)
- R + libraries (ggplot2, tidyverse, optparse)

## Workflow
Let's go through the process: 
1. Download this repo to your computer. 
2. Download the required ribo files and create a new directory to house all the ribo files.
3. Run `python consistency_Ribo-seq-main/script/HeLa_correlation.py --ribonput "input folder of ribo files" --GSMinput "directory for HELA_list.csv" --outdir "you output folder"`. And you will get ribo_hela_cpm.csv in your output directory.
4. Run `Rscript consistency_Ribo-seq-main/script/HeLa_plot.R --CPMinput "directory of step3 output" --GSMinput "directory for HELA_list.csv" --outdir "path/to/output"`. And you will get ribo_HeLa_Spearman.csv and dedup_ribo_HeLa_spearman_cor.pdf in your output directory.
The final plot of the Spearman correlation between HeLa samples is:
![consistency_Ribo-seq](https://github.com/CenikLab/consistency_Ribo-seq/blob/main/processed/dedup_ribo_HeLa_spearman_cor.jpg "compare")
