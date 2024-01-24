# consistency_Ribo-seq
This repo details the consistency and reliability of ribosome profiling data through the analysis of 132 HeLa samples across 20 studies.

## Overview 
There are 2 main folders in the repo:
- script/ - houses all the code to run the pipeline and provides information about 132 HeLa samples (HELA_list.csv)
- processed/ - houses intermediate files and the final Spearman correlation plot

## Getting started
### Files you need
Original ribo files required to reproduce these results are not included in the repo for size reasons and need to be downloaded. What's required to run this out of the box is:
- Downloading the directory of ribo files: zenodo link

Here is a more detailed overview of the files I've worked on in the data directory.

| File | Description | Source | Is Included |
| --- | ----------- |--| --|
| `ribo/*` | Directory of ribo files | download from zenodo link | No (But required!) 
| `HELA_list.csv` | List of 132 HeLa samples | data obtained from GEO | Yes |

### Dependencies
- Python + libraries (pandas, ribopy,  numpy, bioinfokit)
- R + libraries (ggplot2, tidyverse)

## Workflow
Let's go through the process: 
1. Create a new directory to house all the ribo files.
2. Double the sample information in the file `script/HELA_list.csv`.
3. Run `python script/HeLa_extract_count.py --ribinput "input folder of ribo files" --GSMinput "directory for HELA_list.csv" --outdir "you output folder"`; you may also change the directories directly in the script.
4. Run `python ribobase_counts_processing.py -i "ribo output form step3" -m "only"`, you can put the script and the result from step3 in the same folder.
5. Run `correlation_ribo.R` with your CPM result from step4.
The final plot of the Spearman correlation between HeLa samples is:
[embed]https://github.com/CenikLab/consistency_Ribo-seq/blob/main/processed/dedup_ribo_HeLa_spearman_cor.pdf[/embed]
