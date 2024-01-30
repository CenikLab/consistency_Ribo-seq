#!/usr/bin/env python3
"""Processes data from the HELA_list.csv"""

######use dedup only
######CDS percentage, keep 85% reads
######HeLa data
######return raw count
######dynamic cutoff for min & max length

import argparse
from bioinfokit.analys import norm
import numpy as np
from optparse import OptionParser
import os
import pandas as pd
import ribopy
from ribopy import Ribo
import scipy.stats as st

__author__ = "Yue Liu"
__credits__ = ["Yue Liu"]
__license__ = "GPLv3"
__version__ = "1.0.0"
__maintainer__ = "Yue Liu"
__email__ = "yliu5@utexas.edu"
__status__ = "Development"


####parameter check
def current_parameters(ribonput, GSMinput, outdir):
    print("input study data folder: ", ribonput)
    print("input samples list: ", GSMinput)   
    print("output results: ", outdir)


###read input experiment
def read_table(Experiment_file):
    table = pd.read_csv(Experiment_file, usecols=['study_alias', 'experiment_alias'])
    return table


####link study_id and and experiment
def study_folder(Experiment_file):
    experiment_list = read_table(Experiment_file)
    study_dic = experiment_list.groupby('study_alias')['experiment_alias'].apply(list).to_dict()
    return study_dic


####read human gene name
def defalternative_human_alias(x):
    x_pieces = x.split("|")
    return x_pieces[1] + "_" + x_pieces[5]


####read ribo files
def ribo_data(ribo_path):
    ribo_object = Ribo(ribo_path, alias=defalternative_human_alias)
    return ribo_object


####read count file by file    
def CDS_count(mmin, mmax, experiment_id, status):
    ribo_object = ribo_data("%s.ribo" % experiment_id)
    CDS_count_length_sum = ribo_object.get_region_counts(
        region_name="CDS", range_lower=int(mmin), range_upper= int(mmax), sum_lengths=True, sum_references=False,
        alias=True) ### get data sum across
    CDS_count_length_sum_fin = CDS_count_length_sum.add_suffix(status)
    return CDS_count_length_sum_fin


#####dynamic cutoff
def intevl(experiment_id):
    ribo_object = ribo_data("%s.ribo" % experiment_id)
    data_tmp=ribo_object.get_length_dist("CDS")
    data_tmp.reset_index(inplace=True)
    data=data_tmp.iloc[6:26]
    pct_85=sum(data["%s" % experiment_id]) * 0.85
    value=data[data["%s" % experiment_id] == data["%s" % experiment_id].max()]["%s" % experiment_id].values[0]
    mmin=mmax=data[data["%s" % experiment_id] == data["%s" % experiment_id].max()]['read_length'].values[0]
    while value <= pct_85 :
        if mmax < 40 and mmin > 21:
            if data[data['read_length'] == mmax+1]["%s" % experiment_id].values[0] >= \
                    data[data['read_length'] == mmin-1]["%s" % experiment_id].values[0]:
                mmax += 1
                value += data[data['read_length'] == mmax]["%s" % experiment_id].values[0]
            else:
                mmin -= 1
                value += data[data['read_length'] == mmin]["%s"%experiment_id].values[0]
        elif mmax == 40:
            mmin -= 1
            value += data[data['read_length'] == mmin]["%s" % experiment_id].values[0]
        elif mmin == 21:
            mmax += 1
            value += data[data['read_length'] == mmax]["%s" % experiment_id].values[0]
    read_pct = value / sum(data["%s" % experiment_id])
    return mmin, mmax, read_pct


#######ribo only data
def ribo_only(GSMinput,ribonput, status):
    file_dict = study_folder(GSMinput)
    new = pd.DataFrame()
    for k, v in file_dict.items():
        if status == "dedup":
            data_dir = ribonput+"/%s_dedup" % k
        try:
            os.chdir(data_dir)
        except FileNotFoundError as e:
            print(e)
            print(k)
        else:
            for j in v:
                try:
                    print("now processing file:"+j)
                    ribo_data("%s.ribo" % j)
                    mmin, mmax, select_pct = intevl(j)
                    CDS_count_length_sum = CDS_count(mmin, mmax, j, status)
                except FileNotFoundError as e:
                    print(e)
                    print(j)
                else:   
                    try:
                        new = pd.merge(new, CDS_count_length_sum, on='transcript')
                    except KeyError:
                        new = pd.DataFrame(index=CDS_count_length_sum.index)
                        new = pd.merge(new, CDS_count_length_sum, on='transcript')
    new.index = new.index.to_series().str.rsplit('_').str[1]
    return new


####prepare file for ribo only data, counts, CPM and Quantile table
####normalized
def CPM_normalize(df):
    # now, normalize raw counts using CPM method 
    nm = norm()
    nm.cpm(df=df)
    # get CPM normalized dataframe
    cpm_df = nm.cpm_norm
    #df_log=np.log2(cpm_df+1)
    return cpm_df


def data_process(df):
    #df_all = pd.read_csv(df,index_col=0)
    df_all = df
    df_count = df_all.groupby(df_all.index).mean()
    df_count.columns = df_count.columns.str.rstrip("dedup")
    df_cpm = CPM_normalize(df_count)
    return df_count, df_cpm


#for propotional data (CLR, IQLR)check rare expressed genes with CPM
def cutoff_gene_df(df, cpm_cutoff=1 ,overall_cutoff=70):
    row_cut_off = int(overall_cutoff / 100 * len(df.columns))
    df_dummy = df[(df < cpm_cutoff).sum(axis='columns') > row_cut_off]
    dummy_gene = df_dummy.index.to_series()
    return dummy_gene


def combine_cutoff_gene(dummy_gene,df):
    non_dummy = df[~df.index.isin(dummy_gene.index)]
    #dummy_df=df[df.index.isin(dummy_gene.index)]
    #dummy_result=pd.DataFrame(dummy_df.sum())
    #dummy_result.rename(columns={0:'dummy_gene'}, inplace=True)
    #frames=[non_dummy,dummy_result.T]
    #df_count_dummy=pd.concat(frames)
    return non_dummy


def main(args):
    print("#################")
    print("ribo HeLa result")
    current_parameters(args.ribonput, args.GSMinput, args.outdir)
    ribo_dedup = ribo_only(args.GSMinput, args.ribonput, status="dedup")
    #ribo_dedup.to_csv(args2.outdir+"ribo_only_raw_dedup.csv")
    print("####preprocessing raw data####")
    ribo_count_temp, ribo_CPM_temp = data_process(ribo_dedup)
    print("####dummy gene####")
    ribo_dummy_gene = cutoff_gene_df(ribo_CPM_temp, 1, 70)
    print("####generate final table####")
    #ribo_count_dummy=combine_cutoff_gene(ribo_dummy_gene,ribo_count_temp)
    ribo_CPM = ribo_CPM_temp[~ribo_CPM_temp.index.isin(ribo_dummy_gene.index)]
    ribo_CPM.to_csv(os.path.join(args.outdir, "ribo_hela_cpm.csv"))
    #r_script_path = '/path/to/your/script.R'


if __name__ == '__main__':
    """
    usage: python HeLa_correlation.py -h
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--riboinput", help='Analysis folder', default='/scratch/users/yliu5/HELA_eletter/zenodo/ribo_file/')
    parser.add_argument("--GSMinput", help='Input experiment files in csv', default='/scratch/users/yliu5/HELA_eletter/zenodo/HELA_list.csv')
    parser.add_argument("--outdir", help='Output directory', default='/scratch/users/yliu5/HELA_eletter/zenodo/processed/')
    args2 = parser.parse_args()
    main(args2)
