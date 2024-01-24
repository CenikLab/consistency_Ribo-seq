######author Yue Liu, email yliu5@utexas.edu 
######process data from the HELA_list.csv
######use dedup only
######CDS percentage, keep 85% reads
######HeLa data
######return raw count
######dynamic cutoff for min & max length
#!/usr/bin/env python3
import ribopy
from ribopy import Ribo
import os
import pandas as pd
import argparse
import shlex
import numpy as np
import math
import scipy.stats as st

####parameter check
def current_parameters(ribinput, GSMinput,outdir):
    print("input study data folder: ",ribinput)
    print("input samples list: ", GSMinput)   
    print("output results: ", outdir)

###read input experiment
def read_table(Experiment_file):
    table=pd.read_csv(Experiment_file,usecols=['study_alias','experiment_alias'])
    return table

####link study_id and and experiment
def study_folder(Experiment_file):
    study_dic={}
    experiment_list=read_table(Experiment_file)
    study_dic=experiment_list.groupby('study_alias')['experiment_alias'].apply(list).to_dict()
    return study_dic

####read human gene name
def defalternative_human_alias(x):
    x_pieces = x.split("|")
    return x_pieces[1] + "_" + x_pieces[5]

####read ribo files
def ribo_data(ribo_path):
    ribo_object = Ribo(ribo_path,alias = defalternative_human_alias)
    return ribo_object

####read count file by file    
def CDS_count(mmin,mmax,experiment_id,status):
    ribo_object = ribo_data("%s.ribo"%experiment_id)
    CDS_count_length_sum = ribo_object.get_region_counts(region_name = "CDS",
                              range_lower    = int(mmin),
                              range_upper    = int(mmax),
                              sum_lengths    = True,
                              sum_references = False,
                              alias          = True) ### get data sum across
    CDS_count_length_sum_fin=CDS_count_length_sum.add_suffix(status)
    return CDS_count_length_sum_fin

#####dynamic cutoff
def intevl(experiment_id):
    ribo_object = ribo_data("%s.ribo"%experiment_id)
    data_tmp=ribo_object.get_length_dist("CDS")
    data_tmp.reset_index(inplace=True)
    data=data_tmp.iloc[6:26]
    pct_85=sum(data["%s"%experiment_id])*0.85
    value=data[data["%s"%experiment_id]==data["%s"%experiment_id].max()]["%s"%experiment_id].values[0]
    mmin=mmax=data[data["%s"%experiment_id]==data["%s"%experiment_id].max()]['read_length'].values[0]
    while value<=pct_85 :
        if mmax<40 and mmin>21: 
            if data[data['read_length']==mmax+1]["%s"%experiment_id].values[0] >= data[data['read_length']==mmin-1]["%s"%experiment_id].values[0]:
                mmax+=1
                value+=data[data['read_length']==mmax]["%s"%experiment_id].values[0]
            else:
                mmin-=1
                value+=data[data['read_length']==mmin]["%s"%experiment_id].values[0]
        elif mmax==40:
            mmin-=1
            value+=data[data['read_length']==mmin]["%s"%experiment_id].values[0]
        elif mmin==21:
            mmax+=1
            value+=data[data['read_length']==mmax]["%s"%experiment_id].values[0]
    read_pct=value/sum(data["%s"%experiment_id])
    return mmin,mmax,read_pct

#######ribo only data
def ribo_only(GSMinput,ribinput,outdir,status):
    file_dict=study_folder(GSMinput)
    new=pd.DataFrame()
    for k,v in file_dict.items():
        if status == "dedup":
            data_dir=ribinput+"%s_dedup"%(k)
        try:
            os.chdir(data_dir)
        except FileNotFoundError as e:
            print(e)
            print(k)
        else:
            for j in v:
                try:
                    print("now processing file:"+j)
                    ribo_data("%s.ribo"%j)
                    mmin,mmax,select_pct=intevl(j)
                    CDS_count_length_sum=CDS_count(mmin,mmax,j,status)
                except FileNotFoundError as e:
                    print(e)
                    print(j)
                else:   
                    try:
                        new=pd.merge(new,CDS_count_length_sum,on='transcript')                      
                    except KeyError:
                            new=pd.DataFrame(index=CDS_count_length_sum.index)
                            new=pd.merge(new,CDS_count_length_sum,on='transcript')
    new.index = new.index.to_series().str.rsplit('_').str[1]
    return new


def main(args2):
    print("#################")
    print("ribo HeLa result")
    current_parameters(args2.ribinput, args2.GSMinput,args2.outdir)
    ribo_dedup=ribo_only(args2.GSMinput,args2.ribinput,args2.outdir,status="dedup")
    ribo_dedup.to_csv(args2.outdir+"ribo_only_raw_dedup.csv")

if __name__ == '__main__':
    """
    usage: python HeLa_extract_count.py
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--ribinput", help='all analysis folder',default='/scratch/users/yliu5/HELA_eletter/zenodo/ribo_file/')
    parser.add_argument("--GSMinput", help='input experiment files in csv', default='/scratch/users/yliu5/HELA_eletter/zenodo/HELA_list.csv')
    parser.add_argument("--outdir", help='output directory', default='/scratch/users/yliu5/HELA_eletter/zenodo/processed/')
    args2 = parser.parse_args()
    main(args2)