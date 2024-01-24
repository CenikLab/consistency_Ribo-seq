######author Yue Liu, email yliu5@utexas.edu 
######process data from raw reads counts talbe
######return processed counts, CPM
######run python ribobase_counts_processing.py -h for help
######ribo only:ribobase_counts_processing.py -i "ribo_input" -m "only"

from bioinfokit.analys import norm
from optparse import OptionParser
import pandas as pd

####input paramater####
parser = OptionParser()
parser.add_option("-i", "--input_ribo_file", help='ribo raw file', dest='ribof')
parser.add_option("-r", "--input_rna_file", help='rna raw file', dest='rnaf')
parser.add_option("-c", "--cpm_cut_off", help='filter genes step1', dest='cpm_cutoff',default=int(1))
parser.add_option("-a", "--overall_cut_off", help='filter genes step2', dest='overall_cutoff',default=int(70))
parser.add_option("-m", "--mode", help='ribo only or paired', dest='mode')
parser.add_option('-o', '--dir', help='work directory', dest='workdir', default='.')
(options, args) = parser.parse_args()

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
    df_all = pd.read_csv(df,index_col=0)
    df_count=df_all.groupby(df_all.index).mean()
    df_count.columns=df_count.columns.str.rstrip("dedup")
    df_cpm=CPM_normalize(df_count)
    return df_count,df_cpm

#for propotional data (CLR, IQLR)check rare expressed genes with CPM
def dummy_gene_df(df,cpm_cutoff=1,overall_cutoff=70):
    row_cut_off = int(overall_cutoff/100*len(df.columns))
    df_dummy = df[(df<cpm_cutoff).sum(axis='columns') > row_cut_off]
    dummy_gene=df_dummy.index.to_series()
    return dummy_gene

def combine_dummy_gene(dummy_gene,df):
    non_dummy=df[~df.index.isin(dummy_gene.index)]
    dummy_df=df[df.index.isin(dummy_gene.index)]
    dummy_result=pd.DataFrame(dummy_df.sum())
    dummy_result.rename(columns={0:'dummy_gene'}, inplace=True)
    frames=[non_dummy,dummy_result.T]
    df_count_dummy=pd.concat(frames)
    return df_count_dummy

####prepare file for paired data, counts, CPM and Quantile table
if options.mode == "only":
    print("####preprocessing raw data####")
    ribo_count_temp,ribo_CPM_temp = data_process(options.ribof)
    print("####dymmy gene####")
    ribo_dummy_gene=dummy_gene_df(ribo_CPM_temp,options.cpm_cutoff,options.overall_cutoff)
    print("####generate final table####")
    ribo_count_dummy=combine_dummy_gene(ribo_dummy_gene,ribo_count_temp)
    ribo_CPM=ribo_CPM_temp[~ribo_CPM_temp.index.isin(ribo_dummy_gene.index)]
    #ribo_Q=ribo_Q_temp[~ribo_Q_temp.index.isin(ribo_dummy_gene.index)]
    ribo_count_dummy.to_csv(options.workdir + "/ribo_only_count_dummy_%s.csv"%(options.overall_cutoff))
    ribo_CPM.to_csv(options.workdir + "/ribo_only_cpm_dummy_%s.csv"%(options.overall_cutoff))
    #ribo_Q.to_csv(options.workdir + "/ribo_only_quantile_dummy_%s.csv"%(options.overall_cutoff))
