# analyse the GM result and output the features for MRR.
import sys
sys.path.append("..")
from ana_res import pp_utils
import pandas as pd
from matplotlib import pyplot as plt


def stat(gtf_f):
    df = pp_utils.load_gff2( gtf_f)
    grouped_df = df.groupby(['seqid'])
    seq_out_l = [grouped_df.get_group(x) for x in grouped_df.groups]
    rows = []
    for i, seq_df in enumerate(seq_out_l):
        """
        rows.append([seq_df.iloc[0]['seqid'].split(':')[0], \
                    seq_df.iloc[0]['seqid'].split(':')[1], \
                    len(seq_df), seq_df['len'].sum(), seq_df['len'].max(), \
                    seq_df['score'].sum(), seq_df['score'].max()])
        """
        rows.append([seq_df.iloc[0]['seqid'],
                    len(seq_df), seq_df['len'].sum(), seq_df['len'].max(), \
                    seq_df['score'].sum(), seq_df['score'].max()])
        #print(seq_df.iloc[0]['seqid'], len(seq_df), seq_df['len'].sum(), seq_df['len'].max(), seq_df['score'].sum(), seq_df['score'].max())

    #stat_df = pd.DataFrame(rows, columns = ['seqid', 'path_score','n_frag', 'sum_len', 'max_len','sum_score', 'max_score'])
    stat_df = pd.DataFrame(rows, columns = ['seqid', 'n_frag', 'sum_len', 'max_len','sum_score', 'max_score'])
    return df, stat_df

def gtf2csv(gtf_f):
    df, stat_df = stat(gtf_f)
    df.to_csv(gtf_f + ".csv")
    stat_df.to_csv(gtf_f + "stat.csv")

#"/home/runzhouyu2/work/aligngraph-master/data/2_geno_mix/mix_6/grouped_beam_50.gtf"
gtf_f = sys.argv[1]
gtf2csv(gtf_f)