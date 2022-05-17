#################################################################################$$
# Author: Runzhou
# Function: 
# store some commonly used functions to analyse experiment results
# Supported file formats:
# 1. genemark output -gtf 
# 2. blast web output - csv
#################################################################################$$
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
import seaborn as sns

"""
data example:
HIV_1	GeneMark.hmm	CDS	336	1838	56.487819	+	0	gene_id=1

call by:
G1_k20_df = load_gff2("/home/runzhouyu2/work/aligngraph-master/src/BeamSearch/gm_res/G1_k20.gtf")
"""
gff_col = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']

def cal_len(row):
    return abs(row['end'] - row['start'])
def load_gff2(score_file):
    res_list = []
    with open(score_file, 'r') as score_f:
        line = score_f.readline()
        while line:
            if len(line.split('\t')) >= 9 :
                res_list.append(line.split('\t')[0:9])
            line = score_f.readline()
    score_df = pd.DataFrame(res_list, index=None, columns=gff_col)
    score_df['score'] = score_df['score'].astype(float)
    score_df['start'] = score_df['start'].astype(int)
    score_df['end'] = score_df['end'].astype(int)
    score_df['len'] = score_df.apply (lambda row: cal_len(row), axis=1)
    #print("score info:\n",score_df['score'].describe())
    #print(score_df['score'].value_counts())
    return score_df

"""
data example:
1_1	YU2_genome	95.761	8800	287	59	639	9397	8840	86	0	14105

call by:
raw_consen_2HIV = load_bls("blast_res/raw_consen_2HIV_3G.csv")
"""
bls_col = ['que_seq', 'sub_seq', 'identity', 'align_len', 'mismatches', 'gap_opens', 'q_start', 'q_end', 's_start', 's_end', 'evalue', 'bit_score']
#bls_col = ['query', 'subject', 'identity', 'align_len', 'mismatches', 'gap_opens', 'q_start', 'q_end', 's_start', 's_end', 'evalue', 'bit_score']
def load_bls(bls_f):
    return pd.read_csv(bls_f, names = bls_col, header = None)
    #return bls_df.loc[bls_f['align_len'] > 5000]