#################################################################################$$
# Author: Runzhou
# Function: 
# After running GM, this code analyse the GTF outputs from GM and apply MRR to select best seq.
# dafault features for MRR:  
# 1. seq length (descending)
# 2. number of genes (ascending)
# 3. max gene score (descending)
# 4. total gene score (descending)
#################################################################################$$
import argparse
import ana_GM
import pandas as pd
from Bio import SeqIO

if __name__ == '__main__':
    #parse args
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', type=str, required=True, help = "Reads output for further ranking(in fasta format).")
    #parser.add_argument('--GM', type=str, required=False, help = "Beamwidth for diverse beam search (default: 500).")
    args = parser.parse_args()
    reads = args.r
    GM_file = reads + ".gtf"
    final_out_file = reads + "_final.fa"
    GM_stat_df = ana_GM.gtf2csv(GM_file)
    GM_stat_df['seq_len'] = None
    seq_list =  list(SeqIO.parse(reads, "fasta"))
    for record in seq_list:
        #assign the sequence length column.
        if record.description not in GM_stat_df.values: 
            continue
        else:
            GM_stat_df.loc[GM_stat_df['seqid'] == record.description, 'seq_len'] =  len(record.seq)
            #print(len(record.seq), "\tAssigned value: ",  GM_stat_df.loc[GM_stat_df['seqid'] == record.description, 'seq_len'].values)

    GM_stat_df['rank_1'] = GM_stat_df['seq_len'].rank(method = 'dense', ascending = False)
    GM_stat_df['rank_2'] = GM_stat_df['n_frag'].rank(method = 'dense', ascending = True)
    GM_stat_df['rank_3'] = GM_stat_df['max_score'].rank(method = 'dense', ascending = False)
    GM_stat_df['rank_4'] = GM_stat_df['sum_score'].rank(method = 'dense', ascending = False)
    GM_stat_df['MRR_value'] = (1/GM_stat_df['rank_1'] + 1/GM_stat_df['rank_2'] + \
                        1/GM_stat_df['rank_3'] + 1/GM_stat_df['rank_4'])/4
    GM_stat_df['rank_MRR'] = GM_stat_df['MRR_value'].rank(method = 'dense', ascending = False)
    GM_stat_df.to_csv( GM_file + "_MRR.csv")    
    best_seq_id = GM_stat_df.loc[GM_stat_df['rank_MRR'] == 1, 'seqid'].values
    print(best_seq_id, " ranks best.")
    
    for record in seq_list:
        #assign the sequence length column.
        if record.description == best_seq_id: 
            SeqIO.write(record, final_out_file, "fasta")
            

    
    
    