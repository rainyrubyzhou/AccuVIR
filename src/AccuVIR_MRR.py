#################################################################################$$
# Author: Runzhou
# Function: 
# After running GM, this code analyse the GTF outputs from GM and apply MRR to select best seq.
#################################################################################$$
import argparse
import ana_GM

if __name__ == '__main__':
    #parse args
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', type=str, required=True, help = "Reads output for further ranking(in fasta format).")
    #parser.add_argument('--GM', type=str, required=False, help = "Beamwidth for diverse beam search (default: 500).")
    args = parser.parse_args()
    print(args.r)
    reads = args.r
    GM_file = reads + ".gtf"
    GM_stat_df = ana_GM.gtf2csv(GM_file)

    #MRR: 