#################################################################################$$
# Author: Runzhou
# Driver Program for DBS and Branched Sampling. 
# For iteratively constructing the graph, users need to mannually pass in the new backbone sequence.
#################################################################################$$
import sys
import time
import os
import random
import argparse
import utils
import beam_diverse as DBS
import branch_sample as BS

if __name__ == '__main__':
    #parse args
    BeamWidth = 500
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', type=str, required=True, help = "Reads file for graph construction(in fasta format).")
    parser.add_argument('-b', type=str, required=True, help = "Backbone sequence file for graph construction (in fasta format).")
    parser.add_argument('-m', type=str, required=False, default = 3, help = "Select mode for the path searching: '1' for diverse beam search; \
                                                                                                    '2' for branched sampling; \
                                                                                                    '3' for both search module. \
                                                               ('1' is recommended for first round search.)")
    parser.add_argument('--beamwidth', type=int, default = 500,  required=False, help = "Beamwidth for diverse beam search (default: 500).")
    
    args = parser.parse_args()

    print(args.r, args.b)
    ec_reads = args.r
    backbone = args.b
    mode = int(args.m)
    BeamWidth = args.beamwidth    
    
    graph_out_pref = (ec_reads + "_ON_" +\
                     (os.path.splitext(os.path.basename(backbone))[0]) )
    time1 = time.time()
    aln_graph = utils.construct_aln_graph_from_fasta(ec_reads,  backbone, 5000)
    time2 = time.time()
    print("Finish building graph in %s s"%(time2 - time1))
    # generate corresponding dagcon seq.
    aln_graph.merge_nodes()
    aln_graph.generate_consensus()
    consen_str = aln_graph.consensus_str
    consen_score = aln_graph.score(aln_graph.consensus_path)
    with open(graph_out_pref + "_consen.fa", 'w') as consen_f:
        consen_f.write(">consen_" + graph_out_pref + ': ' + str(consen_score) + '\n')
        consen_f.write(consen_str + '\n')
    nx_pack = utils.aln2nx(aln_graph)        
    utils.nx2gfa(nx_pack[0], graph_out_pref + '.graph')
    if mode == 1:
        print("Running mode: only DBS")
        print("Running DBS of beam width: ", BeamWidth)
        DBS.Grouped_beam(nx_pack[0], BeamWidth, graph_out_pref + "_DBS_" + str(BeamWidth) + ".fa")
        
        utils.extract_longest(graph_out_pref + "_DBS_" + str(BeamWidth) + ".fa", graph_out_pref + "_DBS_" + str(BeamWidth) + "_longest.fa")
    if mode == 2:
        #locate all homopolymer regions
        print("Running mode: only sampling")
        print("Scaning homopolymer regions...")
        utils.nx_homo(nx_pack[0], graph_out_pref + "_homo_loc.csv")
        print("Branched sampling...")
        BS.sample(nx_pack[0], graph_out_pref, graph_out_pref + "_sampling.fa")
    if mode == 3:
        print("Running mode: both search module")
        print("Running DBS of beadm width: ", BeamWidth)
        DBS.Grouped_beam(nx_pack[0], BeamWidth, graph_out_pref + "_DBS_" + str(BeamWidth) + ".fa")
    
        print("Scaning homopolymer regions...")
        utils.nx_homo(nx_pack[0], graph_out_pref + "_homo_loc.csv")
        print("Branched sampling...")
        BS.sample(nx_pack[0], graph_out_pref, graph_out_pref + "_sampling.fa")
        # merge the output of two modules and filter out short ones for intermediate output(mainly from DBS).
        #utils.merge_filter(graph_out_pref + "_DBS_" + str(BeamWidth) + ".fa", graph_out_pref + "_sampling.fa", graph_out_pref + "_filtered.fa")
        utils.merge_filter(graph_out_pref + "_DBS_" + str(BeamWidth) + ".fa", graph_out_pref + "_sampling.fa", graph_out_pref + "_merge.fa")
    
    
