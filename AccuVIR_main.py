#run Diverse Beam Search and Branched Sampling on a reads set.
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
    parser.add_argument('-r', type=str, required=True)
    parser.add_argument('-b', type=str, required=True)
    parser.add_argument('--beamwidth', type=int, required=True)
    args = parser.parse_args()

    print(args.r, args.b)
    ec_reads = args.r
    backbone = args.b
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
    print("Running DBS of beadm width: ", BeamWidth)
    #DBS.Grouped_beam(nx_pack[0], BeamWidth, graph_out_pref + "_DBS_" + str(BeamWidth) + ".fa")
    
    #locate all homopolymer regions
    utils.nx_homo(nx_pack[0], graph_out_pref + "_homo_loc.csv")
    BS.sample(nx_pack[0], graph_out_pref)
    
    
