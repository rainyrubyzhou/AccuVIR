import sys
import os
assert("bioinfo" == os.environ['CONDA_DEFAULT_ENV'])
import time
import random
import networkx as nx
import numpy as np
import pandas as pd
from BeamSearchTree import BeamSearchTreeNode
import heapq
import copy


def count_branch(nx_g, thres):
    # Input a path(currently backbone path), a nx_g graph
    # Condition for candicate patch: 
    # 1. begin node: out_degree > 1; weight of second largest ed > threshold(e.g = 3)
    # 2. end node: closest node from the begin node s.t. in_degree > 1 
    import heapq
    candi_loc = []
    topo_l = list(nx.topological_sort(nx_g))
    #print("topo_last: ", topo_l[-1])
    for u in range(1, topo_l[-1]):
        if nx_g.out_degree(u) > 1:
            wei_l = []
            out_l = []
            for u, u_out, attr in nx_g.out_edges(u, data = True):
                wei_l.append(attr['weight'])
                out_l.append(u_out)
            '''
            secon_max = heapq.nlargest(2, wei_l)[1]
            if secon_max > thres * sum(wei_l): 
                candi_loc.append(u)
            '''
            if (sum(wei_l) - max(wei_l)) > thres * sum(wei_l):
                candi_loc.append(u)
            
    #print("num of candi_loc:", len(candi_loc))
    return candi_loc

def assert_candi(nx_g, thres, node):
    # given a graph, a node, a threshold, decide whether it is a candidate node
    flag = False
    if nx_g.out_degree(node) > 1:
        wei_l = []
        out_l = []
        for node, nd_out, attr in nx_g.out_edges(node, data = True):
            wei_l.append(attr['weight'])
            out_l.append(nd_out)
        if (sum(wei_l) - max(wei_l)) > thres * sum(wei_l):
            flag = True #this is a candidate node of given threshold
    return flag

def test_candi_homo(nx_graph, thres, homo_df, topo_list):
    candi_loc = count_branch(nx_graph, thres)
    res = [] 
    for nd in candi_loc:
        for i, row in homo_df.iterrows():
            if topo_list.index(nd) <= topo_list.index(row['homo_end']) and \
                topo_list.index(nd) >= topo_list.index(row['homo_begin']):
                res.append((nd, row))
                break
    #print("# candi :", len(res)) #number of candidate nodes
    #print(res[0:10])
    #return res
#test_candi_homo(nx_graph, 0.35, homo_df) #check how many candidate nodes are in detected homopolyer region
 
def loca_hp(nd, homo_df, topo_list):
    for i, row in homo_df.iterrows():
        if topo_list.index(nd) <= topo_list.index(row['homo_end']) and \
            topo_list.index(nd) >= topo_list.index(row['homo_begin']):
            return row['homo_len']

def sample(nx_graph, graph_out_pref, out_f):
    homo_df = pd.read_csv(graph_out_pref + "_homo_loc.csv")
    topo_list = list(nx.topological_sort(nx_graph))    
    bb_path = list(range(1, topo_list[-1])) #e.g. [0, 1, ..., 9722] as initial backbone path.   
    # generate a bunch of reads for thres_list:
    thres_list = [0.1 + 0.01*i for i in range(41)]  #step up weight from 0.1 to 0.5
    branch_out_log = graph_out_pref + "_sampling.log"
    branch_out_f = out_f
    n_branch = [] #store the number of actual changing sites given different threshold.
    
    with open(branch_out_f, 'w') as out, open(branch_out_log, 'w') as log:   
        for thres in thres_list:
            bb_path = list(range(1, topo_list[-1]))
            candi_on_bb_len_list = []
            patch_list = []
            candi_loc = count_branch(nx_graph, thres)
            patched_bb_list = [] #store multiple pateched paths
            print("Number of initial candidate nodes: ", len(candi_loc), file = log)
            candi_on_bb = []
            for bb_nd in bb_path:
                if bb_nd in candi_loc:
                    candi_on_bb.append(bb_nd)
            #branch for patches:
            candi_rd2 = [] # record candidate nodes for second round branch
            #print("example candidate nodes on backbone:", candi_on_bb[0:10])
            #print("Candidate nodes on backbone path (searching starting sites)", len(candi_on_bb), "\t", candi_on_bb)
            for bb_nd in candi_on_bb:
                patch = [bb_nd]
                patch_len = 1
                # start from a node that diverges with high confidence
                while True:
                    wei_l = []
                    out_l = []
                    cur = patch[-1]
                    if nx_graph.out_degree(cur) == 0: # not the final node
                        break
                    #print("cur: ", cur)
                    for u, u_out, attr in nx_graph.out_edges(cur, data = True):
                        wei_l.append(attr['weight'])
                        out_l.append(u_out)
                    idx = -1 #idx for choosing the next node in out_node list.
                    if patch_len == 1:
                        #iterate the weight list from large to small, to choose the best node not in backbone list
                        #for i,wei in enumerate(heapq.nlargest(2, wei_l)):
                        for i,wei in enumerate(sorted(wei_l, reverse = True)):
                            if out_l[wei_l.index(wei)] not in bb_path: #next node not in backbone
                                idx = wei_l.index(wei)
                                break
                        if idx == -1: #no suitable next node for this candidate
                            break
                    if patch_len > 1: 
                        idx = wei_l.index(max(wei_l))     
                        #for nodes afterwards, if it's in candidate list, record the best for another path
                        """
                        if cur not in candi_loc:
                            idx = wei_l.index(max(wei_l)) #choose the largest node
                        else:
                            candi_rd2.append(cur)                          
                        """
                        #define tolerance using the location of current node (length of homopolymer region)
                        tolerance = 0.5
                        #hp_len = loca_hp(cur, homo_df, topo_list) # locate the length of this homopolymer region
                        hp_len = 3
                        if assert_candi(nx_graph, thres * tolerance, cur):
                            candi_rd2.append(cur) #this node is worth doing second round sampling                        
                    next = out_l[idx]
                    patch.append(next)
                    patch_len += 1
                    #if next <= topo_list[-1] and nx_graph.out_degree(next) == 1:
                    if next <= topo_list[-1]:
                        break
                if len(patch) > 1:
                    patch_list.append(patch)
            patched_bb = bb_path  #initial patched path is the same as backbone path
            cnt = 0 # coubt of patches longer than 2
            #print(thres, "\t", [len(patch) for patch in patch_list])
            for patch in patch_list:
                if len(patch) > 2: #only add the patch longer than before
                    flag = 0
                    cnt += 1
                    if patch[0] in patched_bb and patch[-1] in patched_bb:
                        flag = 1 #This patch is in the output seq
                        idx1 = patched_bb.index(patch[0])
                        idx2 = patched_bb.index(patch[-1])
                        #patched_bb[idx1 - 1 : idx2] = patch
                        patched_bb[idx1 : idx2 + 1] = patch
                    #print(patch, flag) 
            n_branch.append(cnt)
            print("actual changes:", cnt, file = log)
            print("second round sites:", len(candi_rd2), candi_rd2, file = log)
            print("new thres:", thres * tolerance, file = log)
            patched_bb_list.append([str(thres) + "_rd1", patched_bb])

            #second round
            patch_list_rd2 = []
            for bb_nd in candi_rd2:
                candi_rd3 = []
                patch = [bb_nd]
                patch_len = 1
                # start from a node that diverges with high confidence
                while True:
                    wei_l = []
                    out_l = []
                    cur = patch[-1]
                    if nx_graph.out_degree(cur) == 0: # not the final node
                        break
                    #print("cur: ", cur)
                    for u, u_out, attr in nx_graph.out_edges(cur, data = True):
                        wei_l.append(attr['weight'])
                        out_l.append(u_out)
                    idx = -1 #idx for choosing the next node in out_node list.
                    if patch_len == 1:
                        #iterate the weight list from large to small, to choose the best node not in backbone list
                        #for i,wei in enumerate(heapq.nlargest(2, wei_l)):
                        for i,wei in enumerate(sorted(wei_l, reverse = True)):
                            if out_l[wei_l.index(wei)] not in bb_path: #next node not in backbone
                                idx = wei_l.index(wei)
                                break
                        if idx == -1: #no suitable next node for this candidate
                            break
                    if patch_len > 1: 
                        idx = wei_l.index(max(wei_l))     
                        #for nodes afterwards, if it's in candidate list, record the best for another path
                        """
                        if cur not in candi_loc:
                            idx = wei_l.index(max(wei_l)) #choose the largest node
                        else:
                            candi_rd3.append(cur)                          
                        """
                        #define tolerance using the location of current node (length of homopolymer region)
                        tolerance = 0.5
                        if assert_candi(nx_graph, thres * tolerance, cur):
                            candi_rd3.append(cur) #this node is worth doing second round sampling                        
                    next = out_l[idx]
                    patch.append(next)
                    patch_len += 1
                    #if next <= topo_list[-1] and nx_graph.out_degree(next) == 1:
                    if next <= topo_list[-1]:
                        break
                if len(patch) > 1:
                    patch_list_rd2.append(patch)
            patched_bb_2rd = copy.deepcopy(patched_bb)  #initial patched path is the same as patched path in first round"
            cnt = 0 # coubt of patches longer than 2
            #print("patch for rd2:\t", [patch for patch in patch_list_rd2])
            for patch in patch_list_rd2:
                if len(patch) > 2: #only add the patch longer than before
                    flag = 0
                    cnt += 1
                    if patch[0] in patched_bb_2rd and patch[-1] in patched_bb_2rd:
                        flag = 1 #This patch is in the output seq
                        idx1 = patched_bb_2rd.index(patch[0])
                        idx2 = patched_bb_2rd.index(patch[-1])
                        #patched_bb[idx1 - 1 : idx2] = patch
                        print((idx1, idx2), patched_bb_2rd[idx1 : idx2 + 2], len(patched_bb_2rd), file = log)
                        patched_bb_2rd[idx1 : idx2 + 1] = patch
                        print((idx1, idx2), patched_bb_2rd[idx1 : idx2 + 2], len(patched_bb_2rd), file = log)

            n_branch.append(cnt)
            #print("actual changes:", cnt)
            #print("third round sites:", len(candi_rd3), candi_rd3)
            #print(len(patched_bb_2rd), len(patched_bb))        
            patched_bb_list.append([str(thres) + "_rd2", patched_bb_2rd])

            #out put the final
            for i, patched_bb in enumerate(patched_bb_list):
                seq = ''.join([nx_graph.nodes[i]['base'] for i in patched_bb[1]]) #remove the begin and end node
                out.write(">" + patched_bb[0] + '\n') 
                out.write(seq+'\n')
