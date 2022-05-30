#################################################################################$$
# Author: Runzhou
# Using Grouped_beam search to make beam search more diverse.
#################################################################################$$
import sys
import time
import os
import random
import beam_utils
from networkx.algorithms.smetric import s_metric
from BeamSearchTree import BeamSearchTreeNode
#assert("bioinfo" == os.environ['CONDA_DEFAULT_ENV'])
sys.path.append("..")
import networkx as nx
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import utils


class BeamSearch:
    def __init__(self, graph, src, des, wid_k, metric) -> None:
        super().__init__()
        self.graph = graph
        self.src = src
        self.des = des
        self.beam_width = wid_k
        self.metric = 2 #metric = 1:"local"; 2:"global without filtering"; 3. "global with filtering"

    def call_beam(self):
        # build BeamSearching Tree and generate corresponding paths
        tree, leaf_l = self.beam_search_tree(self.src, self.metric)
        path_grt = beam_utils.gen_path_in_tree(tree, leaf_l, self.beam_width)
        #print(tree)
        return path_grt

    def map_cut(self, index, cut_list):
        # tranform the selected index(among all candicates of next layer) into selected nodes in cur_layer
        # input: 0/1/2/5 [1,4,6]
        # output:  (range No.3)
        for i, ceiling in enumerate(cut_list):
            if index < ceiling:
                 return i-1

    def beam_search_tree(self, root, metric):
        # input: a graph, src of search
        # return: root node of a searching graph
        # bfs: 
        #root = BeamSearchTreeNode(self.src)
        root = BeamSearchTreeNode(node_id = self.src)
        leaf_list = []  
        node_layer_que = []
        layer = 0
        des_flag = self.beam_width
        topo_list = list(nx.topological_sort(self.graph))
        backbone_len = topo_list[-1] - topo_list[0]
        # set this threshold to keep path longer than (thres*backbone_len)
        len_thres = 0.8
        len_list = []

        node_layer_que.append([self.src, layer])
        time_1 = time.time()
        cur_layer_tree_nd_l = [root] #store the updated nodes of a layer in every loop;
        cur_layer = 0
        
        if metric != 1:
            while des_flag > 0 and cur_layer_tree_nd_l:
                #print("flag:", len(node_layer_que), des_flag)
                #cur_layer = node_layer_que[0][1]
                cur_layer += 1                    
                next_layer_nd_l = []
                next_layer_ed_l = []
                next_layer_weight_l = []
                next_layer_sum_l = []
                cut_list = []  
                count = 0
                # pop all nodes from the same layer
                #print(node_layer_que, cur_layer)
                for tree_nd in cur_layer_tree_nd_l:
                    cut_list.append(count)
                    for ed in list(self.graph.edges(tree_nd.node_id())):
                        # if the topological order difference > 500, ignore these edge
                        #if (self.graph.get_edge_data(*ed) - tree_nd.node_id()) > 500:
                        #    continue
                        local_wei = self.graph.get_edge_data(*ed)['weight']
                        """
                        if metric == 2:
                            if ed[1] == self.des and (cur_layer < len_thres * backbone_len):
                                print(ed, cur_layer)
                                continue
                        """
                        next_layer_ed_l.append(ed)
                        next_layer_nd_l.append(ed[1])
                        next_layer_weight_l.append(local_wei)
                        count += 1
                
                cut_list.append(count) #append the total count of nodes for next layer
                for i in range(len(next_layer_nd_l)):
                    paren_idx = self.map_cut(i, cut_list)
                    sum_wei = cur_layer_tree_nd_l[paren_idx].sum_weight() + next_layer_weight_l[i]
                    next_layer_sum_l.append(sum_wei)
                #selected_nd_l =  random.choices(next_layer_nd_l, k = 2)
                if self.beam_width < len(next_layer_nd_l):
                    """selecting by local weight distribution"""
                    local_selected_indices = beam_utils.trunc_best_k(np.array(next_layer_weight_l), self.beam_width)
                    """selecting by local weight distribution"""
                    sum_selected_indices = beam_utils.trunc_best_k(np.array(next_layer_sum_l), self.beam_width)
                    #selected_indices = list(set(local_selected_indices) | set(sum_selected_indices))
                    selected_indices = local_selected_indices
                    
                else:
                    # if the beam width > # candidate nodes of next layer, return all the indices
                    selected_indices = range(len(next_layer_nd_l))
                new_layer_tree_nd_l = []
                for index in selected_indices:
                    #add node according to the selected nodes:
                    node_layer_que.append([next_layer_nd_l[index], cur_layer + 1])
                    
                    # tell from the index of the node, to decide the parent node
                    paren_idx = self.map_cut(index, cut_list) 
                    #print(index, cut_list, paren_idx, len(cur_layer_tree_nd_l))
                    parent = cur_layer_tree_nd_l[paren_idx] #choose the parent node to append the child

                    child = BeamSearchTreeNode(node_id = next_layer_nd_l[index], 
                                            local_weight = next_layer_weight_l[index],
                                            sum_weight = next_layer_sum_l[index],
                                            parent_nd = parent)                                    
                    new_layer_tree_nd_l.append(child) #prepare the tree node list for next layer
                    parent.add_child(child)
                    
                    
                    if next_layer_nd_l[index] == self.des: 
                        if metric == 2:
                            # do not filter the length when searching
                            len_list.append(cur_layer)
                            des_flag -= 1 
                            self.beam_width -= 1                        
                            leaf_list.append(child)
                        if metric == 3:
                            # filter the length
                            if (cur_layer >= len_thres * backbone_len):
                                len_list.append(cur_layer)
                                des_flag -= 1 
                                self.beam_width -= 1
                                leaf_list.append(child)
                            """2021.07.21: Limiting length here will not work"""
                #update the layer_tree_node list with selected nodes
                cur_layer_tree_nd_l = new_layer_tree_nd_l
                
                #print("\ncur_layer:",[nd_lay_pair[0] for nd_lay_pair in cur_layer_nd_l])
                #print("candi:", next_layer_ed_l) 
                #print("weigh:", next_layer_weight_l) 
                #print("edges for next layer", [next_layer_ed_l[index] for index in selected_indices]) 
                
        time_2 = time.time()
        #print("Building the searching tree: %s s"%(time_2 - time_1))
        #print("path lengths:", len_list)
        return root, leaf_list


#################################################################################$$
# Function1: 
# spliting the search into G groups at beginning based on content(first base),
# search k/G paths in each group.
# do not combine the dissimilarity of a path with other groups 
# input: 
#def Grouped_beam(graph = nx_g, width = beam_wid, out_f = path_out_f):
def Grouped_beam(nx_g, beam_wid, path_out_f):
    nx_graph = nx_g
    beam_wid = beam_wid
    out_f = path_out_f 
    log_f = out_f.split(".")[0] + ".log"
    # 1. choose 4 groups(starting from A, C, G, T respectively)
    # 2. choose N groups(N = number of out nodes from beginnning node)
    topo_list = list(nx.topological_sort(nx_graph))
    with open(log_f, 'w+') as log:
        print("src and des:", topo_list[0], topo_list[-1], file = log)
        #src_group = [ed[1] for ed in nx_graph.edges(topo_list[0])] # divide the groups by out_edges of beginning node
        src_group = [ed[1] for ed in nx_graph.edges(0)] 
        g_group = len(src_group)
        des_nd = topo_list[-1]
        path_grt_list = []
        
        #src_group = [1]
        cnt = len(src_group)
        for src_nd in src_group: 
            #print("weight of this group beginning", nx_graph.edges[topo_list[0], src_nd]['weight']    )
            #if nx_graph.edges[topo_list[0], src_nd]['weight'] < 1 or src_nd >= 500:
            if nx_graph.edges[0, src_nd]['weight'] < 1 or src_nd >= 500:
                cnt -= 1
        print("totol out_nodes: %s;valid group:%s"%(len(src_group), cnt), file = log)
        sub_width = int(beam_wid/cnt)
        if sub_width == 0:
            sub_width = 20 #if the beam_wid < g_group, make sure we have at least 1 path in each group
        #limit the number of group
        
        print("number of groups(from src):", len(src_group), "k in each=", sub_width, file = log)
        
        for src_nd in src_group:    
            #if nx_graph.edges[topo_list[0], src_nd]['weight'] < 1 or src_nd >= 2000:
            if nx_graph.edges[0, src_nd]['weight'] < 1 or src_nd >= 2000:
                #this does not fit the condition to be a group
                continue
            search = BeamSearch(nx_graph, src = src_nd, des = des_nd, wid_k = sub_width, metric = 3)
            path_grt = search.call_beam()
            path_grt_list.append(path_grt)
    with open(out_f, 'w') as out, open(log_f, 'w+') as log:
        # transform the nodes into k sequences
        for i, path_grt in enumerate(path_grt_list):
            for j, path in enumerate(path_grt):
                #seq = ''.join([nx_graph.nodes[i]['base'] for i in path[::-1]])
                seq = ''.join([nx_graph.nodes[i]['base'] for i in reversed(path[1:-1])]) #remove the begin and end node
                score = utils.nx_path_score(nx_graph, path[::-1])\
                
                out.write(">G%s_k%s: %s"%(i+1,j+1, score)+'\n') # the sequence id ">G1_k2: 30000" stands the 2-nd path in 1-st group with score 30000.
                out.write(seq+'\n')             
    return 

# Function2: 
# wrap general BeamSearch(build the searching tree from src)
def Normal_beam(nx_g, beam_wid, path_out_f):
    nx_graph = nx_g
    beam_wid = beam_wid
    out_f = path_out_f 
    
    topo_list = list(nx.topological_sort(nx_graph))
    print("src and des:", topo_list[0], topo_list[-1])
    src_group = [ed[1] for ed in nx_graph.edges(topo_list[0])] # divide the groups by out_edges of beginning node
    
    src_nd = topo_list[0]
    des_nd = topo_list[-1]
    path_grt_list = []
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
    search = BeamSearch(nx_graph, src = src_nd, des = des_nd, wid_k = beam_wid, metric = 2)
    
    #only search the first half
    #search = BeamSearch(nx_graph, src = 0, des = 8360, wid_k = beam_wid, metric = 2)
    path_grt = search.call_beam()
    path_grt_list.append(path_grt)
    with open(out_f, 'w') as out:
        # transform the nodes into k sequences
        for i, path_grt in enumerate(path_grt_list):
            for j, path in enumerate(path_grt):
                #seq = ''.join([nx_graph.nodes[i]['base'] for i in path[::-1]])
                seq = ''.join([nx_graph.nodes[i]['base'] for i in reversed(path[1:-1])]) #remove the begin and end node
                score = utils.nx_path_score(nx_graph, path[::-1])\
                
                out.write(">G%s_k%s: %s"%(i+1,j+1, score)+'\n') # the sequence id ">G1_k2: 30000" stands the 2-nd path in 1-st group with score 30000.
                out.write(seq+'\n')
                #path_out.write(">G%s_k%s: %s"%(i+1,j+1, score)+'\n' + "\n")
                #path_out.write(str(path[::-1]) + "\n")                
    return 
