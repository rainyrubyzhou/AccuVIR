#################################################################################$$
# Author: Runzhou
# Commonly used functions in for graph.
#################################################################################$$
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord 
import aligngraph
import os
import numpy
import networkx as nx

edge_l = []
node_l = []
node_id_dict = {}

def get_nd_no(nd_id):
    """
    map the node id back to a numerical id if the key exists
    """
    try:
        return node_id_dict[nd_id]
    except:
        print("The id of node not in the dictionary")
        return -1
        
def get_nd_id(nd_no):
    try:
        return [id for id, no in node_id_dict.items() if no == nd_no]
    except:
        print("The no of node not in the dictionary")

def write_lines(lines, out_file):
    with open(out_file, 'w') as out_f:
        for line in lines:
            print(line, file = out_f)

def aln2gfa(aln_graph, out_f = None):
    global edge_l, node_l, node_id_dict
    edge_l = list(aln_graph.edges.values())
    node_l = list(aln_graph.nodes.values())
    for i, nd in enumerate(node_l):
        node_id_dict[nd.ID] = i

    lines = []
    lines.append("H	VN:Z:1.0")
    # add all segment lines
    # segment example: "S	2480	TGATAA"
    for i, nd in enumerate(node_l) :
        if nd.base == 'B' or nd.base == 'E':
            nd.base = "AAA"
        line = "S\t"+ str(i+1) + '\t' + str(nd.base)
        lines.append(line)
    # add all link lines
    # link example: "L	5373	+	5374	+	0M"
    for ed in edge_l:
        src_no = get_nd_no(ed.in_node.ID) + 1
        des_no = get_nd_no(ed.out_node.ID) + 1
        line = "L\t" + str(src_no) + '\t+\t' + str(des_no) + '\t+\t0M'
        lines.append(line)
    #print(lines[0:5], lines[-5:])    
    
    out_file = out_f if out_f else "aln_graph.gfa"
    write_lines(lines, out_file)
    if lines:
        os.system("vg view -Fd %s |dot -Tsvg  -o %s"%(out_file, out_file+".svg"))
        return 0
    else:
        return 1

def nx2gfa(nx_graph, out_f = None):
    lines = []
    n_nd = nx_graph.number_of_nodes()
    n_ed = nx_graph.number_of_edges()
    lines.append("V:\t" + str(n_nd) + "\t" + "E:\t" + str(n_ed))
    for nd, attr in nx_graph.nodes(data=True):
        line = "V\t" + str(nd) + '\t' + attr['base']
        lines.append(line)
    for u, v, attr in nx_graph.edges(data=True):
        line = "E\t" + str(u) + '\t' + str(v) + '\t:' + str(attr['weight'])
        lines.append(line)    
    out_file = out_f if out_f else "aln_graph.gfa"
    write_lines(lines, out_file)

def nx_gfa_output(nx_graph, out_f = None):
    # output a nx_graph to gfa format and visualize it.
    for i, nd in enumerate(node_l):
        node_id_dict[nd.ID] = i

    lines = []
    lines.append("H	VN:Z:1.0")
    # add all segment lines
    # segment example: "S	2480	TGATAA"
    for nd, attr in nx_graph.nodes(data = True):
        if attr['base'] == 'B' or attr['base'] == 'E':
            attr['base'] = "AAA"
        line = "S\t"+ str(nd + 1) + '\t' + str(attr['base'])
        lines.append(line)
    # add all link lines
    # link example: "L	5373	+	5374	+	0M"
    for u, v, attr in nx_graph.edges(data = True):
        line = "L\t" + str(u + 1) + '\t+\t' + str(v + 1) + '\t+\t0M'
        lines.append(line)
    #print(lines[0:5], lines[-5:])    
    out_file = out_f if out_f else "aln_graph.gfa"
    write_lines(lines, out_file)
    #os.system("vg view -Fd %s |dot -Tsvg  -o %s"%(out_file, out_file+".svg"))
    return 

def load_nx_graph(graph_f):
    nx_graph = nx.DiGraph()
    with open(graph_f) as graph_f:
        head_line = graph_f.readline()
        n_nd = int(list(head_line.split())[1])
        n_ed = int(list(head_line.split())[3])
        for i in range(n_nd):
            #load nodes 
            line = graph_f.readline().split()
            nd_id = int(line[1])
            nd_base = line[2]
            nx_graph.add_node(nd_id, base = nd_base)
        for i in range(n_ed):
            #load edges
            line = graph_f.readline().split()
            src = int(line[1])
            des = int(line[2])
            wei = int(line[3].split(':')[-1])
            nx_graph.add_edge(src, des, weight = wei)

    return nx_graph



class Res:
    def __init__(self, path):
        if path == []:
            self.path = []
            self.len = 0
        else:
            self.path = path
            self.len = len(path)
class SearchPara:
    def __init__(self, q_str, aln_graph, src, des):
        self.q_str = q_str
        self.aln_graph = aln_graph
        self.begin = src
        self.end = des
        self.res = []
        self.shortest_len = len(self.q_str)
        self.Fullmatch = False
        self.match = False
        self.max_match = 1

max_len = 0
longest_path = Res([])
max_match = 0 

def print_path(path, aln_graph):
    print("path score: ", aln_graph.score(path), '\tlen:',len(path))
    print([get_nd_no(nd.ID)+1 for nd in path])

def search_util(para):
    global max_len, longest_path, all_path, max_match
    q_str =  para.q_str
    candi_nodes = []
    if para.begin!= para.end:
        #prepare the 
        candi_nodes = [ed.out_node for ed in para.begin._out_edges]

    if q_str[0] != para.begin.base:
        return para
        #current node not match, return the same parameter
    if len(all_path) >= max_match:
        return para
        #there're enough matches

    else:
        #current node matched, update and continue
        para.res.append(para.begin)
        if len(q_str) == 1:
            para.Fullmatch = True
            # This is the last base and is matched
            all_path.append(para.res)
            print_path(para.res, para.aln_graph)
            return para

        else:
            para.match = True
            temp_res = []
            for node in candi_nodes:                
                new_para = SearchPara(q_str[1:], para.aln_graph, node, para.aln_graph.end_node)
                new_para.res = para.res
                res_para = search_util(new_para)
                if res_para.match == False:
                    #The next base of q_str does not match next node.
                    #para.res.pop()
                    continue
                else:
                    # Next base is matched with next node, update the returned value
                    # Get the shortest query length (longest matched length)
                    if res_para.shortest_len < para.shortest_len:
                        para.shortest_len = res_para.shortest_len
                    if res_para.Fullmatch == True:
                        para.Fullmatch = res_para.Fullmatch
                    if len(res_para.res) > len(temp_res):
                        temp_res = res_para.res
            #para.res = temp_res   
    if len(para.res) > 0:
        para.res.pop()
    if max_len < len(para.res):
        longest_path = Res(para.res)
        max_len = len(para.res)

    return para    
        
def search(q_str, aln_graph, k=5, m_match = 3):    

    print("new")
    global edge_l, node_l, node_id_dict, max_match
    max_match = m_match
    edge_l = list(aln_graph.edges.values())
    node_l = list(aln_graph.nodes.values())
    for i, nd in enumerate(node_l):
        node_id_dict[nd.ID] = i
    global all_path
    #iterate each node in the graph as a possible search beginning
    
    # generate some blocks to see whether they exist in the graph.
    # 1. shatter the reference into k pieces
    graph_loc = 0
    for j in range(k):
        all_path = []
        s = j * (int(len(q_str)/k))
        sub_q_str = q_str[s : s + int(len(q_str)/k)]
        print("\nquery:[%s,%s]"%(s, s + int(len(q_str)/k)))
        for i, node in enumerate(aln_graph.get_sorted_nodes()):
            if i < graph_loc:
                continue
            search_para = SearchPara(sub_q_str, aln_graph, node, aln_graph.end_node)
            search_res = search_util(search_para)
        
        if len(all_path) > 0:
            graph_loc += len(sub_q_str)

        #print("# matches for this query:%s"%(len(all_path)))
        #print("Final match flag:%s, shortest_query:%s"%(search_res.Fullmatch, search_res.shortest_len))
        if len(all_path) == 0:
            print("Longest exact match:", len(longest_path.path))
            print_path(longest_path.path, aln_graph)

def nx_info(nx_g):
    # print the basic information of a networkx graph
    print("# nodes", len(nx_g.nodes))
    print("# edges", len(nx_g.edges))

def aln2nx(aln_graph): 
    global edge_l, node_l, node_id_dict
    node_l = list(aln_graph.nodes.values())
    edge_l = list(aln_graph.edges.values())

    nx_g = nx.DiGraph()
    #dict stores the id and the number id of each node
    node_id_dict = {}
    for i, nd in enumerate(node_l):
        node_id_dict[nd.ID] = i
        nx_g.add_node(i, base = nd.base)
        if nd.base == 'B':
            src = i
        elif nd.base == 'E':
            des = i      
    for i, edge in enumerate(edge_l):
        #src = node_id_dict[edge.in_node.ID]
        #des = node_id_dict[edge.out_node.ID]
        src = get_nd_no(edge.in_node.ID)
        des = get_nd_no(edge.out_node.ID)
        if src == -1 or des == -1:
            print(i," missing")
        nx_g.add_edge(src, des, weight = edge.count)
    return [nx_g, src, des]

def sort_nx(nx_graph):
    # input a nx_graph and transform it to sorted(map nodes id)
    topo_list = list(nx.topological_sort(nx_graph))
    nd_order_dict = {}
    nx_sorted = nx.DiGraph()
    if topo_list[0] != -1:        
        for i, old_id in enumerate(topo_list):
            nd_order_dict[old_id] = i + 1
            nx_sorted.add_node(i + 1, base = nx_graph.nodes[old_id]['base'])
    '''elif topo_list[0] == 1:        
        for i, old_id in enumerate(topo_list):
            nd_order_dict[old_id] = i
            nx_sorted.add_node(i, base = nx_graph.nodes[old_id]['base'])'''
    # for every edge in the original graph, add it to new graph using sorted ID;
    for u, v, attr in nx_graph.edges(data=True):
        # (u,v): edge pair; attr: attribute dictionary
        nx_sorted.add_edge(nd_order_dict[u], nd_order_dict[v],  weight = attr['weight'])        
    return nx_sorted

def nx_path_score(nx_graph, path):
    # input: a nx graph with edge property 'weight'; a path
    # return: sum weight of this path
    def path2iter(path):
        for i in range(len(path)-1):
            yield (path[i], path[i+1])
    ed_l = path2iter(path)
    score = sum(nx_graph.get_edge_data(*ed)['weight'] for ed in ed_l )
    return score
def nx_homo(nx_g, out_f):
    #count the homopolymer region and return a list of range
    homo_loc = []
    min_len = 3
    win_len = 10
    topo_list = list(nx.topological_sort(nx_g)) 
    #print(topo_list[0:30])
    for i, nd in enumerate(topo_list[10:]):
        #nd_l = [u for u,v, attr  in nx_g.in_edges(nd, data = True)]
        nd_l = []
        for u,v in nx_g.in_edges(nd):
            nd_l.append(u)
        in_base_l = [nx_g.nodes[nd]['base'].upper() for nd in nd_l]
        char_l = ['A', 'C', 'G', 'T']
        freq_l = []
        freq_l.append(in_base_l.count('A'))
        freq_l.append(in_base_l.count('C'))
        freq_l.append(in_base_l.count('G'))
        freq_l.append(in_base_l.count('T'))
        
        if max(freq_l) >= min_len: #it's a homopolymer region of length min_len
            print(i, nd, nd_l, freq_l)
            homo_end = nd
            homo_begin = 0
            homo_char = char_l[freq_l.index(max(freq_l))] #base of this hp region
            #now locate the beginning of this region
            #win = topo_list[ i - win_len + 1: i + 1]
            win = topo_list[ 10 + i - win_len + 1: 10 + i] #i is the index start from 10 in topo-list
            
            for win_nd in win[::-1]:
                if nx_g.nodes[win_nd]['base'].upper() != homo_char:
                    homo_begin = win_nd
                    homo_len = topo_list.index(homo_end) - topo_list.index(homo_begin) - 1 #exclude the first base that's not in homopolymer region
                    break
            homo_loc.append((homo_begin, homo_end, homo_char, homo_len))
    #print("homopolymer regions:", len(homo_loc),"\n", homo_loc)    
    from matplotlib import pyplot as plt
    import pandas as pd
    homo_df = pd.DataFrame(homo_loc, columns = ['homo_begin', 'homo_end', 'homo_char', 'homo_len'])
    homo_df.hist('homo_len')
    plt.show()
    homo_df.to_csv(out_f)

    return homo_df
