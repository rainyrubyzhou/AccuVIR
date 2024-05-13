#################################################################################$$
# Author: Runzhou
# Commonly used functions in for the whole pipeline.
# Graph building related utils are from pgdagcon's repo.
#################################################################################$$
from math import log, sqrt
import os
import pandas as pd
import aligngraph
from aligngraph import convert_mismatches, AlnGraph
import networkx as nx
from Bio import SeqIO


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
            #print(i, nd, nd_l, freq_l)
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
    homo_df = pd.DataFrame(homo_loc, columns = ['homo_begin', 'homo_end', 'homo_char', 'homo_len'])
    
    homo_df.to_csv(out_f)

    return homo_df
    
class AlignGraphUtilError(Exception):
    pass

def phi_coeff(xvec, yvec):
    nx1 = sum( xvec == 1 )
    nx0 = len(xvec) - nx1
    ny1 = sum( yvec == 1 )
    ny0 = len(yvec) - ny1
    xyvec = (xvec << 1) + yvec
    n11 = sum( xyvec == 3 )
    n10 = sum( xyvec == 2 )
    n01 = sum( xyvec == 1 )
    n00 = sum( xyvec == 0 )
    return (n11 * n00 - n10 * n01) / sqrt( nx1 * nx0 * ny1 * nx0 + 1)

class Simple_Alignment_Hit(object):
    """ A simple class to wrap the output of the blasr "-m 5" option """
    def __init__(self, rm5_line):
        rm5_line = rm5_line.strip().split()
        self.query_id = rm5_line[0]
        self.query_length = int(rm5_line[1])
        self.query_start = int(rm5_line[2])
        self.query_end = int(rm5_line[3])
        self.query_strand = rm5_line[4]
        self.target_id = rm5_line[5]
        self.target_length = int(rm5_line[6])
        self.target_start = int(rm5_line[7])
        self.target_end = int(rm5_line[8])
        self.target_strand = rm5_line[9]
        self.alignedQuery = rm5_line[16]
        self.alignedTarget = rm5_line[18]

def simple_align_hit_iterator(rm1_fn, ref_group = None):
    with open(rm1_fn) as f:
        for l in f:
            ll = l.strip().split()
            if ref_group != None and ll[5] != ref_group:
                continue
            yield Simple_Alignment_Hit(l)

def get_aln_array(aln_iter,
                  max_num_reads = None, 
                  remove_in_del = False, 
                  min_length = None):

    rMap = dict(zip("ACGTacgtNn-", "TGCAtgcaNn-"))

    alns = []
    nread = 0
    i = 0
    bp = {}
    reads = {}

    for h in aln_iter:

        if min_length != None:
            if h.target_end - h.target_start < min_length:
                continue

        nread += 1

        if max_num_reads != None:
            if nread > max_num_reads:
                break

        read_id = h.query_id 

        if read_id in reads: continue #avoid duplicate readId in the data
        
        ts = h.target_start
        te = h.target_end
        qs = h.query_start
        qe = h.query_end

        alnT, alnQ = h.alignedTarget.upper(), h.alignedQuery.upper()
        dqPos = 1
        
        if h.target_strand == '-':
            alnT = "".join([rMap[c] for c in alnT[::-1]])
            alnQ = "".join([rMap[c] for c in alnQ[::-1]])
            dqPos = -1
            qs, qe = qe - 1, qs - 1

        if remove_in_del:
            aln_pair = zip(alnQ, alnT)
            new_aln_pair = []
            for b1, b2 in aln_pair:
                if b1 == "-":
                    b1 = b2
                if b2 != "-":
                    new_aln_pair.append( (b1,b2) )
            alnQ, alnT = zip(*new_aln_pair)
            alnQ = "".join(alnQ)
            alnT = "".join(alnT)

        alnQ, alnT = convert_mismatches(alnQ,alnT)

        if alnQ[0] == "-" or alnT[0] == "-":
            continue
        if alnQ[-1] == "-" or alnT[-1] == "-":
            continue

        #print h.target_strand 
        #print "r %10d %10d" % (qs, qe), alnQ
        #print "t %10d %10d" % (ts, te), alnT
        #print

        bp[read_id] = ( (qs, qe), ts)
        alns.append(  ( ( qs, qe, alnQ ), ( ts, te, alnT ) , read_id ) )
        tPos = ts
        qPos = qs
        assert len(alnQ) == len(alnT)
        reads[read_id] = alnQ.replace("-","")
        i += 1

    return alns

def construct_aln_graph_from_fasta(read_fasta_fn, 
                                    backbone_fasta_fn, 
                                    max_num_reads = None, 
                                    remove_in_del = False, 
                                    ref_group = None,
                                    min_length = None):

    import os

    os.system("blasr %s %s -m 5 --nproc 16 --out %s" % (read_fasta_fn, backbone_fasta_fn, read_fasta_fn+".aln_unsorted"))
    os.system("cat %s | sort > %s" % (read_fasta_fn+".aln_unsorted", read_fasta_fn+".aln" ))
    
    aln_hit_iterator = simple_align_hit_iterator(read_fasta_fn+".aln", ref_group = ref_group)

    alns = get_aln_array(aln_hit_iterator, 
                         max_num_reads = max_num_reads, 
                         remove_in_del = remove_in_del, 
                         min_length = min_length)
        
    backboneSeq = open(backbone_fasta_fn).read()
    backboneSeq = "".join(backboneSeq.split('\n')[1:])

    g = AlnGraph(backboneSeq)

    i = 0
    for aln in alns:
        rId = aln[2]
        rId = rId.split("/")[0]
        aln = aln[0:2]
        g.add_alignment( aln, "%s" % rId)

    return g

def sorted_nodes(g):
    return g.get_sorted_nodes()

def read_node_vector(g, entropy_th = 2.5):
    read_to_nodes, high_entropy_nodes = g.get_read_node_vector(entropy_th=entropy_th)
    return read_to_nodes, high_entropy_nodes

def clustering_read( read_to_nodes, high_entropy_nodes, k_cluster = 2, random_seed = 42, cleanup_th = 0.5):
    
    import random
    random.seed(random_seed)

    cluster = {}
    read_to_binary_vector = {}
    count = 0
    for r in read_to_nodes:   
        read_to_binary_vector[r] = np.array([ 1 if c != "-" else -1 for c in read_to_nodes[r] ])
        k = random.randint(0, k_cluster-1)
        #k = count % k_cluster
        cluster.setdefault(k,[])
        cluster[k].append(r)
        count += 1
    
    n_iteration = 0
    cluster_vec = {}
    while 1:
        for k in range(k_cluster):
            new_vec = np.zeros(len(high_entropy_nodes),dtype=np.float)
            #print k, len(cluster[k])
            for r in cluster[k]:
                new_vec += read_to_binary_vector[r]
            new_vec /=  (len(cluster[k])+1)
            cluster_vec[k] = np.array([1 if v>=0 else -1 for v in new_vec])

        #avoid degnerated cluster    
        for k in range(k_cluster):
            for j in range(k+1, k_cluster):
                if sum(cluster_vec[j] == cluster_vec[k]) > 0.5*len(high_entropy_nodes):
                    cluster_vec[j] = np.array([random.choice([-1,1]) for v in new_vec])

        cluster = {}
        for r in read_to_nodes:   
            distances = []
            for k in range(k_cluster):
                cluster.setdefault(k,[])
                distances.append( (sum(read_to_binary_vector[r] * cluster_vec[k]), k) )

            distances.sort()

            cluster[distances[-1][1]].append(r)
        
        n_iteration += 1
        if n_iteration  > 10:
            break
            
    cluster = {}
    for r in read_to_nodes:   
        distances = []
        for k in range(k_cluster):
            cluster.setdefault(k,[])
            distances.append( (sum(read_to_binary_vector[r] * cluster_vec[k]), k) )

        distances.sort()
        if distances[-1][0] > cleanup_th*len(high_entropy_nodes):
            cluster[distances[-1][1]].append(r)
    return cluster, cluster_vec

def sorted_node_data(aln_graph, entropy_th = 0.5, interval = None):
    ne, hne = aln_graph.get_high_entropy_nodes(coverage_th=0)
    node_to_entropy = dict( [ (v[1],v[2]) for v in ne ] ) 
    read_ids = set()
    for n in aln_graph.nodes.values():
        for r in n.info:
            read_ids.add(r)
    read_id_to_pos = dict(( (x[1],x[0]) for x in enumerate(list(read_ids))) )

    backbone_node_to_pos =  aln_graph.backbone_node_to_pos 

    data = []
    for n in sorted_nodes(aln_graph):
        s = [" "] * len(read_id_to_pos)
        for r in n.info:
            s[read_id_to_pos[r]] = n.base

        if n.base not in ["B", "E"]:
            bpos = backbone_node_to_pos[n.backbone_node]
            if interval != None and (bpos < interval[0] or bpos > interval[1]):
                continue
        ent = node_to_entropy[n] if n in node_to_entropy else 0
        if n.base not in ["B","E"] and ent >= entropy_th:
            data.append ( (  backbone_node_to_pos[n.backbone_node],\
                  "+" if n in aln_graph.consensus_path else "-",\
                  "+" if n.is_backbone == True else "-",\
                  n.base, "".join(s), len(n.info),\
                  n.backbone_node.coverage,\
                  node_to_entropy[n] if n in node_to_entropy else "-" ) )
    return data

def detect_missing(aln_graph, entropy_th = 0.66, interval = None):
    data = sorted_node_data(aln_graph, entropy_th = 0, interval = interval)
    s = []
    for d in data:
        if d[7] > entropy_th and d[2] == "-":
            s.append(d[3].lower())
            continue
        if d[2] =="+" and d[7] > entropy_th:
            #pass
            s.append(d[3].lower())
        elif d[2] == "+":
            s.append(d[3].upper())
    return "".join(s)

def mark_lower_case_base(aln_graph, entropy_th = 0.66, interval = None):
    data = sorted_node_data(aln_graph, entropy_th = 0, interval = interval)
    s = []
    for d in data:
        #if d[7] > entropy_th and d[2] == "-":
        #    s.append(d[3].lower())
        #    continue
        if d[2] =="+" and d[7] > entropy_th:
            #pass
            s.append(d[3].lower())
        elif d[2] == "+":
            s.append(d[3].upper())
    return "".join(s)
    
def merge_filter(dbs_file, samp_file, filter_file):
    # merge the output of two modules and filter out short ones (mainly from DBS).
    max_len = 0
    thres = 100
    for record in SeqIO.parse(dbs_file, "fasta"):
        if len(record.seq) > max_len:
            max_len = len(record.seq)
    for record in SeqIO.parse(samp_file, "fasta"):
        if len(record.seq) > max_len:
            max_len = len(record.seq)

    #print("longest seq = ", max_len)
    filtered_seqs = []
    for record in SeqIO.parse(dbs_file, "fasta"):
        if len(record.seq) > max_len - thres:
            filtered_seqs.append(record)
    for record in SeqIO.parse(samp_file, "fasta"):
        if len(record.seq) > max_len - thres:
            filtered_seqs.append(record)
    SeqIO.write(filtered_seqs, filter_file, "fasta")

    return

def merge(dbs_file, samp_file, merge_file):
    import shutil
    with open(merge_file,'w') as merge_f:
        for f in [dbs_file, samp_file]:
            with open(f,'r') as fd:
                shutil.copyfileobj(fd, merge_f)

def get_median_length(sequences):
    lengths = [len(seq) for seq in sequences]
    sorted_lengths = sorted(lengths)
    mid = len(sorted_lengths) // 2
    if len(sorted_lengths) % 2 == 0:
        median_length = (sorted_lengths[mid - 1] + sorted_lengths[mid]) / 2.0
    else:
        median_length = sorted_lengths[mid]
    return median_length

def extract_median(input_file, output_file):
    sequences = list(SeqIO.parse(input_file, 'fasta'))
    median_length = get_median_length(sequences)
    sequences_with_median_length = [seq for seq in sequences if len(seq) == median_length]
    print("Median_length of DBS results:", median_length, "# sequences of median length: ", len(sequences_with_median_length))
    with open(output_file, 'w') as f:
        SeqIO.write(sequences_with_median_length[0], f, 'fasta')

def get_longest_length(sequences):
    lengths = [len(seq) for seq in sequences]
    sorted_lengths = sorted(lengths)
    longest_length = sorted_lengths[-1]
    return longest_length

def extract_longest(input_file, output_file):
    sequences = list(SeqIO.parse(input_file, 'fasta'))
    longest_length = get_longest_length(sequences)
    sequences_with_longest_length = [seq for seq in sequences if len(seq) == longest_length]
    print("Longest_length of DBS results:", longest_length, "# sequences of longest length: ", len(sequences_with_longest_length))
    with open(output_file, 'w') as f:
        SeqIO.write(sequences_with_longest_length[0], f, 'fasta')    
    