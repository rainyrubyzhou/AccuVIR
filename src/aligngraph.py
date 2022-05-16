import re
import textwrap
import collections
import time
import operator
from random import shuffle

rev_table = {"A": "T", "T": "A", "C": "G", "G": "C", "-": "-"}


def reverse_graph(aln_graph):
    """
    reverse the graph
    :param aln_graph: graph we reverse
    :return:
    """
    graph_list = aln_graph.get_sorted_nodes()

    origin_nodes_to_edge = aln_graph.nodes_to_edge
    in_node_list = {}
    for node in graph_list:

        if node.base == "B":
            node.base = "E"

        elif node.base == "E":
            node.base = "B"

        else:
            node.base = rev_table[node.base]

        in_node_list[node] = []
        for in_edge in node._in_edges:
            in_node = in_edge.in_node
            in_node_list[node].append(in_node)

    del_list = []
    for edge in aln_graph.edges:
        del_list.append(aln_graph.edges[edge])

    for edge in del_list:
        aln_graph.delete_edge(edge)

    for node in graph_list:

        for in_node in in_node_list[node]:
            edge_add = AlnEdge(node, in_node)
            edge_add.count = origin_nodes_to_edge[(in_node.ID, node.ID)].count
            edge_add.score = origin_nodes_to_edge[(in_node.ID, node.ID)].score
            aln_graph.add_edge(edge_add)
    
    node_temp = aln_graph.begin_node
    aln_graph.begin_node = aln_graph.end_node
    aln_graph.end_node = node_temp
    
        
    rev_back_bone_nodes = {}
    forward_keys = list(aln_graph.backbone_nodes.keys())        
    for node_index, index in enumerate(reversed(forward_keys)):
        new_key = forward_keys[index]
        rev_back_bone_nodes[new_key] = aln_graph.backbone_nodes[node_index]
            
    aln_graph.backbone_nodes = rev_back_bone_nodes
    
    old_start = aln_graph.consensus_start
    old_end = aln_graph.consensus_end
    if aln_graph.consensus_start != None:
        aln_graph.consensus_start = forward_keys[-1] - old_end
    if aln_graph.consensus_end != None:
        aln_graph.consensus_end = forward_keys[-1] - old_start
    # print aln_graph.consensus_start, aln_graph.consensus_end, forward_keys[-1]

def convert_mismatches(alnQ, alnT):
    qs = []
    ts = []
    assert len(alnQ) == len(alnT)
    alnQ = list(alnQ)
    alnT = list(alnT)
    for i in range(len(alnQ)):
        qb = alnQ[i]
        tb = alnT[i]
        if qb != tb and qb != "-" and tb != "-":
            qs.append("-")
            ts.append(tb)
            qs.append(qb)
            ts.append("-")
        else:
            qs.append(qb)
            ts.append(tb)
    alnQ = qs
    alnT = ts

    # push gaps to the right
    for i in range(len(alnQ) - 1):
        if alnT[i] == "-":
            j = i
            while 1:
                j += 1
                c = alnT[j]
                if c != "-" or j >= len(alnQ) - 1:
                    break
            if c == alnQ[i]:
                alnT[i] = c
                alnT[j] = "-"

        if alnQ[i] == "-":
            j = i
            while 1:
                j += 1
                c = alnQ[j]
                if c != "-" or j >= len(alnT) - 1:
                    break
            if c == alnT[i]:
                alnQ[i] = c
                alnQ[j] = "-"
    qs2 = []
    ts2 = []
    for i in range(len(alnQ)):
        if alnQ[i] != '-' or alnT[i] != '-':
            qs2.append(alnQ[i])
            ts2.append(alnT[i])

    return "".join(qs2), "".join(ts2)


class AlnEdge(object):
    def __init__(self, in_node, out_node):
        self.ID = id(self)
        self.in_node = in_node
        self.out_node = out_node
        in_node.addout_edge(self)
        out_node.add_in_edge(self)
        self.visited = False
        self.count = 0
        self.score = 0

    def increase_count(self):
        self.count += 1

    def set_score(self, s):
        self.score = s

    def add_to_score(self, s):
        self.score += s

    def __repr__(self):
        return "(edge_ID:%d, in_node:%s, out_node:%s)" % (self.ID, self.in_node.__repr__(), self.out_node.__repr__())


class AlnNode(object):
    def __init__(self, base):
        self.ID = id(self) # inbuilt id() function; 
        self.base = base
        self._in_edges = []
        self._out_edges = []
        self.is_backbone = False
        self.coverage = 0
        self.weight = 0
        self.backbone_node = None
        self.best_in_edge = None
        self.best_out_edge = None
        self.info = []

    def add_in_edge(self, in_edge):
        self._in_edges.append(in_edge)

    def addout_edge(self, out_edge):
        self._out_edges.append(out_edge)

    def increase_weight(self, w=1):
        self.weight += w

    def in_edges(self):
        return self._in_edges

    def __repr__(self):
        return "(node_id:%d, base:%s, b:%s, w:%d, c:%d)" % (
            self.ID, self.base, self.is_backbone, self.weight, self.coverage)


class AlnGraph(object):
    def __init__(self, backbone_seq):
        """
        AlnGraph is instantiated by giving a "backbone" sequence
        >>> g = AlnGraph("ACGTCTAT")
        """

        self.nodes = {}
        self.edges = {}
        self.backbone_nodes = {}
        self.backbone_to_reads = {}
        self.nodes_to_edge = {}
        self._max_node_id = 0
        self._max_edge_id = 0
        self.consensus_path = None
        self.consensus_str = None
        self.consensus_start = None
        self.consensus_end = None
        self.path_count = 0

        self.begin_node = AlnNode("B")
        self.begin_node.backbone_node = self.begin_node
        self.add_node(self.begin_node)
        self.read_range = {}
        last_node = self.begin_node

        # set up a barebone linear graph for the seed sequence
        for pos in range(len(backbone_seq)):
            node = AlnNode(backbone_seq[pos])
            self.backbone_nodes[pos] = node
            node.is_backbone = True
            node.backbone_node = node
            node.increase_weight()
            self.add_node(node)

            edge = AlnEdge(last_node, node)
            edge.set_score(0)
            edge.increase_count()
            self.add_edge(edge)
            last_node = node

        self.end_node = AlnNode("E")
        self.end_node.backbone_node = self.end_node
        self.add_node(self.end_node)
        edge = AlnEdge(last_node, self.end_node)
        edge.set_score(0)
        edge.increase_count()
        self.add_edge(edge)

        self.backbone_node_to_pos = dict(list(zip(list(self.backbone_nodes.values()),
                                             list(self.backbone_nodes.keys()))))

    def add_alignment(self, aln, rId=None):

        read_pos, read_end_pos, read_aln_seq = aln[0]
        backbone_pos, backbone_end_pos, backbone_aln_seq = aln[1]

        if rId != None:
            self.read_range[rId] = (backbone_pos, backbone_end_pos)

        last_node = self.begin_node

        for aln_pos in range(len(read_aln_seq)):
            read_base = read_aln_seq[aln_pos]
            backbone_base = backbone_aln_seq[aln_pos]

            if read_base == backbone_base:
                node = self.backbone_nodes[backbone_pos]
                node.increase_weight()
                node.backbone_node.coverage += 1

                if (last_node.ID, node.ID) not in self.nodes_to_edge:
                    edge = AlnEdge(last_node, node)
                    self.add_edge(edge)
                    edge.increase_count()
                else:
                    self.nodes_to_edge[(last_node.ID, node.ID)].increase_count()

                read_pos += 1
                backbone_pos += 1
                last_node = node
                if rId != None:
                    node.info.append(rId)

            elif read_base == "-" and backbone_base != "-":
                node = self.backbone_nodes[backbone_pos]
                node.backbone_node.coverage += 1
                backbone_pos += 1

            elif read_base != "-" and backbone_base == "-":
                node = AlnNode(read_base)
                node.increase_weight()
                node.backbone_node = self.backbone_nodes[backbone_pos]

                self.add_node(node)

                edge = AlnEdge(last_node, node)
                self.add_edge(edge)
                edge.increase_count()

                read_pos += 1
                last_node = node
                if rId != None:
                    node.info.append(rId)

            self.backbone_to_reads.setdefault(backbone_pos, set())

            if rId != None:
                self.backbone_to_reads[backbone_pos].add(rId)

        if (last_node.ID, self.end_node.ID) not in self.nodes_to_edge:
            edge = AlnEdge(last_node, self.end_node)
            edge.set_score(0)
            edge.increase_count()
            self.add_edge(edge)

    def add_edge(self, edge):
        edge.ID = id(edge)
        self.edges[edge.ID] = edge
        self.nodes_to_edge[(edge.in_node.ID, edge.out_node.ID)] = edge

    def add_node(self, node):
        self.nodes[node.ID] = node

    def delete_edge(self, edge):
        n1 = edge.in_node
        n2 = edge.out_node
        n1._out_edges.remove(edge)
        n2._in_edges.remove(edge)
        del self.edges[edge.ID]
        del edge

    def delete_node(self, node):
        for in_edge in node._in_edges:
            self.delete_edge(in_edge)
        for out_edge in node._out_edges:
            self.delete_edge(out_edge)
        del self.nodes[node.ID]
        del node

    def merge_in_nodes(self, node):

        node_groups = {}
        for in_edge in node._in_edges:
            in_node = in_edge.in_node
            if len(in_node._out_edges) == 1:
                node_groups.setdefault(in_node.base, [])
                node_groups[in_node.base].append(in_node)

        for b in node_groups:
            if len(node_groups[b]) <= 1:
                continue

            nodes = node_groups[b]

            for n in nodes[1:]:
                nodes[0]._out_edges[0].count += n._out_edges[0].count
                nodes[0].increase_weight(n.weight)
                nodes[0].info.extend(n.info)

            for n in nodes[1:]:
                for in_edge in n._in_edges:
                    n1 = in_edge.in_node
                    if (n1.ID, nodes[0].ID) not in self.nodes_to_edge:
                        e = AlnEdge(n1, nodes[0])
                        self.add_edge(e)
                        e.count = self.nodes_to_edge[(n1.ID, n.ID)].count
                        e.visited = self.nodes_to_edge[(n1.ID, n.ID)].visited
                    else:
                        e = self.nodes_to_edge[(n1.ID, nodes[0].ID)]
                        e.count += self.nodes_to_edge[(n1.ID, n.ID)].count
                        e.add_to_score(self.nodes_to_edge[(n1.ID, n.ID)].score)
                self.delete_node(n)

            self.merge_in_nodes(nodes[0])

    def merge_out_nodes(self, node):

        node_groups = {}

        for out_edge in node._out_edges:
            out_node = out_edge.out_node
            if len(out_node._in_edges) == 1:
                node_groups.setdefault(out_node.base, [])
                node_groups[out_node.base].append(out_node)

        for b in node_groups:
            if len(node_groups[b]) <= 1:
                continue

            nodes = node_groups[b]

            for n in nodes[1:]:
                nodes[0]._in_edges[0].count += n._in_edges[0].count
                nodes[0].increase_weight(n.weight)
                nodes[0].info.extend(n.info)

            for n in nodes[1:]:
                for out_edge in n._out_edges:
                    n2 = out_edge.out_node
                    if (nodes[0].ID, n2.ID) not in self.nodes_to_edge:
                        e = AlnEdge(nodes[0], n2)
                        self.add_edge(e)
                        e.count = self.nodes_to_edge[(n.ID, n2.ID)].count
                        e.visited = self.nodes_to_edge[(n.ID, n2.ID)].visited
                    else:
                        e = self.nodes_to_edge[(nodes[0].ID, n2.ID)]
                        e.count += self.nodes_to_edge[(n.ID, n2.ID)].count
                        e.add_to_score(self.nodes_to_edge[(n.ID, n2.ID)].score)

                self.delete_node(n)

    def merge_nodes(self):
        for e in list(self.edges.values()):
            e.visited = False

        seed_nodes = []

        for node_id, node in list(self.nodes.items()):
            if len(node._in_edges) == 0:
                seed_nodes.append(node)

        while 1:
            if len(seed_nodes) == 0:
                break
            node = seed_nodes.pop(0)
            self.merge_in_nodes(node)
            self.merge_out_nodes(node)

            for out_edge in node._out_edges:
                out_node = out_edge.out_node
                out_edge.visited = True
                in_edges = [e for e in out_node._in_edges if e.visited == False]
                if len(in_edges) == 0:
                    seed_nodes.append(out_node)

    def find_best_path(self):
        seed_nodes = []
        node_score = {}
        best_node_score_edge = {}

        for node_id, node in list(self.nodes.items()):
            if len(node._out_edges) == 0:
                seed_nodes.append(node)
                node_score[node.ID] = 0
                node.best_out_edge = None
                node.best_in_edge = None

        for edge in list(self.edges.values()):
            edge.visited = False

        while 1:

            if len(seed_nodes) == 0:
                break

            node = seed_nodes.pop(0)
            best_edge = None
            best_score = None
            for out_edge in node._out_edges:
                out_node = out_edge.out_node

                score = node_score[out_node.ID]

                backbone_node_coverage = out_node.backbone_node.coverage
                if out_node.is_backbone == True and out_node.weight == 1:
                    new_score = score - 10
                else:
                    new_score = out_edge.count - backbone_node_coverage * 0.5 + score
                if best_score == None or new_score > best_score:
                    best_edge = out_edge
                    best_score = new_score

            if best_edge != None:
                best_node_score_edge[node.ID] = best_edge
                node_score[node.ID] = best_score

            for in_edge in node._in_edges:
                in_node = in_edge.in_node
                in_edge.visited = True
                out_edges = [e for e in in_node._out_edges if e.visited == False]
                if len(out_edges) == 0:
                    seed_nodes.append(in_node)

        node = self.begin_node
        consensus_path = []
        
        con_start = len(self.backbone_nodes)
        con_end = 0

        while 1:
            consensus_path.append(node)
            if node.ID not in best_node_score_edge:
                break
            else:
                best_out_edge = best_node_score_edge[node.ID]
                node.best_out_edge = best_out_edge
                best_node_score_edge[node.ID].out_node.best_in_edge = best_node_score_edge[node.ID]
                node = best_node_score_edge[node.ID].out_node
                node.best_in_edge = best_out_edge

        # print consensus_path
        con_start = len(self.backbone_nodes)
        con_end = 0
        
        for node in self.backbone_nodes:
            # print self.backbone_nodes[node]
            if self.backbone_nodes[node] in consensus_path:
                # print "hhaa"
                if int(node) < con_start:
                    con_start = node
                if int(node) > con_end: 
                    con_end = node

        self.consensus_start = con_start
        self.consensus_end = con_end
        # print self.consensus_start
        # print self.consensus_end
        print("path score of consensus_path:", self.score(consensus_path))
        return consensus_path

    def generate_consensus(self, min_cov=0):

        self.merge_nodes()

        if self.consensus_path == None:
            self.consensus_path = self.find_best_path()

        s = []
        c = []
        cov = []    # a list for the coverage of all node in the path
        cov_good = []
        for n in self.consensus_path:
            if n not in [self.begin_node, self.end_node]:
                s.append(n.base)
                cov.append(n.backbone_node.coverage)
                if len(n.info) >= min_cov:
                    # print len(n.info), n.info
                    cov_good.append("1")
                else:
                    cov_good.append("0")

                rn = n.weight
                if n.is_backbone == True:
                    rn -= 1
                if n.best_out_edge != None:
                    en = n.best_out_edge.count
                else:
                    en = 0
                if n.best_in_edge != None:
                    en2 = n.best_in_edge.count
                else:
                    en2 = 0
                c.append((rn, en2, en))

        s = "".join(s)
        cov_good = "".join(cov_good)

        p = re.compile("1+")
        best_range = (0, 0)
        best_l = 0
        for m in p.finditer(cov_good):
            b, e = m.start(), m.end()
            if e - b > best_l:
                best_range = (b, e)
                best_l = e - b
        b, e = best_range
        # print "best range", b, e
        if e - b != 0:
            self.consensus_str = s[b:e]
            c = c[b:e]
            cov = cov[b:e]
        else:
            self.consensus_str = ""
            c = []

        return self.consensus_str, c, cov, self.consensus_start, self.consensus_end

    def get_high_entropy_nodes(self,
                               ignore_backbone=False,
                               coverage_th=20,
                               overlap_th=10,
                               entropy_th=0.5):
        from math import log

        if self.consensus_path == None:
            self.generate_consensus()

        backbone_node_to_pos = self.backbone_node_to_pos

        node_entropy = []

        for node_id in self.nodes:

            node = self.nodes[node_id]

            if node.backbone_node.coverage < coverage_th: continue
            if ignore_backbone and node.is_backbone: continue

            p = 1.0 * (len(node.info) + 1) / (node.backbone_node.coverage + 1)
            if abs(p - 1) < 1e-5 or p > 1:
                ent = 0
            else:
                ent = - p * log(p) - (1 - p) * log(1 - p)

            node_entropy.append([node_id, node, ent])

        node_entropy.sort(key=lambda x: -x[2])

        high_entropy_nodes = [n for n in node_entropy if n[2] > entropy_th]

        return node_entropy, high_entropy_nodes

    def get_sorted_nodes(self):
        seed_nodes = []
        # first loop should only add the "Start" node into seed_nodes
        for node_id, node in list(self.nodes.items()):
            if len(node._in_edges) == 0:
                seed_nodes.append(node)

        # all flagged as unvisited
        for edge in list(self.edges.values()):
            edge.visited = False

        sorted_nodes = []

        while 1:
            if len(seed_nodes) == 0:
                break

            # pop(index) removes the item at the given index from the list and returns the removed item.
            node = seed_nodes.pop(0)

            for in_edge in node._in_edges:
                in_node = in_edge.in_node
            # why this step? in_node not used
            sorted_nodes.append(node)

            for out_edge in node._out_edges:
                out_node = out_edge.out_node
                out_edge.visited = True
                in_edges = [e for e in out_node._in_edges if e.visited == False]
                if len(in_edges) == 0:
                    seed_nodes.append(out_node)

        return sorted_nodes

    def get_edge_count(self, src, des):
        find = False
        for edge in src._out_edges:
            if edge.out_node == des:
                flag = True
                return edge.count
        if not find:
            return 0
    
    def score(self, res):
        '''
        sum_score = 0
        for i in range(len(res)-2):
            sum_score += self.get_edge_count(res[i], res[i+1])
        return sum_score
        '''
        return sum([self.get_edge_count(res[i],res[i+1]) for i in range(len(res)-2)])

    def all_topo_sort_util(self, res, nodes_visited, src_node, des_node, out):
        '''
        node_dic = self.nodes
        for node_id, node in node_dic.items():
        '''
        # mark current src node as visited and store in path
        nodes_visited[src_node.ID] == True
        res.append(src_node)
        if src_node == des_node:
            if self.path_count >= 1000:
                return
            else:
                #print(">score_", self.score(res))
                base_path = [node.base for node in res]
                self.path_count += 1
                print(''.join(base_path[1:-1]), self.score(res), file=out)
            
        else:
            out_degree = sum([out_e.count for out_e in list(src_node._out_edges)])
            out_nodes_l = [out_e.out_node for out_e in list(src_node._out_edges)] 
            selected_nodes_l =[]           
            
            # only traverse to path that has weight higher than 1/(number of out edges):
            '''for out_e in src_node._out_edges:
                if float(out_e.count/out_degree) >= 1/len(out_nodes_l):
                    selected_nodes_l.append(out_e.out_node)
            #for next_node in out_nodes_l:
            for next_node in selected_nodes_l:
                if nodes_visited[next_node.ID] == False:
                    self.all_topo_sort_util(res, nodes_visited, next_node, des_node, out)    
            '''            
            # !! only traverse to path that has weight higher than 1/2 will fail to when node is diverged

            # sort the out_edge first and traverse from the most possible node
            out_edges_sorted = sorted(src_node._out_edges, key=operator.attrgetter('count'))
            # selected_nodes_l.append(out_edges_sorted[0].out_node)
            for out_e in out_edges_sorted:
                next_node = out_e.out_node
                if nodes_visited[next_node.ID] == False:
                    self.all_topo_sort_util(res, nodes_visited, next_node, des_node, out)
            
        res.pop()
        nodes_visited[src_node.ID] == False


    # call the all_topo_sort_util() recursively
    def all_topo_sort(self,out):
        '''
        edges_cp = self.edges.values()
        for edge in list(edges_cp):
            edge.visited = False
        '''
        s_t = time.time()
        nodes_visited = {node_id: False for node_id in list(self.nodes.keys())}
        res = []
        with open(out+"_path.fa", 'w+') as path_out:
            self.all_topo_sort_util(res, nodes_visited, self.begin_node, self.end_node, path_out)
        print("there are %d path found in the graph"%(self.path_count))
        print("Taken time: %s seconds"%(time.time() - s_t))

    def random_walk_util(self, res, nodes_visited, src_node, des_node, out):
        '''
        node_dic = self.nodes
        for node_id, node in node_dic.items():
        '''
        # mark current src node as visited and store in path
        nodes_visited[src_node.ID] == True
        res.append(src_node)
        if self.path_count >= 2000:
            return
        if src_node == des_node:
            if self.path_count >= 2000:
                return
            else:
                #print(">score_", self.score(res))
                base_path = [node.base for node in res]
                self.path_count += 1
                print(''.join(base_path[1:-1]), self.score(res), file=out)
            
        else:
            out_nodes_l = [out_e.out_node for out_e in list(src_node._out_edges)] 
            shuffle(out_nodes_l)       
            for next_node in out_nodes_l:
                if nodes_visited[next_node.ID] == False:
                    self.all_topo_sort_util(res, nodes_visited, next_node, des_node, out)            
        res.pop()
        nodes_visited[src_node.ID] == False

    def random_walk_k(self):
        k = 1000
        with open("random_walk_k.fa", 'w') as path_out:
            for i in range(k):
                rd_path = []
                rd_path.append(self.begin_node)
                while(rd_path[-1].base != 'E'):
                    out_nodes_l = [out_e.out_node for out_e in list(rd_path[-1]._out_edges)] 
                    shuffle(out_nodes_l) 
                    rd_path.append(out_nodes_l[0])
                base_path = [node.base for node in rd_path]
                print(''.join(base_path[1:-1]), self.score(rd_path), file=path_out)
            print("finish %d random walking"%i)
        return 0

    # call the all_topo_sort_util() recursively
    def random_walk(self,out):
        '''
        edges_cp = self.edges.values()
        for edge in list(edges_cp):
            edge.visited = False
        '''
        out_N = 1000
        s_t = time.time()
        nodes_visited = {node_id: False for node_id in list(self.nodes.keys())}
        res = []
        with open(out+"_random_path.fa", 'w+') as path_out:
            self.random_walk_util(res, nodes_visited, self.begin_node, self.end_node, path_out)
        print("there are %d path found in the graph"%(self.path_count))
        print("Taken time: %s seconds"%(time.time() - s_t))

    def get_read_node_vector(self, entropy_th=0.5):

        ne, hne = self.get_high_entropy_nodes(coverage_th=0, entropy_th=entropy_th)
        node_to_entropy = dict([(v[1], v[2]) for v in ne])
        read_ids = set()

        for n in list(self.nodes.values()):
            for r in n.info:
                read_ids.add(r)

        backbone_node_to_pos = self.backbone_node_to_pos

        sn = self.get_sorted_nodes()
        high_entropy_nodes = [(n,
                               backbone_node_to_pos[n.backbone_node],
                               node_to_entropy[n]) for n in sn if node_to_entropy[n] > entropy_th]

        read_to_nodes = {}
        for r_id in read_ids:
            read_to_nodes[r_id] = ["-"] * len(high_entropy_nodes)
            for i, n in enumerate(high_entropy_nodes):
                if r_id in n[0].info:
                    read_to_nodes[r_id][i] = n[0].base

        return read_to_nodes, high_entropy_nodes

    def output_consensus_fasta(self, f_name, rID):
        f = open(f_name, "a")
        print(">%s" % rID, file=f)
        for i in range(0, len(self.consensus_str), 60):
            print(self.consensus_str[i:i + 60], file=f)
        f.close()

    def filter_graph(self, threshold = 8):
        edge_to_remove = []
        for edge in self.edges:
            in_node = self.edges[edge].in_node
            out_node = self.edges[edge].out_node

            if in_node.backbone_node.coverage < out_node.backbone_node.coverage:
                coverage = in_node.backbone_node.coverage
            else:
                coverage = out_node.backbone_node.coverage
    
            if self.edges[edge].count <= 1 and coverage > threshold:
                edge_to_remove.append(self.edges[edge])

        for edge in edge_to_remove:
            self.delete_edge(edge)
        
        node_to_remove = []
        #  and (self.nodes[node].is_backbone == False)
        for node in self.nodes:
            if (not(self.nodes[node]._in_edges) and not(self.nodes[node]._out_edges)):
                node_to_remove.append(self.nodes[node])
            elif not(self.nodes[node]._in_edges) and (self.nodes[node] != self.begin_node):
                edge_add = AlnEdge(self.begin_node, self.nodes[node])
                edge_add.count = 1
                edge_add.score = 0
                self.add_edge(edge_add)
            elif not(self.nodes[node]._out_edges) and (self.nodes[node] != self.end_node):
                edge_add = AlnEdge(self.nodes[node], self.end_node)
                edge_add.count = 1
                edge_add.score = 0
                self.add_edge(edge_add)  

        for node in node_to_remove:
            self.delete_node(node)
        

    def jsOutput(self):
        """returns a list of strings containing a a description of the graph for viz.js, http://visjs.org"""
        # get the backbone of the graph
        backbone_path = self.backbone_nodes
        pathdict = {}
        for i, nodeID in enumerate(list(backbone_path.values())):
            pathdict[nodeID] = i*150

        lines = ['var nodes = [']

        ni = self.nodes.values()
        count = 0
        for node in ni:
            line = '    {id:'+str(node.ID)+', label: "'+node.base+'"'
            if node.ID in pathdict and count % 5 == 0:
                line += ', allowedToMoveX: false, x: ' + str(pathdict[node.ID]) + ', y: 0 , allowedToMoveY: true },'
            else:
                line += '},'
            lines.append(line)

        lines[-1] = lines[-1][:-1]
        lines.append('];')

        lines.append(' ')

        lines.append('var edges = [')
        ni =  self.nodes.values()
        for node in ni:
            nodeID = str(node.ID)
            for edge in node._out_edges:
                #target = str(edge)
                #target =  "(%d) -> (%d) " % (edge.in_node.ID, edge.out_node.ID)
                target =  "%d" % (edge.out_node.ID)
                #weight = str(len(node._outEdges[edge].labels)+1)
                weight = str(edge.count)
                lines.append('    {from: '+nodeID+', to: '+target+", arrows:'to', value: "+weight+'},')
            '''for alignededge in node.alignedTo:
                # These edges indicate alignment to different bases, and are
                # undirected; thus make sure we only plot them once:
                if node.ID > alignededge:
                    continue
                target = str(alignededge)
                lines.append('    {from: '+nodeID+', to: '+target+', value: 1, style: "dash-line"},')'''
        lines[-1] = lines[-1][:-1]
        lines.append('];')
        return lines
    
    def htmlOutput(self, outfile):
        header = """
                  <!doctype html>
                  <html>
                  <head>
                    <title>POA Graph Alignment</title>
                    <style type="text/css">
                        html, body {
                        font: 10pt arial;
                        }
                        #mynetwork {
                        width: 1800px;
                        height: 600px;
                        border: 1px solid lightgray;
                        }
                    </style>
                    <script type="text/javascript" src="https://visjs.github.io/vis-network/standalone/umd/vis-network.min.js"></script>
                  </head>
                  <body>
                  <div id="mynetwork"></div>
                  <script type="text/javascript">
                    // create a network
                    var nodes = null;
                    var edges = null;
                    var network = null;
                    function draw() {
                  """
        outfile.write(textwrap.dedent(header[1:]))
        lines = self.jsOutput()
        for line in lines:
            outfile.write(line+'\n')
        footer = """
                      // Instantiate our network object.
                    var container = document.getElementById('mynetwork');
                    var data = {
                    nodes: nodes,
                    edges: edges
                    };
                    var options = {
                    nodes: {
                        shape: 'dot',
                        scaling: {
                        customScalingFunction: function (min,max,total,value) {
                            return value/total;
                        },
                        min:5,
                        max:150
                        }
                    }
                    };
                    network = new vis.Network(container, data, options);
                }
                </script>
                <body onload="draw()">
                </body>
                </html>
                """
        outfile.write(textwrap.dedent(footer))