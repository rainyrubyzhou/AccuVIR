# define the BeamSearchTree here; this tree is built during the Beam search process.
from networkx.algorithms.shortest_paths import weighted
from networkx.readwrite.json_graph.node_link import node_link_data


class BeamSearchTreeNode(object):
    
    def __init__(self, node_id = -1, local_weight = 0, sum_weight = 0, parent_nd = None) -> None:
        # para: edge; out_node as the new node; 
        self._node_id = node_id 
        self._local_weight = local_weight
        self._sum_weight = sum_weight  #store the accumulated weight
        self._children = []
        self._parent = parent_nd

    def node_id(self):
        return self._node_id
    
    def local_weight(self):
        return self._local_weight

    def sum_weight(self):
        return self._sum_weight

    def children(self):
        return self._children
    
    def parent(self):
        return self._parent
    
    def add_child(self, child):
        self._children.append(child)  
    
    def __str__(self, level=0):
        # print the tree line by line
        ret = "\t"*level+repr(self._node_id)+"\n"
        for child in self.children():
            ret += child.__str__(level+1)
        return ret
