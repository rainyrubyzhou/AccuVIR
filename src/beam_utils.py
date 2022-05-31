#################################################################################$$
# Author: Runzhou
# Some commonly used functions in Beam search
#################################################################################$$
import time
import numpy as np
def gen_path_in_tree(tree, leaf_l, k = 1):
    # input: root node of a tree; list of leaf nodes;
    # return: generator of all paths from root to leaf(follow the 'child->parent' pointer)
    time_1 = time.time()
    for nd in leaf_l:
        path = [nd.node_id()]
        while(nd._parent != None):
            nd = nd._parent
            path.append(nd._node_id)
        yield path
        # no need to return when using yield
    time_2 = time.time()
    #print("Constructing paths and generate sequences: %s s"%(time_2 - time_1))

def trunc_best_k(trunc_array, k):
    # input: a numpy array;  k to be selected
    # return: indices of k max element (indices in ascending order)
    indices = np.argpartition(trunc_array, -k)[-k :]
    return indices[np.argsort(trunc_array[indices])]