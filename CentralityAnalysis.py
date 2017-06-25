import DCMGenerator as dcm_g
import numpy as np
import networkx as nx
import scipy.stats as st
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
"""
# generate a list of DC graphs
DCMlist = [dcm_g.DCMGenerator(1, 1, 3, 1000, 'Erased') for i in range(100)]

# test for the graph centrality
graph_centra = [g.graph_centrality for g in DCMlist]

st.describe(graph_centra)

plt1 = plt.figure(1)
plt1.hist(graph_centra)
plt1.title("Histogram for graph btweeness centrality")
patch = mpatches.Patch(label='Mean: 0.105 Variance: 0.0034')
plt1.legend(handles=patch)
"""

# define a remove impact function
## tool function
def max_dic(dic):
    v = list(dic.values())
    k = list(dic.keys())
    return k[v.index(max(v))]

def graph_centra(graph, type = "btw"):
    N = nx.number_of_nodes(graph)
    if type == "btw":
        bc_dic = nx.betweenness_centrality(graph)
        bc_list = list(bc_dic.values())
        bcmax = max(bc_list)
        return (N * bcmax - sum(bc_list)) / (N-1)


def imp_remove(dcm, n, rule = "pagerank"):
    # copy the graph for manipulation
    graph_copy = dcm.graph
    # the list of centrality
    central = [dcm.graph_centrality]
    if rule == "pagerank":
        rank_copy = dict(dcm.page_rank)
    elif rule == "btwcentrality":
        rank_copy = dict(dcm.betweenness_centrality)

    elim = []  # eliminated list
    for i in range(n):
        node_lab = max_dic(rank_copy)
        elim += [node_lab]
        rank_copy.pop(node_lab)

        # remove the node
        graph_copy.remove_nodes_from([node_lab])

        # add the graph centrality
        central += [graph_centra(graph_copy)]

    return central


central_1 = imp_remove(DCMlist[0],20)
central_10 = imp_remove(DCMlist[0], 20, rule = "btwcentrality")

# plot
plt2 = plt.figure(2)
plt2.title('Graph centrality after eliminating highest ranking node Correlation: 0.968')
plt.xlabel('Eliminate nodes number')
plt.ylabel('Graph betweenness centrality')
plt.plot(central_10, color = 'blue', marker = 'o', markersize = "3")
plt.plot(central_1, color = 'red', marker = 'v', markersize = "3")
plt.legend(['Betweenness centrality', 'Pagerank'])

# check the correlation
st.pearsonr(central_1, central_10)
