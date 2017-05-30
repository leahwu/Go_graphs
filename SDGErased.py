import networkx as nx
import matplotlib.pyplot as plt
import pylab
from collections import Counter

# generate the simple directed configuration graph from the sample bi-degree sequence
def gen_simple_DCM(bi_degree):
    # convert the bi-seq array to list
    d_in = bi_degree[0].tolist()
    d_out = bi_degree[1].tolist()

    # generate the multigraph
    dcm = nx.directed_configuration_model(d_in, d_out)

    # remove parallel edges
    dcm = nx.DiGraph(dcm)
    # remove self-loops
    dcm.remove_edges_from(dcm.selfloop_edges())

    # get the simple directed configuration graph
    return dcm


def plot_graph(G, title):
    pos = nx.spring_layout(G)
    nx.draw(G, pos, arrows = True, node_size=0.5)
    plt.title(title)
    plt.savefig(title + '.eps')



def plot_hist(G, title):
    in_degrees = G.in_degree()  # dictionary node:degree
    in_values = sorted(set(in_degrees.values()))
    in_hist = [Counter(in_degrees.values())[x] for x in in_values]
    plt.plot(in_values, in_hist,'ro-')  # in-degree

    out_degrees = G.out_degree()
    out_values = sorted(set(out_degrees.values()))
    out_hist = [Counter(out_degrees.values())[x] for x in out_values]
    plt.plot(out_values, out_hist, 'bv-')  # out-degree

    plt.legend(['In-degree', 'Out-degree'])
    plt.xlabel('Degree')
    plt.ylabel('Number of nodes')
    plt.title(title)
    plt.savefig(title + '.eps')

