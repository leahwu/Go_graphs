import numpy as np
import PowerLawDistribution as pld
from scipy.stats import poisson
import networkx as nx
import scipy.stats as st
import matplotlib.pyplot as plt


def gene_one_pair(c, beta):
    """
    method to generate one sample (d+, d-)
    :return: tuple with a pair of d+ and d-
    """
    w_sample = pld.generate_w(c, beta)
    # derive the sample degree from W+ and W-
    d = int(poisson(mu=w_sample).rvs())
    return d

def rvs(n, c, beta):
    d_seq = np.zeros(n, dtype=np.int)
    for i in range(0, n):
        d_seq[i] = gene_one_pair(c,beta)

    # resample to get even sum
    while sum(d_seq)%2 == 1:
        d_seq[0] = gene_one_pair(c,beta)

    return d_seq

def Graph_generate(seq):
    # generate the multigraph
    cm = nx.configuration_model(seq)

    # remove parallel edges
    cm = nx.DiGraph(cm)

    # remove self-loops
    cm.remove_edges_from(cm.selfloop_edges())
    return cm

# exogeneous parameters
c = 4.33
beta = 30.74
n = 1000

def ranking_spear_corr(graph):
    pr = list(nx.pagerank(graph).values())
    bc = list(nx.betweenness_centrality(graph).values())
    corr, pvalue =st.spearmanr(pr, bc)

    return corr, pvalue


rankcorr = []
for i in range(100):
    graph_temp = Graph_generate( rvs(n,c,beta) )
    rankcorr += [ranking_spear_corr(graph_temp)[0]]


mean_rankcorr = np.mean(rankcorr)