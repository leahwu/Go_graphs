from collections import Counter
import SDGErased as sdg_e
import PowerLawDistribution as pld
import matplotlib.pyplot as plt
import SDGRepeated as sdg_r
import networkx as nx
import operator
from ValidDegree import *
import scipy.stats as st


class CDMGenerator(object):

    def __init__(self, a, alpha, beta, n, algorithm, b=2):
        fg = pld.PowerLaw(a, alpha, beta, b)
        self.bi_seq = directed_gen(alpha, beta, fg, n)

        if algorithm == 'Erased':
            self.graph = sdg_e.gen_simple_DCM(self.bi_seq)
        if algorithm == 'RevisedRepeated':
            self.graph = sdg_r.gen_simple_DCM2(self.bi_seq)
        if algorithm == 'Repeated':
            self.graph = sdg_r.gen_simple_DCM(self.bi_seq)

        # return to a dictionary
        self.page_rank = nx.pagerank(self.graph)
        self.betweenness_centrality = nx.betweenness_centrality(self.graph)

        self.din = list(self.graph.in_degree().values())
        self.dout = list(self.graph.out_degree().values())



    def degrees_plot(self):
        """
        method that plots the degree distribution of bi-sequence and the generated simple graph   
        """
        plt.figure(1)
        G_in_values = sorted(set(self.din))
        G_in_hist = [Counter(self.din)[x] for x in G_in_values]
        plt.plot(G_in_values, G_in_hist, 'ro-', markersize=5)  # in-degree

        G_out_values = sorted(set(self.dout))
        G_out_hist = [Counter(self.dout)[x] for x in G_out_values]
        plt.plot(G_out_values, G_out_hist, 'bv-', markersize=5)  # out-degree

        # plot the sequence degree distribution
        in_degrees = self.bi_seq[0]
        in_values = sorted(set(in_degrees))
        in_hist = [Counter(in_degrees)[x] for x in in_values]
        plt.plot(in_hist, 'yo-', markersize=3)  # in-degree

        out_degrees = self.bi_seq[1]
        out_values = sorted(set(out_degrees))
        out_hist = [Counter(out_degrees)[x] for x in out_values]
        plt.plot(out_hist, 'cv-', markersize=3)  # out-degree

        plt.legend(['Graph In-degree', 'Graph Out-degree', 'Sequence In-degree', 'Sequence Out-degree'])
        plt.xlabel('Degree')
        plt.ylabel('Number of nodes')
        plt.xlim([0,40])


    def wilx_test(self):
        in_seq = self.bi_seq[0]
        out_seq = self.bi_seq[1]

        return [st.wilcoxon(in_seq,self.din), st.wilcoxon(out_seq, self.dout)]

    # plot the nodes in the decreasing order of pagerank, and then plot their corresponding
    # betweenness_centrality value
    def pr_vs_bc_plot(self):
        plt.figure(2)

        pr = self.page_rank
        bc = self.betweenness_centrality

        pr_sort = sorted(pr.items(), key=operator.itemgetter(1), reverse=True)

        # nodes = [str(node[0]) for node in pr_sort ]
        pr_scores = [node[1] for node in pr_sort]
        bc_scores = [bc[node[0]] for node in pr_sort]
        bc_scaled_scores = [elem / sum(bc_scores) for elem in bc_scores]
        plt.plot(pr_scores, 'ro', markersize=3)
        plt.plot(bc_scaled_scores, 'bx', markersize=3)

        plt.legend(['Page_Rank', 'Betweeness_Centrality'])
        plt.xlabel('Node')
        plt.ylabel('Ranking')
        plt.title('Comparison between Pagerank and Betweenness centrality')


    def spearman_test(self):
        pr = list(self.page_rank.values())
        bc = list(self.betweenness_centrality.values())
        corr, pvalue =st.spearmanr(pr, bc)

        return corr, pvalue

