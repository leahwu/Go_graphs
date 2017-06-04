from collections import Counter
from scipy.stats import wilcoxon
import SDGErased as sdg_e
import PowerLawDistribution as pld
import matplotlib.pyplot as plt
import SDGRepeated as sdg_r
import networkx as nx
import operator
from ValidDegree import *


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

        self.graph_bi_seq = ([val for val in self.graph.in_degree().values()], [val for val in self.graph.out_degree().values()])
        # betweenness_centrality as attribute
        self.betweeness_centrality = nx.betweenness_centrality(self.graph)
        # page-rank
        self.page_rank = nx.pagerank(self.graph)


    def degrees_plot(self):
        """
        method that plots the degree distribution of bi-sequence and the generated simple graph   
        """

        G_in_degrees = self.graph_bi_seq[0]
        G_in_values = sorted(set(G_in_degrees.values()))
        G_in_hist = [Counter(G_in_degrees.values())[x] for x in G_in_values]
        plt.plot(G_in_values, G_in_hist, 'ro-', markersize=5)  # in-degree

        G_out_degrees = self.graph_bi_seq[1]
        G_out_values = sorted(set(G_out_degrees.values()))
        G_out_hist = [Counter(G_out_degrees.values())[x] for x in G_out_values]
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

        graph_in = self.graph_bi_seq[0]
        graph_out = self.graph_bi_seq[1]

        return [wilcoxon(in_seq,graph_in), wilcoxon(out_seq, graph_out)]

    def pk_vs_bc_plot(self):
        pr = self.page_rank
        bc = self.betweeness_centrality

        pr_sort = sorted(pr.items(), key=operator.itemgetter(1), reverse=True)

        nodes = [str(node[0]) for node in pr_sort ]
        pr_scores = [node[1] for node in pr_sort]
        bc_scores = [bc[node[0]] for node in pr_sort]
        bc_scaled_scores = [elem / sum(bc_scores) for elem in bc_scores]
        plt.plot(pr_scores, 'ro', markersize=3)
        plt.plot(bc_scaled_scores, 'bx', markersize=3)

        plt.legend(['Page_Rank', 'Betweeness_Centrality'])
        plt.xlabel('Node')
        plt.ylabel('Ranking')
        plt.title('Comparison between Pagerank and Betweeness centrality')











