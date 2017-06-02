from collections import Counter
import SDGErased as sdg_e
import PowerLawDistribution as pld
import matplotlib.pyplot as plt
import SDGRepeated as sdg_r
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

        self.graph_bi_seq = (self.graph.in_degree(), self.graph.out_degree())


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

        def stat_


