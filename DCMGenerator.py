from collections import Counter
import PowerLawDistribution as pld
import matplotlib.pyplot as plt
import networkx as nx
import operator,math
import ValidDegree as vd
import scipy.stats as st
import DCMRevised as dcm_r
import numpy as np
from statsmodels.distributions.empirical_distribution import ECDF


class DCMGenerator(object):

    def __init__(self, a, d, beta, n, algorithm):


        if algorithm == 'Erased':
            self.fg = pld.PowerLaw(a, d, beta)
            degree_seq = vd.directed_gen(d, beta, self.fg, n)

            # after modifying the degree sequence to make the sum(d_in) = sum(d_out) using alg 2.1
            self.mean_equal_sum_in_seq = np.mean(degree_seq[0])
            self.mean_equal_sum_out_seq = np.mean(degree_seq[1])
            self.d_in = degree_seq[0].tolist()
            self.d_out = degree_seq[1].tolist()

            # original degree sequence generated by fg-distribution
            self.mean_original_in_seq = np.mean(degree_seq[2])
            self.mean_original_out_seq = np.mean(degree_seq[3])
            self.d_in_original = degree_seq[2].tolist()
            self.d_out_original = degree_seq[3].tolist()

            # generate the multigraph
            dcm = nx.directed_configuration_model(self.d_in, self.d_out)

            # remove parallel edges
            dcm = nx.DiGraph(dcm)
            # remove self-loops
            dcm.remove_edges_from(dcm.selfloop_edges())

            # get the simple directed configuration graph
            self.graph = dcm

        if algorithm == 'Repeated':
            flag = False
            while not flag:
                self.fg = pld.PowerLaw(a, d, beta)
                degree_seq = vd.directed_gen(d, beta, self.fg, n)

                # after modifying the degree sequence to make the sum(d_in) = sum(d_out) using alg 2.1
                self.d_in = degree_seq[0].tolist()
                self.d_out = degree_seq[1].tolist()

                # original degree sequence generated by fg-distribution
                self.d_in_original = degree_seq[2].tolist()
                self.d_out_original = degree_seq[3].tolist()

                (model, flag) = dcm_r.directed_configuration_model_revised(self.d_in, self.d_out)
            self.graph = nx.DiGraph(model)

        # return to a dictionary
        self.page_rank = nx.pagerank(self.graph)
        self.betweenness_centrality = nx.betweenness_centrality(self.graph)

        self.graph_din = list(self.graph.in_degree().values())
        self.graph_dout = list(self.graph.out_degree().values())
        self.size = len(self.graph_din)
        self.mean_in_degree = sum(self.graph_din) / self.size
        self.mean_out_degree = sum(self.graph_dout) / self.size

        BC = list(self.betweenness_centrality.values())
        bcmax = max(BC)
        N = len(BC)
        self.graph_centrality = (N * bcmax - sum(BC)) / (N-1)



    def __str__(self):
        s = "The params are:\n"
        s += "Expectation of W^minus is " + repr(self.fg.e_w_minus) + '\n'
        s += "Expectation of W^plus is " + repr(self.fg.e_w_plus) + '\n'
        s += '\n'

        s += "Mean of original in-degree sequence is " + repr(self.mean_original_in_seq) + '\n'
        s += "Mean of original out-degree sequence is " + repr(self.mean_original_out_seq) + '\n'
        s += '\n'

        s += "After being modified by Algorithm 2.1, \n"
        s += "Mean of equal-sum in-degree sequence is " + repr(self.mean_equal_sum_in_seq) + '\n'
        s += "Mean of equal-sum out-degree sequence is " + repr(self.mean_equal_sum_out_seq) + '\n'
        s += '\n'

        s += "After removing self-loops and parallel edges:\n"
        s += "Mean of in-degree sequence is " + repr(self.mean_in_degree) + '\n'
        s += "Mean of out-degree sequence is " + repr(self.mean_out_degree) + '\n'
        s += '\n'

        return s

    def plot_helper(self, seq, c, m, ms):
        values = sorted(set(seq))
        hist = [Counter(seq)[x] for x in values]
        plt.plot(values, hist, color=c, marker=m, markersize=ms)  # in-degree


    def test_equal_sum_algorithm(self):
        """
        method that plots the degree distribution of 
        1) bi-sequence, generated by fg-distribution
        2) bi-sequence, modified by algorithm 2.1 to make equal sum   
        """
        fig = plt.figure()
        self.plot_helper(self.d_in_original, 'blue', 'o', 5)
        self.plot_helper(self.d_out_original, 'green', 'v', 5)
        self.plot_helper(self.d_in, 'red', 'o', 5)
        self.plot_helper(self.d_out, 'cyan', 'v', 5)

        plt.legend(['Original In-degree Sequence', 'Original Sequence Out-degree Sequence',
                    'Equal-sum In-degree Sequence', 'Equal-sum Out-degree Sequence'])
        plt.xlabel('Degree')
        plt.ylabel('Number of nodes')
        plt.xlim([0, 40])

        txt = ''
        for para in self.fg.params.items():
                    txt += para[0] + ' = ' + "%0.2f" % para[1] + ' '

        plt.title(txt)

        plt.show()

    def degrees_plot(self):
        """
        method that plots the degree distribution of 
        1) bi-sequence, modified by algorithm 2.1 to make equal sum
        2) generated simple graph   
        """

        fig = plt.figure()
        self.plot_helper(self.d_in, 'blue','o', 5)
        self.plot_helper(self.d_out, 'green', 'v', 5)
        self.plot_helper(self.graph_din, 'red', 'o', 5)
        self.plot_helper(self.graph_dout, 'cyan','v', 5)

        plt.legend(['Equal-sum In-degree Sequence ', 'Equal-sum Out-degree Sequence', 'Graph In-degree Sequence ',
                    'Graph Out-degree Sequence'])
        plt.xlabel('Degree')
        plt.ylabel('Number of nodes')
        plt.xlim([0,40])

        txt = ''
        for para in self.fg.params.items():
            txt += para[0] + ' = ' + "%0.2f" % para[1] + ' '

        plt.title(txt)

        plt.show()


    def wilx_test(self):
        return [st.wilcoxon(self.d_in, self.d_in_original), st.wilcoxon(self.d_out, self.d_out_original)]

    # plot first k nodes in the decreasing order of pagerank, and then plot their corresponding
    # betweenness_centrality value
    def pr_vs_bc_plot(self, k=None):
        if k is None:
            k = self.size

        fig = plt.figure()

        pr = self.page_rank
        bc = self.betweenness_centrality

        pr_sort = sorted(pr.items(), key=operator.itemgetter(1), reverse=True)

        # nodes = [str(node[0]) for node in pr_sort ]
        pr_scores = [node[1] for node in pr_sort]
        bc_scores = [bc[node[0]] for node in pr_sort]
        bc_scaled_scores = [elem / sum(bc_scores) for elem in bc_scores]
        plt.plot(bc_scaled_scores[0: k], 'bx', markersize=3)
        plt.plot(pr_scores[0: k], 'ro', markersize=1)


        plt.legend(['Betweeness_Centrality', 'Page_Rank'])
        plt.xlabel('Node')
        plt.ylabel('Ranking')
        txt = ''
        for para in self.fg.params.items():
            txt += para[0] + ' = ' + "%0.2f" % para[1] + ' '

        plt.title(txt)
        plt.title('First '+repr(k) + ' nodes in decreasing order of Page Rank \n' + txt)


        # plot first k nodes in the decreasing order of bc, and then plot their corresponding pagerank
    def bc_vs_pr_plot(self, k=None):
        if k is None:
            k = self.size

        plt.figure()

        pr = self.page_rank
        bc = self.betweenness_centrality

        bc_sort = sorted(bc.items(), key=operator.itemgetter(1), reverse=True)

        # nodes = [str(node[0]) for node in pr_sort ]
        bc_scores = [node[1] for node in bc_sort]
        pr_scores = [pr[node[0]] for node in bc_sort]
        bc_scaled_scores = [elem / sum(bc_scores) for elem in bc_scores]
        plt.plot(pr_scores[0: k], 'ro', markersize=1)
        plt.plot(bc_scaled_scores[0: k], 'bx', markersize=3)

        plt.legend(['Page_Rank', 'Betweenness_Centrality'])
        plt.xlabel('Node')
        plt.ylabel('Ranking')
        txt = ''
        for para in self.fg.params.items():
            txt += para[0] + ' = ' + "%0.2f" % para[1] + ' '

        plt.title(txt)
        plt.title('First ' + repr(k) + ' nodes in decreasing order of Betweenness centrality\n' + txt)



    def spearman_test(self):
        pr = list(self.page_rank.values())
        bc = list(self.betweenness_centrality.values())
        corr, pvalue =st.spearmanr(pr, bc)

        return corr, pvalue




    # return to the number and the percentage of overlapping nodes in top k page-rakned
    # and top k betweenness-centrality nodes
    def overlaps(self, k):
        bc_sort = sorted(self.betweenness_centrality.items(), key=operator.itemgetter(1), reverse=True)
        pr_sort = sorted(self.page_rank.items(), key=operator.itemgetter(1), reverse=True)

        bc_topk = bc_sort[0: k]
        pr_topk = pr_sort[0: k]

        set_bc = {tup[0] for tup in bc_topk}
        set_pr = {tup[0] for tup in pr_topk}

        overlap_number = len(set.intersection(set_bc, set_pr))
        overlap_percentage = overlap_number / k

        return overlap_number, overlap_percentage

    def corr_in_and_out(self):
        """
        Calculate the correlation between in-degree sequence and out-degree sequence.
        While executing, plt the graph of in_seq and out_seq
        :return: the correlation
        """
        corr, p_value = st.pearsonr(self.d_in, self.d_out)

        return corr, p_value

    def plot_in_and_out(self):
        corr, p_value = st.pearsonr(self.d_in, self.d_out)

        plt.figure()

        self.plot_helper(self.d_in, 'red', 'o', 5)
        self.plot_helper(self.d_out, 'cyan', 'v', 5)

        plt.legend(['In-degree sequence', 'Out-degree sequence'])
        plt.xlabel('Degree')
        plt.ylabel('Number of nodes')
        plt.xlim([0, 40])

        txt = ''
        for para in self.fg.params.items():
            txt += para[0] + ' = ' + "%0.2f" % para[1] + ' '

        plt.title("Correlation: " + repr(corr) + " p-value: " + repr(p_value) + '\n' + txt)

    def plot_tail_dist(self, d, name):
        """
        Plot the tail distribution of data
        1-F(x), where F(x) is the empirical distribution of data
        :param d: self.page_rank or self.betweenness_centrality, type: dictionary
        :param name: name of d   
        :return: void
        """

        data = list(d.values())
        cdf = ECDF(data)

        # filter out the zero term
        for x in cdf.x:
            if x <= 0:
                cdf.x = cdf.x[1:]
                cdf.y = cdf.y[1:]
        # eliminating the cdf = 1 term
        cdf.x = cdf.x[:-1]
        cdf.y = cdf.y[:-1]

        plt.plot([math.log(elem) for elem in cdf.x], [math.log(1 - elem) for elem in cdf.y], label=name, marker='<', markerfacecolor='none', markersize=1)

    def bc_vs_pr_dist(self):
        """
        Plot the tail distribution of page rank and betweenness centrality
        :return: Void
        """
        plt.figure()
        self.plot_tail_dist(self.betweenness_centrality, 'betweenness centrality')
        self.plot_tail_dist(self.page_rank, 'page rank')
        plt.legend()

        txt = ''
        for para in self.fg.params.items():
            txt += para[0] + ' = ' + "%0.2f" % para[1] + ' '

        plt.title(txt)
        plt.title('Log Tail distribution\n' + txt)

        plt.show()


