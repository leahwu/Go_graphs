import networkx as nx
from ValidDegree import *
from PowerLawDistribution import *
import SimpleDirectedGraph as SDG
import matplotlib.pyplot as plt


# assign the parameters
a = 1 # the floor for W+
b = 2 # the floor for W-
alpha = 2 # the power for W+
beta = 2 # the power for W+

n = 2000 # simulation times

# the Power Law Distribution
fg = PowerLaw(a, b, alpha, beta)

# generate the sample of power law distribution
bi_seq = directed_gen(alpha, beta, fg, n)

# generate simple directed configuration graph
D = SDG.gen_simple_DCM(bi_seq)

# plot the degree distribution
plt.figure(1)
SDG.plot_hist(D, 'ErasedAlg_Power_law_degree_distribution')

# plot the graph
plt.figure(2)
SDG.plot_graph(D, 'ErasedAlg_Graph')

