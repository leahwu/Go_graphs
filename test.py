import networkx as nx
from ValidDegree import *
from PowerLawDistribution import *
import SimpleDirectedGraph as sdg


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
D = sdg.gen_simple_DCM(bi_seq)
sdg.plot_hist(D)

# plot the graph
draw_exp = directed_gen(alpha, beta, fg, 100)
draw_graph = sdg.gen_simple_DCM(draw_exp)
# plot the graph
sdg.plot_graph(draw_graph)

