import networkx as nx
import ValidDegree as VD
import PowerLawDistribution as PLD
import SDGRepeated as SDGr
import matplotlib.pyplot as plt


# assign the parameters
a = 1 # the floor for W+
alpha = 2 # the power for W+
beta = 2 # the power for W+
b = 2
n = 2000 # simulation times

# the Power Law Distribution
fg = PLD.PowerLaw(a,alpha, beta, b)

# generate the sample of power law distribution
bi_seq = VD.directed_gen(alpha, beta, fg, n)

# generate simple directed configuration graph
D = SDGr.gen_simple_DCM2(bi_seq)

# plot the degree distribution
plt.figure(1)
SDGr.plot_hist(D, 'test21')

# plot the graph
plt.figure(2)
SDGr.plot_graph(D, 'test22')

