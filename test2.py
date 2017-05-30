import networkx as nx
import ValidDegree as VD
import PowerLawDistribution as PLD
import SDGRepeated as SDGr
import matplotlib.pyplot as plt



# assign the parameters
a = 1 # the floor for W+
alpha = 2 # the power for W+
beta = 2.5 # the power for W+
b = 2
n = 2000 # simulation times

# the Power Law Distribution
fg = PLD.PowerLaw(a,alpha, beta, b)

# generate the sample of power law distribution
bi_seq = VD.directed_gen(alpha, beta, fg, n)

# generate simple directed configuration graph
D = SDGr.gen_simple_DCM2(bi_seq)

# plot the bi-degree distribution with generated simple DCM degree distribution
fig =plt.figure(1,figsize=(6,4))
fig.subplots_adjust(bottom=0.025, left=0.025, top = 0.975, right=0.975)

# plot the graph degree distribution
sub1 = fig.add_subplot(2,2, (1,2))
SDGr.plot_hist(D, 'RevisedRepeatedAlg_Power_law_degree_distribution_2')

# plot the sequence degree distribution
sub2 = fig.add_subplot(2,2,3)
plt.hist(bi_seq[0], bins ="auto", color = 'r' )
plt.hist(bi_seq[1], bins ="auto", color = 'b')
plt.legend(['In-degree', 'Out-degree'])

# plot the sequence degree bivariate distribution
sub3 = fig.add_subplot(2,2,4)
plt.hist2d(bi_seq[0],bi_seq[1])
cbar = plt.colorbar()
cbar.ax.set_ylabel('Counts')
plt.xlabel('In-degree')
plt.ylabel('Out-degree')

# plot the graph
plt.figure(2)
SDGr.plot_graph(D, 'RevisedRepeatedAlg_Graph_2')


