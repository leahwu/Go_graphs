import DCMGenerator as cdm_g
import SDGRepeated as SDGr
import matplotlib.pyplot as plt
import time

start = time.clock()

# assign the parameters
a = 1 # the floor for W+
alpha = 2 # the power for W+
beta = 2 # the power for W+
b = 2
n = 2000 # simulation times

# generate simple directed configuration model
dcm = cdm_g.CDMGenerator(a,alpha, beta, n, 'RevisedRepeated')
D = dcm.graph

# plot the bi-degree distribution with generated simple DCM degree distribution
fig =plt.figure(1,figsize=(6,4))
fig.subplots_adjust(bottom=0.025, left=0.025, top = 0.975, right=0.975)

# plot the graph degree distribution
sub1 = fig.add_subplot(2,2, (1,2))
SDGr.plot_hist(D, 'Case2RevisedRepeatedAlg_Power_law_degree_distribution')

# plot the sequence degree distribution
sub2 = fig.add_subplot(2,2,3)
plt.hist(dcm.bi_seq[0], bins ="auto", color = 'r' )
plt.hist(dcm.bi_seq[1], bins ="auto", color = 'b')
plt.legend(['In-degree', 'Out-degree'])

# plot the sequence degree bivariate distribution
sub3 = fig.add_subplot(2,2,4)
plt.hist2d(dcm.bi_seq[0],dcm.bi_seq[1])
cbar = plt.colorbar()
cbar.ax.set_ylabel('Counts')
plt.xlabel('In-degree')
plt.ylabel('Out-degree')

# plot the graph
plt.figure(2)
SDGr.plot_graph(D, 'Case2RevisedRepeatedAlg_Graph')


elapsed = (time.clock() - start)
print("Time used:",elapsed)