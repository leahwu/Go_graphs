## Case when alpha not equals to beta
from ValidDegree import *
import SimpleDirectedGraph as SDG
import matplotlib.pyplot as plt
from PowerLawDistribution import *
from ValidDegree import *

## Case of different power alpha not equal to beta
# assign the parameters
a = 1 # the floor for W+
alpha = 2 # the power for W+
beta = 2.5 # the power for W+, setting beta to 2.5
n = 1000 # simulation times


fg_beta = PowerLaw(1, alpha, beta)
bi_seq_beta = directed_gen(alpha, beta, fg_beta, n)
D_beta = SDG.gen_simple_DCM(bi_seq_beta)

# plot the graph degree distribution
fig =plt.figure(2,figsize=(6,4))
fig.subplots_adjust(bottom=0.025, left=0.025, top = 0.975, right=0.975)
sub4 = fig.add_subplot(2,2, (1,2))
SDG.plot_hist(D_beta, 'Power_law_degree_different_power')

# plot the sequence degree distribution
sub5 = fig.add_subplot(2,2,3)
plt.hist(bi_seq_beta[0], bins ="auto", color = 'r' )
plt.hist(bi_seq_beta[1], bins ="auto", color = 'b')
plt.legend(['In-degree', 'Out-degree'])

# plot the sequence degree bivariate distribution
sub6 = fig.add_subplot(2,2,4)
plt.hist2d(bi_seq_beta[0],bi_seq_beta[1])
cbar = plt.colorbar()
cbar.ax.set_ylabel('Counts')
plt.xlabel('In-degree')
plt.ylabel('Out-degree')

# test for linear dependent a = 1.01
fg_a1 = PowerLaw(1.01,alpha, beta)
bi_seq_a1 = directed_gen(alpha, beta, fg_a1, n)
D_a1 = SDG.gen_simple_DCM(bi_seq_a1)
plt.figure(4)
SDG.plot_hist(D_a1, 'ErasedAlg_Power_law_degree_distribution_a=1.01 ')

