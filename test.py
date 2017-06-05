from collections import Counter
import SDGErased as SDG
import DCMGenerator as cdm_g
import matplotlib.pyplot as plt
import time

start = time.clock()

# assign the parameters
# [Special Case [Alpha = Beta]]
a = 2  # the lower bound for W+
alpha = 4  # the power for W+
beta = 3  # the power for W+
b = 2.5

n = 2000 # simulation times



# generate simple directed configuration model
dcm = cdm_g.CDMGenerator(a,alpha, beta, n, 'Erased')

# plot the bi-degree distribution to compare with generated simple DCM degree distribution
dcm.degrees_plot()

#
dcm.pr_vs_bc_plot()


# calculate the rank corr of the dcm
corr, pvalue = dcm.spearman_test()