
import DCMGenerator as dcm_g
import time

start = time.clock()

# assign the parameters
# [Special Case [Alpha = Beta]]
a = 1.2  # the lower bound for W+

alpha = 3  # the power for W+
beta = 3.5  # the power for W+

d = beta / alpha

n = 2000 # simulation times



# generate simple directed configuration model
dcm = dcm_g.DCMGenerator(a, d, beta, n, 'Erased')

# plot the bi-degree sequence before and after adjusting to make equal sum
dcm.test_equal_sum_algorithm()

# plot the bi-degree distribution to compare with generated simple DCM degree distribution
dcm.degrees_plot()

#
dcm.pr_vs_bc_plot()
dcm.bc_vs_pr_plot()

dcm.pr_vs_bc_plot(200)
dcm.bc_vs_pr_plot(200)

# calculate the rank corr of the dcm
corr, pvalue = dcm.spearman_test()