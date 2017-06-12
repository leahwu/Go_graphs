import DCMGenerator as dcm_g
# tail distribution of page rank and betweenness centrality

# the correlation between in-degree distribution and out-degree distribution

# assign the parameters
    # [Special Case [Alpha = Beta]]
a = 1.5  # the lower bound for W+

alpha = 3  # the power for W+
beta = 2.5  # the power for W+

d = beta / alpha

n = 2000 # graph size

# generate simple directed configuration model
model = dcm_g.DCMGenerator(a, d, beta, n, 'Erased')

model.bc_vs_pr_dist()

corr, pvalue = model.corr_in_and_out()

model.bc_vs_pr_plot(200)
model.bc_vs_pr_plot(2000)

model.pr_vs_bc_plot(200)
model.pr_vs_bc_plot(2000)