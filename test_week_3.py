import DCMGenerator as dcm_g
import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt
# tail distribution of page rank and betweenness centrality

# the correlation between in-degree distribution and out-degree distribution

# assign the parameters
    # [Special Case [Alpha = Beta]]

def test1():
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

def test2():
    """
    We need beta > 2, beta > 2 * d, beta > d + 1 for finite variance
    :return: qualified models
    """

    f = open('week3_2.txt', 'w+')

    a_range = np.arange(0.5, 4, 0.5)  # 7 choices

    beta_range = np.arange(2.6, 5.0, 0.4) # 7 choices
    d_range = np.arange(0.6, 2.5, 0.4)  # 6 choices

    n = 2000

    models = {}
    s = 'ALl the graph\'s size is ' + repr(n) + '.\n'

    rank_corr = []
    degree_corr = []
    i = 1
    for a in a_range:
        for beta in beta_range:
            for d in d_range:
                if not (beta > 2 * d and beta > d + 1):
                    continue

                model = dcm_g.DCMGenerator(a, d, beta, n, 'Erased')
                if model.mean_in_degree < 1:
                    continue

                print(i)
                s += "Model" + repr(i) + ":\n"

                s += '\n'

                s += "The params are:\n"
                for para in model.fg.params.items():
                    s += para[0] + ' = ' + "%0.5f" % para[1] + '\n'
                s += '\n'

                s += "Expectation of W^minus is " + repr(model.fg.e_w_minus) + '\n'
                s += "Expectation of W^plus is " + repr(model.fg.e_w_plus) + '\n'
                s += '\n'

                s += "Mean of in-degree sequence is " + repr(model.mean_in_degree) + '\n'
                s += "Mean of out-degree sequence is " + repr(model.mean_out_degree) + '\n'
                s += '\n'

                s += "Spearsman's rank correlation test:\n"
                corr, pvalue = model.spearman_test()
                rank_corr.append(corr)
                s += "correlation = " + repr(corr) + ", pvalue = " + repr(pvalue) + "\n"
                s += '\n'

                s += "Correlation between in-degree sequence and out-degree sequence is: \n"
                corr, pvalue = model.corr_in_and_out()
                degree_corr.append(corr)
                s += "corr = " + repr(corr) + ", pvalue = " + repr(pvalue) + "\n"

                models[i] = model

                s += '\n'
                i += 1

    corr = st.pearsonr(rank_corr, degree_corr)
    s += "Correlation between rank_corr and degree_corr is:\n"
    s += repr(corr)

    f.write(s)
    f.close()

    return models, rank_corr, degree_corr


models, rank_corr, degree_corr = test2()