import DCMGenerator as dcm_g
import time
import numpy as np


def test1():
    # assign the parameters
    # [Special Case [Alpha = Beta]]
    a = 0.5  # the lower bound for W+

    alpha = 3  # the power for W+
    beta = 2.5  # the power for W+

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

    # check the overlapping in first k-ranked nodes
    overlap_percentage = []
    overlap_number = []
    for k in range(10, 2000, 100):
        num, per = dcm.overlaps(k)
        overlap_number.append(num)
        overlap_percentage.append(per)
        print(k, per)


def test2():
    f = open("results1", 'w+')

    a_range = np.arange(0.8, 1.3, 0.1) # 5 choices

    #beta_range = np.arange(3.2,3.5,0.2) # 2.8, 3.0, 3.2, 3.4
    #d_range = np.arange(1.1, 1.3, 0.2) # 0.6, 0.8, 1.0, 1.2

    beta_range = np.arange(2.8, 3.3, 0.2) # 3 choices
    d_range = np.arange(0.8, 1.3, 0.2) # 3 choices
    n = 2000
# 5*3*3=45
    models = []
    s = 'ALl the graph\'s size is' + repr(n) + '.\n'
    i = 1
    for a in a_range:
        for beta in beta_range:
            for d in d_range:
                print(i)
                s += "Model"+ repr(i) +":\n"
                model = dcm_g.DCMGenerator(a, d, beta, n, 'Erased')
                s += '\n'

                s += "The params are:\n"
                for para in model.fg.params.items():
                    s += para[0] + ' = ' + "%0.5f" % para[1] + '\n'
                s += '\n'

                s += "Expectation of W^minus is " + repr(model.fg.e_w_minus) +'\n'
                s += "Expectation of W^plus is " + repr(model.fg.e_w_plus) + '\n'
                s += '\n'

                s += "Mean of original in-degree sequence is " + repr(model.mean_original_in_seq) + '\n'
                s += "Mean of original out-degree sequence is " + repr(model.mean_original_out_seq) + '\n'
                s += '\n'

                s += "After being modified by Algorithm 2.1, \n"
                s += "Mean of equal-sum in-degree sequence is " + repr(model.mean_equal_sum_in_seq) +'\n'
                s += "Mean of equal-sum out-degree sequence is " + repr(model.mean_equal_sum_out_seq) + '\n'
                s += '\n'

                s += "After removing self-loops and parallel edges:\n"
                s += "Mean of in-degree sequnces is " + repr(model.mean_in_degree) + '\n'
                s += "Mean of out-degree sequnces is " + repr(model.mean_out_degree) + '\n'
                s += '\n'

                s += "The percentage of overlapping nodes in top k ranked nodes by bc, and pr are:\n"
                s += "k\t : \t percentage of overlapping nodes \n"
                for k in range(50, 2000, 50):
                    s += repr(k) + '\t : \t' + repr(model.overlaps(k)[1]) + '\n'
                s += '\n'

                s += "Spearsman's rank correlation test:\n"
                corr, pvalue = model.spearman_test()
                s += "correlation = " + repr(corr) + ", pvalue = " + repr(pvalue) + "\n"

                models.append(model)

                s += '\n'
                i += 1

    f.write(s)
    f.close()

