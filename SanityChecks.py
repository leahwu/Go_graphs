import DCMGenerator as dcm_g
import PowerLawDistribution as pld
import  validate as vd
import numpy as np
from math import sqrt
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
from statsmodels.distributions.empirical_distribution import ECDF
from math import log
from scipy import stats

def test1(a, alpha, beta, m=2000):
    """
    Test EW+ = EW- ~ average of W Samples
    EW+ = b * alpha / (alpha - 1) = EW- = c * beta / (beta - 1)
    :return: 
    """
    model = pld.PowerLaw(a, alpha, beta)
    d = model.d

    expectation = model.alpha * model.b / (model.alpha - 1)
    expectation_2 = beta * model.c / (beta - 1)

    e_w_minus_square = beta * model.c**2 / (beta - 2)
    e_w_plus_square = model.alpha * model.b ** 2 / (model.alpha - 2)

    w_plus, w_minus = model.gene_iid_pairs_of_w(m)

    ave_plus = w_plus.mean()
    ave_minus = w_minus.mean()

    ave_w_plus_square = np.dot(w_plus, w_plus) / m
    ave_w_minus_square = np.dot(w_minus, w_minus) / m

    print(expectation, ave_plus)
    print(expectation_2, ave_minus)
    print(e_w_plus_square, ave_w_plus_square)
    print(e_w_minus_square, ave_w_minus_square)

    return w_plus, w_minus


def test2(c, beta, m=2000):
    w = [pld.generate_w(c, beta) for i in range(0, m)]
    w = np.array(w)

    expectation = beta * c / (beta - 1)
    ave_w = w.mean()

    e_w_square = beta * c**2 / (beta - 2)
    ave_w_square = np.dot(w, w) / m

    print(expectation, ave_w)
    print(e_w_square, ave_w_square)

    return w


def test3(a, alpha, beta, b=10, n=2000):
    """
    test for degree correlation
    """

    model = pld.PowerLaw(a, alpha, beta,b)
    a = model.a
    b = model.b
    c = model.c
    d = model.d
    beta = model.beta
    alpha = model.alpha

    d_in, d_out = model.rvs(n)

    e_d_in = b * alpha / (alpha - 1)
    ave_d_in = d_in.mean()
    print("ED+ =", e_d_in, 'Sample mean =', ave_d_in)

    e_d_out = model.c * beta / (beta - 1)
    ave_d_out = d_out.mean()
    print("ED- =", e_d_out, 'Sample mean = ', ave_d_out)

    e_prod = a * beta * c**(d+1) / (beta - d - 1)
    ave_prod = np.dot(d_in, d_out) / n
    print("ED+D- = ", e_prod, "Sample mean =", ave_prod)

    e_w_plus_square = alpha * b ** 2 / (alpha - 2)
    e_w_minus_square = beta * c ** 2 / (beta - 2)

    var_d_in = e_d_in + e_w_plus_square - e_d_in**2
    var_d_out = e_d_out + e_w_minus_square - e_d_out**2

    sample_var_d_in = d_in.var()
    sample_var_d_out = d_out.var()

    print("var_d_in = ", var_d_in, 'sample = ', sample_var_d_in)
    print("var_d_outm = ", var_d_out, 'sample = ', sample_var_d_out)

    corr = (e_prod - e_d_in * e_d_out) / sqrt(var_d_in * var_d_out)
    sample_corr = pearsonr(d_in, d_out)
    print('correlation(d_in, d_out) = ', corr, "sample corr = ", sample_corr[0])


def test4(a, alpha, beta, b=None):
    """
    test for correlation of the graph
    :param a: 
    :param alpha: 
    :param beta: 
    :param b: 
    :return: 
    """

    model = dcm_g.DCMGenerator(a, alpha, beta, 2000, 'Erased', b=b)

    corr1 = pearsonr(model.d_in_original, model.d_out_original)[0]
    # corr2 = pearsonr(model.d_in, model.d_out)[0]
    corr3 = pearsonr(model.graph_din, model.graph_dout)[0]

    corr = get_corr(a, alpha, beta, b)

    rank_corr = model.spearman_test()[0]

    return corr, corr1, corr3, rank_corr, model


def get_b(a, alpha, beta):
    """for alpha != beta"""
    if alpha == beta:
        raise ValueError('Compute the value of b for alpha != beta only!')
    return (alpha / (alpha - 1) * (beta - 1) / beta * a ** (alpha / beta)) ** (beta / (alpha - beta))


def get_corr(a, alpha, beta, b=None):

    if b is None:
        b = get_b(a, alpha, beta)

    d = beta / alpha

    c = (b / a) ** (alpha / beta)

    e_d_in = b * alpha / (alpha - 1)
    e_d_out = c * beta / (beta - 1)

    e_w_plus_square = alpha * b ** 2 / (alpha - 2)
    e_w_minus_square = beta * c ** 2 / (beta - 2)

    e_prod = a * beta * c ** (d + 1) / (beta - d - 1)
    var_d_in = e_d_in + e_w_plus_square - e_d_in ** 2
    var_d_out = e_d_out + e_w_minus_square - e_d_out ** 2

    corr = (e_prod - e_d_in * e_d_out) / sqrt(var_d_in * var_d_out)

    return corr


def get_expected_degree(alpha, b):
    return b * alpha / (alpha - 1)


def test5():
    corr = []
    corr1 = []
    corr2 = []
    models = []
    rank_corr = []

    alpha = 5
    beta = 5
    a = 1
    # set up b_range

    m = 50  # test for page rank

    b_range = np.arange(6, 56, 5)
    for b in b_range:
        rank_corr_i = []
        c1_i = []
        c2_i = []
        model = None
        c = 0
        for i in range(0, m):
            c, c1, c2, r_c, model = test4(a, alpha, beta, b=b)
            rank_corr_i.append(r_c)
            c1_i.append(c1)
            c2_i.append(c2)

        c1_i = np.array(c1_i)
        c2_i = np.array(c2_i)
        rank_corr_i = np.array(rank_corr_i)

        corr.append(c)
        corr1.append(c1_i)
        corr2.append(c2_i)
        rank_corr.append(rank_corr_i)
        models.append(model)

    corr = np.array(corr)

    return b_range, corr, corr1, corr2, rank_corr, models


def test6(a, alpha, beta, b=None):
    if b is None:
        b = get_b(a, alpha, beta)
    print(b, get_corr(a, alpha, beta, b), get_expected_degree(alpha, b))


def test7():
    corr = []
    corr1 = []
    corr2 = []
    models = []
    rank_corr = []

    alpha = 4
    beta = 3
    a_range = np.arange(1.3, 3.3, 0.3)
    # set up b_range

    m = 10  # simulation times

    for a in a_range:
        rank_corr_i = []
        c1_i = []
        c2_i = []
        model = None
        c = 0
        for i in range(0, m):
            c, c1, c2, r_c, model = test4(a, alpha, beta)
            rank_corr_i.append(r_c)
            c1_i.append(c1)
            c2_i.append(c2)

        c1_i = np.array(c1_i)
        c2_i = np.array(c2_i)
        rank_corr_i = np.array(rank_corr_i)

        corr.append(c)
        corr1.append(c1_i)
        corr2.append(c2_i)
        rank_corr.append(rank_corr_i)
        models.append(model)

    print('finish loops')

    corr = np.array(corr)

    return a_range, corr, corr1, corr2, rank_corr, models


def get_mean(seqs):
    """
    seqs should be like: [[],[],[]]
    :param seqs: 
    :return: 
    """
    mean_seq = []
    for seq in seqs:
        mean_seq.append(seq.mean())

    return mean_seq


def get_std(seqs):
    std_seq = []
    for seq in seqs:
        std_seq.append(seq.std())

    return std_seq


def get_loglog(seq):
    """
    help get the tail distribution logx, logy of the seq
    :param seq: 
    :return: 
    """
    cdf = ECDF(seq)
    for x in cdf.x:
        if x <= 0:
            cdf.x = cdf.x[1:]
            cdf.y = cdf.y[1:]
    # eliminating the cdf = 1 term
    cdf.x = cdf.x[:-1]
    cdf.y = cdf.y[:-1]

    logx = [log(x) for x in cdf.x]
    logy = [log(1 - y) for y in cdf.y]

    return logx, logy


def abline(slope, intercept):
    """Plot a line from slope and intercept"""
    axes = plt.gca()
    x_vals = np.array(axes.get_xlim())
    y_vals = intercept + slope * x_vals
    plt.plot(x_vals, y_vals, '--')


def plot_tail(model):
    model.plot_tail_dist_loglog(list(model.page_rank.values()), 'page rank')
    model.plot_tail_dist_loglog(list(model.betweenness_centrality.values()), 'betweenness centrality')
    model.plot_tail_dist_loglog(model.graph_din, 'in-degree sequence')
    model.plot_tail_dist_loglog(model.graph_dout, 'out-degree sequence')
    plt.legend()


def linear_fit(x, y):
    # seq : list-like// graph.din, list(page_rank.values())
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    return slope, intercept


def d_tolist(d):
    return list(d.values())


