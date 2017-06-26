import DCMGenerator as dcm_g
import PowerLawDistribution as pld
import  validate as vd
import numpy as np
from math import sqrt
from scipy.stats import pearsonr

def test1(a, d, beta, m=2000):
    """
    Test EW+ = EW- ~ average of W Samples
    EW+ = b * alpha / (alpha - 1) = EW- = c * beta / (beta - 1)
    :return: 
    """
    model = pld.PowerLaw(a, d, beta)
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


def test3(a ,d, beta, n=2000):
    """
    test for degree correlation
    :param a: 
    :param d: 
    :param beta: 
    :param m: 
    :return: 
    """

    model = pld.PowerLaw(a, d, beta)
    a = model.a
    b = model.b
    c = model.c
    d = model.d
    beta =model.beta
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

