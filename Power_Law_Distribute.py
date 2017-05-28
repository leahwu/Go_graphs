import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
from math import exp, factorial
from scipy.stats import rv_continuous, rv_discrete, poisson, describe, histogram


class power_gen(rv_continuous):
    "Power Tail Law"
    def _pdf(self, x, ceil, tao ):
        return (tao/ceil) * (x / ceil) ** (-tao-1)

def bi_degree(a, d, alpha, c, n):
    """
    a function that generates sample pair of bidegree 
    the function returns a tuple  (in-degree, out-degree)
    
    Input Parameters
    -----------------
    a : coefficient
    d : power term
    alpha : power for W+ CDF
    c : W- ceiling
    n: simulation sample number
    
    """
    b = a * c ** d   # W+ ceiling
    beta = alpha * d   # power for W- CDF

    # generate the power law distribution specifically for W-
    power_law_minus = power_gen( a = c, name="power_law_minus")


    W_minus = power_law_minus( ceil = c, tao = beta).rvs( size = n)
    W_plus = a * W_minus ** d

    # generate the bi-degree sequence sample from W+ and W-
    D_minus = np.array([])
    D_plus = np.array([])

    for w in zip(W_plus, W_minus):
        d_plus = poisson(mu = w[0]).rvs(size = 1)
        d_minus = poisson(mu = w[1]).rvs(size = 1)
        D_plus = np.append(D_plus, d_plus)
        D_minus = np.append(D_minus, d_minus)
    return (D_plus, D_minus)

def stat_test( arry ):
    print(describe(arry))
    plt.hist(arry, bins = 'auto')
    plt.tittle("Degree Histogram")
    plt.show()

seq = bi_degree(1, 1.5, 2, 2, 10000)
stat_test( seq[0] )
stat_test( seq[1] )



