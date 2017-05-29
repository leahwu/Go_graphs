import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import random
from math import exp, factorial
from scipy.stats import rv_continuous, rv_discrete, poisson, describe, histogram, randint


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
    plt.title("Degree Histogram")
    plt.show()

def Directed_gen(a, d, alpha, c, n):
    """
    a function that generates same-sum sample pair of bidegree using the algorithm
    the function returns a tuple  (in-degree, out-degree)

    Input Parameters
    -----------------
    a : coefficient
    d : power term
    alpha : power for W+ CDF
    c : W- ceiling
    n: simulation sample number

    """
    beta = alpha * d
    # derive k
    k = min(1 - 1/alpha, 1 - 1/beta, 1/2)
    delta0 = k/2  # take k/2 as fixed

    Delta_n = sum(bi_degree(a, d, alpha, c, n)[0]) - sum(bi_degree(a, d, alpha, c, n)[1])

    # repeat until Delta_n is small enough to be "negligible"
    while Delta_n > n**(1 - k + delta0):
        samesum_bi_degree = bi_degree(a, d, alpha, c, n)
        Delta_n = sum(samesum_bi_degree[0]) - sum(samesum_bi_degree[1])

    if abs(Delta_n) > 0:
        # derive random sample nodes deltai
        delat_i = random.sample(range(0, n), Delta_n)
        for i in delat_i:
            if Delta_n > 0:
                samesum_bi_degree[1][i] + 1 # D_minus random degree increase 1
            else:
                samesum_bi_degree[0][i] + 1 # D_plus random degree increse 1
        print("lxh")

    return samesum_bi_degree







