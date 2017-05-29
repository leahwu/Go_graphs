from PowerLawDistribution import *
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import random

def Directed_gen(a, d, alpha, c, n):
    """
    a function that generates same-sum sample pair of bi-degree using the algorithm
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