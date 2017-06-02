import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import poisson, describe


def generate_w(c, beta):
    # generate r.v. w ~ (x/c)^{-\beta}
    u = np.random.uniform(0, 1)
    w = c * (1 - u)**(-1.0/beta)
    return w


class PowerLaw:
    """
    a distribution class to specifically generate the bi-degree sequence 
    from double power law distribution with intrinsic dependence.
    
    """

    def __init__(self, a, alpha, beta, b = 2):
        # special case when alpha equals to beta, then a must be 1 and W+ equals to W-
        # b equals to c, could be assigned with any positive values, we default it as b = 2
        if alpha == beta:
            self.a = 1
            self.b = b
        else:
            self.a = a # floor for W+
            # calculating the floor for W+ b using alpha and beta
            self.b = (alpha / (alpha - 1) * (beta - 1) / beta * a ** (alpha / beta)) ** (beta / (alpha - beta))

        self.alpha = alpha  # power for W+
        self.beta = beta  # power for W-

        self.d = beta / alpha  # power dependence between W+ and W-
        self.c = (b / a) ** (alpha / beta)  # W- ceiling

    def gene_one_pair(self):
        """
        method to generate one sample (d+, d-)
        :return: tuple with a pair of d+ and d- 
        """
        w_minus = generate_w(self.c, self.beta)
        w_plus = self.a * w_minus ** self.d
        # derive the sample degree from W+ and W-
        d_plus = int(poisson(mu=w_plus).rvs())
        d_minus = int(poisson(mu=w_minus).rvs())
        return d_plus, d_minus

    def rvs(self, n):
        """
        # method to generate n sample pairs from the distribution
        :param n: number of sample pairs
        :return: tuple of n sample sequence pairs 
        """
        d_out = np.zeros(n, dtype=np.int)
        d_in = np.zeros(n, dtype=np.int)
        for i in range(0, n):
            pair = self.gene_one_pair()
            d_in[i] = pair[0]
            d_out[i] = pair[1]

        return d_in, d_out

def difference(bi_degree):
    """
    # a function that returns the difference between in-degree sequence and out-degree sequence
    :param bi_degree: the sample of bidegree sequence
    :return: the difference 
    """
    in_sum = sum(bi_degree[0])
    out_sum = sum(bi_degree[1])
    # denote delta_n as the difference
    delta_n = in_sum - out_sum
    return delta_n


# test for the convergence of bi_degree sequence
N = range(1000,6000,1000)
a = 1
alpha = 3

def test_converg_speed(alpha, beta, N,  rep = 10, a = 1):
    power_law = PowerLaw(a, alpha, beta)
    relative_diff = []
    for n in N:
        # generate the sequence with sample size n
        diff_arr = np.zeros(rep)
        for i in range(rep):
            bi_seq = power_law.rvs(n)
            diff_arr[i] = difference(bi_seq)

        # get the average sequence degree difference
        diff = np.average(diff_arr)
        # use the relative with n to measure the convergence
        relative_diff += [diff/n]
    return relative_diff




