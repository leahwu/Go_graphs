import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import poisson



def generate_w(c, beta):
    # generate r.v. w ~ (x/c)^{-\beta}
    u = np.random.uniform(0, 1)
    w = c * (1 - u)**(-1.0/beta)
    return w



class PowerLaw:
    """
    a distribution class to specifically generate the bi-degree sequence 
    from double power law distribution with intrinsic dependence.
    
    We need beta > 2, beta > 2 * d, beta > d + 1 for finite variance, and beta > 1 for zero mean.
    
    """

    def __init__(self, a, d, beta):
        # special case when alpha equals to beta, then a must be 1 and W+ equals to W-
        # b equals to c, could be assigned with any positive values, we default it as b = 2
        self.a = a
        self.d = d
        self.beta = beta
        self.b = a**((beta-1)/(beta-1-(beta-d)*d)) * ((beta - d)/((beta-1)*a))**(d/(beta -1 - (beta-d)*d))

        self.c = (self.b / a)**(1/d)
        self.alpha = beta / d
        print("The parameters are:")
        print("a = ", self.a)
        print("b = ", self.b)
        print("c = ", self.c)
        print("d = ", self.d)
        print("alpha = ", self.alpha)
        print("beta = ", self.beta)


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
    # print("in_sum", end=":")
    # print(in_sum)
    out_sum = sum(bi_degree[1])
    # print("out_sum", end=":")
    # print(out_sum)
    # denote delta_n as the difference
    delta_n = in_sum - out_sum
    return delta_n


# test for the convergence of bi_degree sequence
#N = range(1000,6000,1000)
#a = 2
#alpha = 3

def test_converg_speed(a, d, beta, N,  rep = 10):
    power_law = PowerLaw(a, d, beta)
    relative_diff = []
    for n in N:
        # generate the sequence with sample size n
        diff_arr = np.zeros(rep)
        for i in range(rep):
            bi_seq = power_law.rvs(n)
            diff_arr[i] = difference(bi_seq)
            print(i, diff_arr[i])

        # get the average sequence degree difference
        diff = np.average(diff_arr)
        # use the relative with n to measure the convergence
        print(diff, n)
        relative_diff += [diff/n]
    return relative_diff

def test():
    a = 2  # the lower bound for W+
    alpha = 5  # the power for W+
    beta = 6  # the power for W+
    b = 2.5  # if alpha != beta, we could choose any b = c
    d = beta / alpha

    n = 2000  # graph size

    test_converg_speed(a, d, beta,[2000])





