import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import poisson
import DCMGenerator as dcm_g


def generate_w(c, beta):
    # generate r.v. w ~ (x/c)^{-\beta} as tail distribution
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
        self.alpha = beta / d

        if d == 1:
            self.b = 10 # set default values
        else:
            self.b = (self.alpha / (self.alpha - 1) * (beta - 1) / beta * a ** (self.alpha / beta)) ** (beta / (self.alpha - beta))

        self.c = (self.b / self.a) ** (self.alpha / beta)

        print("The parameters are:")
        print("a = ", self.a)
        print("b = ", self.b)
        print("c = ", self.c)
        print("d = ", self.d)
        print("alpha = ", self.alpha)
        print("beta = ", self.beta)

        self.e_w_minus = self.c * beta/ (self.beta - 1)
        self.e_w_plus = a*beta*self.c**d / (beta - d)

        self.params = {'a': self.a, 'd': self.d, 'beta': self.beta, 'alpha': self.alpha,
                       'b': self.b, 'c': self.c}

        # print("E[W^minus] = ", self.e_w_minus)
        # print("E[W^plus] = ", self.e_w_plus)

    def gene_iid_pairs_of_w(self, n):

        w_minus = np.zeros(n)
        w_plus = np.zeros(n)

        for i in range(0, n):
            w_minus[i] = generate_w(self.c, self.beta)
            w_plus[i] = self.a * w_minus[i] ** self.d

        return w_plus, w_minus

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

        # print("mean of in-degree sequence: ", np.mean(d_in))
        # print("mean of out-degree sequence: ", np.mean(d_out))

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
    a = 1.2  # the lower bound for W+
    alpha = 3  # the power for W+
    beta = 3.5  # the power for W+

    d = beta / alpha

    n = 5000  # graph size

    pld = PowerLaw(a, d, beta)
    din, dout = pld.rvs(n)

    w_minus, w_plus = pld.gene_iid_pairs_of_w(n)
    print("mean of w_minus is ", np.mean(w_minus))
    print("mean of w_plus is ", np.mean(w_plus))

    print("mean of in-degree sequence is ", np.mean(din))
    print("mean of out-degree sequence is ", np.mean(dout))


ww = []
poisson_rvs = []
for i in range(0, 2000):
    w = generate_w(2, 3)
    rv = poisson.rvs(mu=w)
    ww.append(w)
    poisson_rvs.append(rv)

ww = np.array(ww)
poisson_rvs = np.array(poisson_rvs)



