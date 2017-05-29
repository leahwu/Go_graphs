import numpy as np
from scipy.stats import rv_continuous, rv_discrete, poisson, describe, histogram, randint


class GenerateW(rv_continuous):
    # generate w
    def _pdf(self, x, floor, tao ):
        return (tao/floor) * (x / floor) ** (-tao-1)


class PowerLaw:
    a = 0
    d = 0
    c = 0
    d = 0
    beta = 0
    wgenerator = ""

    def __init__(self, a, b, alpha, beta):
        self.a = a
        self.b = b
        self.beta = beta
        self.d = beta / alpha  # power dependence between W+ and W-
        self.c = (b / a) ** (alpha / beta)  # W- ceiling
        self.wgenerator = GenerateW(a=self.c, name="wnegerator")

    # generate one (d+, d-)
    def gene_one_pair(self):
        w_minus = self.wgenerator(floor=self.c, tao=self.beta).rvs(size=1)
        w_plus = self.a * w_minus ** self.d
        d_plus = int(poisson(mu=w_plus).rvs(size=1))
        d_minus = int(poisson(mu=w_minus).rvs(size=1))
        return d_plus, d_minus

    def iid(self, n):
        d_out = np.zeros((n,), dtype=np.int)
        d_in = np.zeros((n,), dtype=np.int)
        for i in range(0, n):
            pair = self.gene_one_pair()
            d_in[i] = pair[0]
            d_out[i] = pair[1]
            return d_in, d_out


