import numpy as np
from scipy.stats import rv_continuous, rv_discrete, poisson, describe, histogram, randint


class GenerateW(rv_continuous):
    # generate W power-law distribution
    def _pdf(self, x, floor, tao ):
        return (tao/floor) * (x / floor) ** (-tao-1)


class PowerLaw:
    """
    a distribution class to specifically generate the bi-degree sequence 
    from double power law distribution with intrinsic dependence.
    
    """

    def __init__(self, a, b, alpha, beta):
        self.a = a # floor for W+
        self.b = b # floor for W-
        self.alpha = alpha # power for W+
        self.beta = beta # power for W-

        self.d = beta / alpha  # power dependence between W+ and W-
        self.c = (b / a) ** (alpha / beta)  # W- ceiling
        self.wgenerator = GenerateW( a=self.c, name="wnegerator")

    # method to generate one sample (d+, d-)
    def gene_one_pair(self):
        w_minus = self.wgenerator(floor=self.c, tao=self.beta).rvs()
        w_plus = self.a * w_minus ** self.d
        # derive the sample degree from W+ and W-
        d_plus = int(poisson(mu=w_plus).rvs())
        d_minus = int(poisson(mu=w_minus).rvs())
        return d_plus, d_minus

    # method to generate n sample pairs from the distribution
    def iid(self, n):
        d_out = np.zeros(n, dtype=np.int)
        d_in = np.zeros(n, dtype=np.int)
        for i in range(0, n):
            pair = self.gene_one_pair()
            d_in[i] = pair[0]
            d_out[i] = pair[1]

        return d_in, d_out


