import SDGErased as sdg_e
import PowerLawDistribution as pld
import ValidDegree as vd
import matplotlib.pyplot as plt
import SDGRepeated as sdg_r
from ValidDegree import *


class CDMGenerator(object):

    def __init__(self, a, alpha, beta, n, algorithm, b=2):
        fg = pld.PowerLaw(a, alpha, beta, b)
        self.bi_seq = directed_gen(alpha, beta, fg, n)
        if algorithm == 'Erased':
            self.graph = sdg_e.gen_simple_DCM(self.bi_seq)
        if algorithm == 'RevisedRepeated':
            self.graph = sdg_r.gen_simple_DCM2(self.bi_seq)
        if algorithm == 'Repeated':
            self.graph = sdg_r.gen_simple_DCM(self.bi_seq)


