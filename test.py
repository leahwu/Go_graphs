from ValidDegree import *
from PowerLawDistribution import *
import SimpleDirectedGraph as sdg

a = 1
b = 2
alpha = 2
beta = 2

n = 2000

fg = PowerLaw(a, b, alpha, beta)


bi_seq = directed_gen(alpha, beta, fg, n)

D = sdg.gen_simple_DCM(bi_seq)
sdg.plot_graph(D)

