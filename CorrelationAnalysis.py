import DCMGenerator as dcm
import numpy as np



alpha = 5
beta = 6
a_lst = list(np.arange(1,10.1,0.1))
cor_lst = []

for i in range(len(a_lst)):
    a = a_lst[i]
    graph_a = dcm.DCMGenerator(a, alpha, beta, n = 5000, algorithm = 'Erased')
    cor_lst.append(graph_a.corr_in_and_out()[0])


