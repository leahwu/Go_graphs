from collections import Counter
import SDGErased as SDG
import DCMGenerator as cdm_g
import matplotlib.pyplot as plt
import time

start = time.clock()

# assign the parameters
# [Special Case [Alpha = Beta]]
a = 1 # the floor for W+
alpha = 3 # the power for W+
beta = 3 # the power for W+
b = 2.5

n = 2000 # simulation times



# generate simple directed configuration model
dcm = cdm_g.CDMGenerator(a,alpha, beta, n, 'Erased')
D = dcm.graph

# plot the bi-degree distribution with generated simple DCM degree distribution
dcm.degrees_plot()