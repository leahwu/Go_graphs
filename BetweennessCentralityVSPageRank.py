import CDMGenerator as cdm_g
import networkx as nx
import operator

# define comparison function
def bc_vs_pk(D):
    """
    
    :param D: the simple DCM graph 
    """
    bc = nx.betweenness_centrality(D)
    pr = nx.pagerank(D)

    bc_sorted = sorted(bc.items(), key=operator.itemgetter(1), reverse=True)
    pr_sorted = sorted(pr.items(), key=operator.itemgetter(1), reverse=True)

    bc_sorted_top20 = bc_sorted[0: 20]
    pr_sorted_top20 = pr_sorted[0: 20]

    print(bc_sorted_top20)
    print(pr_sorted_top20)

# Homogenuous bi-degree
a = 1
b = 2
alpha = 2
beta = 2

n = 2000

dcm_era = cdm_g.CDMGenerator(a,alpha, beta, n,'Erased')
dcm_rep = cdm_g.CDMGenerator(a,alpha, beta, n,'Erased')

D_era = dcm_era.graph
D_rep = dcm_rep.graph

bc_vs_pk(D_era)

# Inhomogenuous bi-degree
beta = 2.5  # setting beta to 2.5

indcm_era = cdm_g.CDMGenerator(a,alpha, beta, n,'Erased')

inD_era =  indcm_era.graph

bc_vs_pk(inD_era)

# consider the enhanced out-degree case (a = 1.1 & beta = 2.5)
indcm2_era = cdm_g.CDMGenerator(1.01,alpha, beta, n,'Erased')

inD2_era =  indcm2_era.graph

bc_vs_pk(indcm2_era)