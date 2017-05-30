import CDMGenerator as cdm_g
import networkx as nx
import operator

a = 1
b = 2
alpha = 2
beta = 2

n = 2000

dcm = cdm_g.CDMGenerator(a,alpha, beta, n,'Erased')

D = dcm.graph

bc = nx.betweenness_centrality(D)
pr = nx.pagerank(D)

bc_sorted = sorted(bc.items(), key=operator.itemgetter(1),reverse=True)
pr_sorted = sorted(pr.items(), key=operator.itemgetter(1), reverse=True)

bc_sorted_top20 = bc_sorted[0 : 20]
pr_sorted_top20 = pr_sorted[0 : 20]

print(bc_sorted_top20)
print(pr_sorted_top20)