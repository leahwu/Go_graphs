import DCMGenerator as dg
import networkx as nx
import matplotlib.pyplot as plt
import operator

dcm = dg.CDMGenerator(1, 3, 3.5, 2000, 'Erased')
dcm.pr_vs_bc_plot()
bc = dcm.betweeness_centrality
pr = dcm.page_rank
pr_sort = sorted(pr.items(), key=operator.itemgetter(1), reverse=True)
pr_scores = [node[1] for node in pr_sort]
bc_scores = [bc[node[0]] for node in pr_sort]

# scale the bc_scores that it sums to 1 in comparison to page-rank
bc_scaled_scores = [elem / sum(bc_scores) for elem in bc_scores]

plt.figure(2)
bc_100 = bc_scaled_scores[0:100]
pr_100 = pr_scores[0:100]
plt.plot(pr_100, 'ro', markersize=2)
plt.plot(bc_100, 'bx', markersize=2)
plt.legend(['Page_Rank', 'Betweeness_Centrality'])
plt.xlabel('Node')
plt.ylabel('Ranking')
plt.title('Top100 nodes comparison between Pagerank and Betweeness centrality')