import pandas as pd
import numpy as np
import networkx as nx
from AnalysisToolFunction import test_shortpath_marginal_2, plot_ave

data = pd.read_csv('wiki-Vote.txt', header=3, sep='\t')
data.columns=['FromNodeId', 'ToNodeId']
edges = []
l = len(data.FromNodeId)
for i in range(0, l):
    edges.append((data.FromNodeId[i], data.ToNodeId[i]))

graph = nx.DiGraph()
graph.add_edges_from(edges)


result = test_shortpath_marginal_2(graph, 5)