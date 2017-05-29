import networkx as nx
import matplotlib.pyplot as plt



def gen_simple_DCM( bi_degree ):
    d_in = bi_degree[0].tolist()
    d_out = bi_degree[1].tolist()
    dcm = nx.directed_configuration_model(d_in, d_out)

    # remove parallel edges
    dcm = nx.DiGraph(dcm)
    # remove self-loops
    dcm.remove_edges_from(dcm.selfloop_edges())

    return dcm


def plot_graph(graph):
    degree_sequence=sorted(nx.degree(graph).values(), reverse=True)
    # plot the degree rank
    dmax = max(degree_sequence)
    plt.loglog(degree_sequence, 'b-', marker='o')
    plt.title("Degree rank plot")
    plt.ylabel("degree")
    plt.xlabel("rank")

    plt.axes([0.45, 0.45, 0.45, 0.45])
    Gcc=sorted(nx.connected_component_subgraphs(graph), key=len, reverse=True)[0]
    pos=nx.spring_layout(Gcc)
    plt.axis('off')
    nx.draw_networkx_nodes(Gcc,pos,node_size=20)
    nx.draw_networkx_edges(Gcc,pos,alpha=0.4)
    # plot the graph
    plt.axes([0.45,0.45,0.45,0.45])
    Gcc=sorted(nx.connected_component_subgraphs(graph), key=len, reverse=True)[0]
    pos=nx.spring_layout(Gcc)
    plt.axis('off')
    nx.draw_networkx_nodes(Gcc, pos, node_size=20)
    nx.draw_networkx_edges(Gcc, pos, alpha=0.4)
