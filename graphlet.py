from matplotlib import pyplot as plt
import networkx as nx
import numpy
import random
import scipy as sp


def get_two_node_hash_table():
    print("hash function")

    # 2-node graphlet binary edge vectors
    # see diagram for reference
    a_1 = (1, 0, 0)
    a_2 = (1, 0, 0)
    a_hash = hash(a_1 + a_2)

    b_1 = (1, 1, 0)
    b_2 = (1, 0, 1)
    b_hash = hash(b_1 + b_2)

    c_1 = (1, 1, 1)
    c_2 = (1, 1, 1)
    c_hash = hash(c_1 + c_2)

    d_1 = (0, 1, 0)
    d_2 = (0, 0, 1)
    d_hash = hash(d_1 + d_2)

    e_1 = (0, 1, 1)
    e_2 = (0, 1, 1)
    e_hash = hash(e_1 + e_2)

    two_node_hash_table = {a_hash: 0, b_hash: 0, c_hash: 0, d_hash: 0, e_hash: 0}
    return two_node_hash_table


def generate_random_multi_graph(n, edge_probability=0.3, edge_label_probability=0.4):
    G = nx.MultiDiGraph()
    G.add_nodes_from(range(n))
    for i in range(n):
        for j in range(i + 1, n):
            if random.random() < edge_probability:
                if random.random() < edge_label_probability:
                    G.add_edge(i, j, label="ppi")
                else:
                    G.add_edge(i, j, label="reg")
            # want to generate a graph that is connected
            elif G.degree(j) == 0:
                if random.random() < edge_label_probability:
                    G.add_edge(j, i, label="ppi")
                else:
                    G.add_edge(j, i, label="reg")
    return G


def draw_labeled_multigraph(G, attr_name, ax=None):
    """
    Length of connectionstyle must be at least that of a maximum number of edges
    between pair of nodes. This number is maximum one-sided connections
    for directed graph and maximum total connections for undirected graph.
    """

    pos = nx.shell_layout(G)
    nx.draw_networkx_nodes(G, pos, ax=ax)
    nx.draw_networkx_labels(G, pos, font_size=20, ax=ax)
    nx.draw_networkx_edges(G, pos, edge_color="grey", ax=ax)

    labels = {
        tuple(edge): f"{attr_name}={attrs[attr_name]}"
        for *edge, attrs in G.edges(keys=True, data=True)
    }
    nx.draw_networkx_edge_labels(
        G,
        pos,
        labels,
        label_pos=0.3,
        font_color="blue",
        bbox={"alpha": 0},
        ax=ax,
    )


def main():
    two_node_hash_table = get_two_node_hash_table()
    G = generate_random_multi_graph(50)
    G_adj_matrix = nx.adjacency_matrix(G)

    adj_matrix = [
        [[0, 0, 0] for _ in range(len(G.nodes()))]
        for _ in range(len(G.nodes()))
    ]

    for edge in G.edges(data=True):
        # print(edge)
        i = edge[0]
        j = edge[1]
        edge_type = edge[2]["label"]

        if edge_type == "ppi":
            adj_matrix[i][j][0] = 1
            adj_matrix[j][i][0] = 1
        elif edge_type == "reg":
            adj_matrix[i][j][1] = 1
            adj_matrix[j][i][2] = 1

    for i in range(len(adj_matrix)):
        for j in range(len(adj_matrix[0])): 
            vector = adj_matrix[i][j] + adj_matrix[j][i]
            if hash(tuple(vector)) in two_node_hash_table:
                two_node_hash_table[hash(tuple(vector))]+=1

    print(two_node_hash_table)
    draw_labeled_multigraph(G, "label")
    plt.show()


if __name__ == "__main__":
    main()
