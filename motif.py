import networkx as nx
import matplotlib.pyplot as plt
import random

def generate_random_multi_graph(n, edge_probability=0.3, edge_label_probability = 0.4):

    G = nx.MultiDiGraph()
    G.add_nodes_from(range(n))
    for i in range(n):
        for j in range(i + 1, n):  # Avoid self-loops
            if random.random() < edge_probability:
                if random.random() <edge_label_probability:
                    G.add_edge(i, j, label="ppi")
                else:
                    G.add_edge(i, j, label="reg")

    return G


def main():
    print("Hello there")
    edge_color_map = {
        "ppi": "blue",
        "reg": "red"
    }
    P = generate_random_multi_graph(50)

    plt.figure(figsize=(10, 8))

    pos = nx.spring_layout(P)
    nx.draw_networkx_nodes(P, pos, node_color="lightblue", node_size=100)

    for edge_type, color in edge_color_map.items():
        edges = [(u, v) for u, v, d in P.edges(data=True) if d["label"] == edge_type]
        nx.draw_networkx_edges(P, pos, edgelist=edges, edge_color=color, label=edge_type, width=2)
    plt.show()


if __name__ == "__main__":
    main()