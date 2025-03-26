import networkx as nx
import matplotlib.pyplot as plt
import random
from itertools import combinations
from collections import Counter

def generate_random_multi_graph(n, edge_probability=0.3, edge_label_probability = 0.4):

    G = nx.MultiDiGraph()
    G.add_nodes_from(range(n))
    for i in range(n):
        for j in range(i + 1, n):
            if random.random() < edge_probability:
                if random.random() <edge_label_probability:
                    G.add_edge(i, j, label="ppi")
                else:
                    G.add_edge(i, j, label="reg")

    return G

def enumerate_labeled_subgraphs(graph, size):

    # extremely inefficient
    nodes = list(graph.nodes)
    subgraphs = []
    for node_set in combinations(nodes, size):
        subgraph = graph.subgraph(node_set).copy()
        if len(subgraph.edges) != 0:
            for u, v, key, data in graph.edges(keys=True, data=True):
                if u in node_set and v in node_set:
                    subgraph.add_edge(u, v, key=key, **data)
            subgraphs.append(subgraph)
    return subgraphs

def classify_labeled_graphlet(graphlet):
    edges = len(graphlet.edges)
    nodes = len(graphlet.nodes)

    # need complete list of classifiers of graphlets
    if nodes == 2:
        edge_labels = [data["label"] for _, _, data in graphlet.edges(data=True)]
        return f"Edge-{edge_labels[0]}"
    elif nodes == 3:
        edge_labels = [data["label"] for _, _, data in graphlet.edges(data=True)]
        labels = "-".join(sorted(edge_labels)) 
        if edges == 3:
            return f"Triangle-{labels}"
        elif edges == 2:
            return f"Path-{labels}"
        elif edges == 1:
            return f"Disconnected-{labels}"
        else:  # edges == 0
            return "Isolated"

def calculate_labeled_graphlet_distribution(graph, max_size=3):
    distribution = Counter()
    for size in range(2, max_size + 1):
        subgraphs = enumerate_labeled_subgraphs(graph, size)
        print(f"{size}-node graphlets count:{len(subgraphs)}")
        for subgraph in subgraphs:
            graphlet_type = classify_labeled_graphlet(subgraph)
            distribution[graphlet_type] += 1
    return distribution

def main():
    edge_color_map = {
        "ppi": "blue",
        "reg": "red"
    }
    mixed_network = generate_random_multi_graph(50)

    distribution = calculate_labeled_graphlet_distribution(mixed_network)

    labels, values = zip(*distribution.items())
    plt.figure(figsize=(10, 5))
    plt.bar(labels, values, color='skyblue')
    plt.xticks(rotation=90)
    plt.title("Labeled Graphlet Distribution (2-3 Node Motifs)")
    plt.xlabel("Graphlet Type with Labels")
    plt.ylabel("Count")
    plt.tight_layout()
    plt.show()

    plt.show()


if __name__ == "__main__":
    main()