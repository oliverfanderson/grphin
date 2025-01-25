import networkx as nx
import csv

def process_edges(
    file_path, G, visited_nodes, label
):
    """Helper function to process edges and add them to the graph."""
    print(f"currently processing {label} edges")
    with open(file_path, "r") as file:
        csv_reader = csv.reader(file)
        next(csv_reader)  # Skip header
        node_count = len(visited_nodes)
        edge_count = G.number_of_edges()

        for row in csv_reader:
            id1 = row[0]
            id2 = row[1]

            if id1 not in visited_nodes:
                visited_nodes.add(id1)
                node_count += 1

            if id2 not in visited_nodes:
                visited_nodes.add(id2)
                node_count += 1

            G.add_edge(id1, id2, label=label)
            edge_count += 1

def read_csv(
    ppi_path,
    reg_path,
):
    """Reads CSV files and constructs a graph with edges labeled as 'ppi' or 'reg'."""
    G = nx.MultiDiGraph()
    visited_nodes = set()

    process_edges(
        ppi_path,
        G,
        visited_nodes,
        "ppi"
    )
    process_edges(
        reg_path,
        G,
        visited_nodes,
        "reg"
    )
    print()

    return G

def simplify_graph_to_undirected(G):
    """
    Simplify the given graph to a new undirected graph with only one edge per connected node pair.

    Parameters:
    G (networkx.Graph): The original graph (can be directed or undirected, multigraph or simple).

    Returns:
    G_prime (networkx.Graph): The simplified undirected graph.
    """
    # Create an empty undirected graph
    G_prime = nx.Graph()

    # Add all nodes to G_prime
    G_prime.add_nodes_from(G.nodes())

    # Add a single undirected edge for every connected node pair in G
    for u, v in G.edges():
        G_prime.add_edge(u, v)

    return G_prime

species_dict = {
    "txid6239" : "elegans",
    "txid7227" : "fly",
    "txid7955" : "drerio",
    "txid224308" : "bsub",
    "txid55929" : "cerevisiae"
}

def main():
    txid = "txid224308"  # Replace with your species ID (e.g., txid6239 for E. coli, txid72
    species = species_dict[txid]

    ppi_path = f"data/{species}_ppi.csv"
    reg_path = f"data/{species}_reg.csv"
    output_dir = f"data/oxidative_stress/{txid}/"
    restart_nodes_file = f"{output_dir}{txid}-stress-proteins.csv"
    out_file = f"{output_dir}{txid}_pageRank.txt"


    G = read_csv(
        ppi_path,
        reg_path
    )

    G_prime = simplify_graph_to_undirected(G)
    
    # Example: Print basic graph stats
    print(f"Original graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    print(f"Simplified graph: {G_prime.number_of_nodes()} nodes, {G_prime.number_of_edges()} edges")

    # Read the list of restart nodes (proteins) from the file, skipping the header
    with open(restart_nodes_file, "r") as f:
        restart_nodes = [
            line.strip().strip('"') for i, line in enumerate(f) 
            if i > 0 and line.strip()  # Skip the first line (header) and empty lines
        ]

    print(f"Restart nodes: {restart_nodes}")

    # Create the personalization vector
    personalization = {node: (1 if node in restart_nodes else 0) for node in G_prime.nodes()}
    
    # Normalize the personalization vector to sum to 1
    total = sum(personalization.values())
    if total > 0:
        personalization = {node: value / total for node, value in personalization.items()}
    else:
        raise ValueError("Personalization vector is all zeros. Ensure restart_nodes contains valid nodes from the graph.")
    
    # Compute personalized PageRank
    pagerank_scores = nx.pagerank(G_prime, alpha=0.85, personalization=personalization)
    
    # Write PageRank scores to a file
    with open(out_file, "w") as f:
        for node, score in sorted(pagerank_scores.items(), key=lambda x: -x[1]):
            f.write(f"{node}\t{score}\n")
    
    print(f"PageRank scores written to {out_file}")


if __name__ == "__main__":
    main()