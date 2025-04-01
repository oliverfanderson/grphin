import networkx as nx
import csv
from matplotlib import pyplot as plt
from pyvis.network import Network

species_dict = {
    "txid6239" : "elegans",
    "txid7227" : "fly",
    "txid7955" : "drerio",
    "txid224308" : "bsub",
    "txid559292" : "cerevisiae"
}

def process_edges(
    file_path, G, visited_nodes, label
):
    """Helper function to process edges and add them to the graph."""
    # print(f"currently processing {label} edges")
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

def threshold_stress_recovery(pagerank, stress_proteins, threshold):
    """
    Thresholds PageRank results by percent of stress proteins recovered.
    
    Parameters:
    pagerank (dict): Dictionary of node IDs and their PageRank scores.
    stress_proteins (list): List of protein IDs for stress response.
    threshold (float): The threshold value (e.g., 0.8 for 80%).

    Returns:
    pagerank_out (dict): A dictionary of nodes with PageRank scores greater than or equal to the threshold score.
    """
    # Filter PageRank scores for stress proteins only
    stress_scores = {node: pagerank[node] for node in stress_proteins if node in pagerank}

    # Sort stress proteins by PageRank score in descending order
    sorted_stress_scores = sorted(stress_scores.items(), key=lambda x: -x[1])

    # Determine the number of stress proteins to include (top 80%)
    total_stress_proteins = len(stress_scores)
    cutoff_index = int(total_stress_proteins * threshold) - 1

    # Get the threshold score from the top 80% of stress proteins
    if cutoff_index >= 0 and cutoff_index < len(sorted_stress_scores):
        threshold_score = sorted_stress_scores[cutoff_index][1]
    else:
        raise ValueError("Threshold calculation failed. Ensure the stress proteins list is non-empty.")
    
    print(f"Threshold score for top {threshold * 100:.0f}% stress proteins: {threshold_score:.4f}")

    # Filter all nodes in PageRank by the threshold score
    pagerank_out = {node: score for node, score in pagerank.items() if score >= threshold_score}

    return pagerank_out

def get_induced_subnetwork(G, pagerank_out):
    """
    Extracts the induced subnetwork based on the given PageRank scores.

    Parameters:
    G (networkx.Graph): The original graph.
    pagerank_out (dict): The PageRank scores for nodes to include in the subnetwork.
    """
    # Create an empty graph to store the induced subnetwork
    G_induced = nx.MultiDiGraph()

    # Add nodes with PageRank scores above the threshold to the induced subnetwork
    G_induced.add_nodes_from(pagerank_out.keys())

    # Add edges between nodes with PageRank scores above the threshold to the induced subnetwork
    for u, v, key, data in G.edges(data=True, keys=True):
        if u in pagerank_out and v in pagerank_out:
            G_induced.add_edge(u, v, key=key, **data)

    return G_induced

def split_to_csv(G_induced, output_dir):
    """
    Writes the graph to CSV files based on edge labels.

    Parameters:
    G_induced (networkx.Graph): The induced graph.
    ppi_input (set): A set of edges for PPI (tuples of (u, v)).
    reg_input (set): A set of edges for REG (tuples of (u, v)).
    output_dir (str): Directory where output CSV files will be saved.
    """

    # File paths for output
    ppi_file = f"{output_dir}stress_ppi.csv"
    reg_file = f"{output_dir}stress_reg.csv"

    # Write edges to CSV files
    with open(ppi_file, "w", newline="") as ppi_out, open(reg_file, "w", newline="") as reg_out:
        ppi_writer = csv.writer(ppi_out, quotechar='"', quoting=csv.QUOTE_ALL)
        reg_writer = csv.writer(reg_out, quotechar='"', quoting=csv.QUOTE_ALL)
        
        # Write CSV headers
        ppi_writer.writerow(["id1", "id2"])
        reg_writer.writerow(["id1", "id2"])
        
        # Iterate over edges
        for u, v, key, data in G_induced.edges(data=True, keys=True):
            # print(f"Edge ({u}, {v}, {key}): {data}")
            label = data.get("label", None)
            if label == "ppi":
                ppi_writer.writerow([u, v])
            elif label == "reg":
                reg_writer.writerow([u, v])

    print(f"PPI edges written to: {ppi_file}")
    print(f"Reg edges written to: {reg_file}")


def main():
   # List of taxon IDs to process
    taxon_ids = ["txid6239", "txid7227", "txid7955", "txid224308", "txid559292"]
    # taxon_ids = ["txid224308"]
    threshold = 0.8

    for txid in taxon_ids:
        species = species_dict.get(txid)

        ppi_path = f"data/{species}_ppi.csv"
        reg_path = f"data/{species}_reg.csv"
        output_dir = f"data/oxidative_stress/{txid}/"
        restart_nodes_file = f"{output_dir}{txid}-stress-proteins.csv"
        out_file = f"{output_dir}{txid}_pageRank.txt"
        threshold_file = f"{output_dir}{txid}_pageRank{threshold*100}.txt"


        G = read_csv(
            ppi_path,
            reg_path
        )

        G_prime = simplify_graph_to_undirected(G)
        
        # # Example: Print basic graph stats
        # print(f"Original graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
        # print(f"Simplified graph: {G_prime.number_of_nodes()} nodes, {G_prime.number_of_edges()} edges")

        # Read the list of restart nodes (proteins) from the file, skipping the header
        with open(restart_nodes_file, "r") as f:
            restart_nodes = {
                line.strip().strip('"') for i, line in enumerate(f) 
                if i > 0 and line.strip()  # Skip the first line (header) and empty lines
            }

        # Convert the set back to a list if needed
        restart_nodes = list(restart_nodes)

        G_sub = nx.induced_subgraph(G, restart_nodes)
        G_sub = nx.Graph(G_sub)
        print(f"Number of connected components without RWR: {nx.number_connected_components(G_sub)}")

        print(f"Total Stress Proteins: {len(restart_nodes)}")

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

        # Threshold the PageRank results by percent of stress proteins recovered
        pagerank_out = threshold_stress_recovery(pagerank_scores, restart_nodes, threshold)  # Example usage: threshold at 80% of stress proteins recovered 
        
        # Write PageRank scores to a file
        with open(threshold_file, "w") as f:
            for node, score in sorted(pagerank_out.items(), key=lambda x: -x[1]):
                f.write(f"{node}\t{score}\n")
        
        print(f"Thresholded PageRank scores written to {threshold_file}")

        G_induced = get_induced_subnetwork(G, pagerank_out)
        
        print(f"Induced subnetwork: {G_induced}")

        print(f"Number of connected components after RWR: {nx.number_connected_components(G_induced)}")

        split_to_csv(G_induced, output_dir)


        # print(f"Total degree of stress subnetwork: {sum(dict(G_induced.degree()).values())}")
        # Have to subtract the above number - len("stress_ppi.csv") to get mixed degree

        # nt = Network('500px', '600px', directed=True)
        # G_induced.remove_edges_from(nx.selfloop_edges(G_induced))

        # nt.barnes_hut()  # Improve layout

        # for node in G_induced.nodes():
        #     nt.add_node(node)
        #     nt.get_node(node)['color'] = '#cfe2f3'
        #     nt.get_node(node)['size'] = 100  # Adjust size if needed

        # # Track added edges to prevent duplication
        # seen_ppi_edges = set()

        # # Iterate over edges
        # for u, v, data in G_induced.edges(data=True):
        #     label = data.get("label", None)  # Get label from edge data
        #     if label == "reg":
        #         nt.add_edge(u, v, color="red", width=20, arrows="to", arrowStrikethrough=True)  # Directed red edge
        #     elif label == "ppi":
        #         if (u, v) not in seen_ppi_edges and (v, u) not in seen_ppi_edges:
        #             nt.add_edge(u, v, color="black", width=20, arrows="")  # Undirected black edge
        #             seen_ppi_edges.add((u, v))  # Mark the pair as seen
        #             seen_ppi_edges.add((v, u))  # Ensure (v, u) isn't added separately

        # nt.set_edge_smooth('dynamic')
        # nt.write_html(f'{output_dir}nx.html')

if __name__ == "__main__":
    main()