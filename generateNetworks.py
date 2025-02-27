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

def split_to_csv(G_prime, out_ppi_path, out_reg_path):
    """
    Writes the graph to CSV files based on edge labels.

    Parameters:
    G_induced (networkx.Graph): The induced graph.
    ppi_input (set): A set of edges for PPI (tuples of (u, v)).
    reg_input (set): A set of edges for REG (tuples of (u, v)).
    output_dir (str): Directory where output CSV files will be saved.
    """

    # Write edges to CSV files
    with open(out_ppi_path, "w", newline="") as ppi_out, open(out_reg_path, "w", newline="") as reg_out:
        ppi_writer = csv.writer(ppi_out, quotechar='"', quoting=csv.QUOTE_ALL)
        reg_writer = csv.writer(reg_out, quotechar='"', quoting=csv.QUOTE_ALL)
        
        # Write CSV headers
        ppi_writer.writerow(["id1", "id2"])
        reg_writer.writerow(["id1", "id2"])
        
        # Iterate over edges
        for u, v, key, data in G_prime.edges(data=True, keys=True):
            # print(f"Edge ({u}, {v}, {key}): {data}")
            label = data.get("label", None)
            if label == "ppi":
                ppi_writer.writerow([u, v])
            elif label == "reg":
                reg_writer.writerow([u, v])

    print(f"PPI edges written to: {out_ppi_path}")
    print(f"Reg edges written to: {out_reg_path}")


def main():
   # List of taxon IDs to process
    # taxon_ids = ["txid6239", "txid7227", "txid7955", "txid224308", "txid559292"]
    taxon_ids = ["txid224308"]

    for txid in taxon_ids:
        ppi_path = f"data/oxidative_stress/{txid}/stress_ppi.csv"
        reg_path = f"data/oxidative_stress/{txid}/stress_reg.csv"
        output_dir = f"data/oxidative_stress/{txid}/randomized_networks/"
        out_ppi_path = f"{output_dir}stress_ppi{i}.csv"
        out_reg_path = f"{output_dir}stress_reg{i}.csv"

        G = read_csv(
            ppi_path,
            reg_path
        )

        G_prime = shuffle_edges(G)
        
        # # Example: Print basic graph stats
        # print(f"Original graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
        
        print(f"Shuffled subnetwork: {G_prime}")

        split_to_csv(G_prime, out_ppi_path, out_reg_path)


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