import networkx as nx
import random
import csv
from matplotlib import pyplot as plt

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
    # Remove self-loops
    self_loops = list(nx.selfloop_edges(G))
    G.remove_edges_from(self_loops)

    return G

def label_edges(G):
    """
    Assigns new edge labels based on the combination of PPI and Reg edges between node pairs.
    Ensures that each node pair is processed only once.

    Parameters:
    G (networkx.MultiDiGraph): Input graph with edges labeled as 'ppi' or 'reg'.

    Returns:
    G_prime (networkx.MultiDiGraph): A new graph with relabeled edges.
    """
    G_prime = nx.MultiDiGraph()
    G_prime.add_nodes_from(G.nodes(data=True))  # Preserve node attributes

    # Step 1: Count edge types between each node pair
    edge_counts = {}
    for u, v, d in G.edges(data=True):
        edge_label = d["label"]

        # Ensure (u, v) key is in dictionary (sorted to avoid duplicates)
        key = tuple(sorted((u, v)))  # This ensures undirected edges are counted only once
        if key not in edge_counts:
            edge_counts[key] = {"ppi": 0, "reg": 0}

        # Count edges correctly
        if edge_label == "ppi":
            edge_counts[key]["ppi"] += 1  # PPI is undirected, count once
        elif edge_label == "reg":
            edge_counts[key]["reg"] += 1  # Reg is directed, count normally

    # Step 2: Assign new labels based on rules, avoiding duplicate edges
    processed_pairs = set()
    for (u, v), counts in edge_counts.items():
        if (u, v) in processed_pairs:
            continue  # Skip if already processed

        num_ppi = counts["ppi"]
        num_reg = counts["reg"]

        # Determine new edge label
        if num_ppi > 0 and num_reg == 0:
            new_label = "ppi"  # Only PPI
        elif num_ppi == 0 and num_reg > 0:
            new_label = "reg"  # Only Reg
        elif num_ppi > 0 and num_reg == 1:
            new_label = "mix"  # One Reg + One PPI
        elif num_ppi == 0 and num_reg > 1:
            new_label = "coreg"  # Two directed Reg edges
        elif num_ppi > 0 and num_reg > 1:
            new_label = "coreg_ppi"  # Two directed Reg + One PPI
        else:
            continue  # Shouldn't happen

        # Add the relabeled edge once
        G_prime.add_edge(u, v, label=new_label)
        processed_pairs.add((u, v))  # Mark as processed

    return G_prime

def swap_edges(G_prime, num_swaps=100):
    """
    Swap edge labels between randomly chosen compatible edges in the graph.

    Parameters:
    G_prime (networkx.MultiDiGraph): The labeled graph.
    num_swaps (int): The number of swap attempts.

    Returns:
    G_shuffled (networkx.MultiDiGraph): A new graph with shuffled edges.
    """
    if G_prime is None or len(G_prime.nodes) == 0:
        print("Error: G_prime is empty or None!")
        return None

    G_random = nx.MultiDiGraph()
    G_random.update(G_prime)

    if len(G_random.nodes) == 0 or len(G_random.edges) == 0:
        print("Error: G_random is empty after copying G_prime!")
        return G_random

    print(f"Initial Graph - Nodes: {len(G_random.nodes)}, Edges: {len(G_random.edges)}")

    edges = list(G_random.edges(keys=True, data=True))  # List of all edges with keys
    if len(edges) < 2:
        print("Not enough edges to perform swaps.")
        return G_random

    edge_dict = {(u, v, k): d["label"] for u, v, k, d in edges}  # Store labels

    swap_count = 0

    for _ in range(num_swaps):
        if len(edges) < 2:
            print("Not enough edges remaining for swaps.")
            break

        try:
            (u, v, k1, data1), (x, y, k2, data2) = random.sample(edges, 2)
        except ValueError:
            print("Error: Not enough unique edges available for sampling.")
            break

        # Ensure (u, v) and (x, y) are the same edge type
        label_uv = data1["label"]
        label_xy = data2["label"]

        if label_uv != label_xy:
            continue  # Skip if edge labels don't match

        # Check (u, y) and (x, v)
        uy_exists = (u, y) in G_random.edges
        xv_exists = (x, v) in G_random.edges

        if uy_exists and xv_exists:
            # Both edges exist, check if they have the same label
            k_uy = list(G_random[u][y].keys())[0]  # Get an arbitrary key
            k_xv = list(G_random[x][v].keys())[0]
            label_uy = G_random[u][y][k_uy]["label"]
            label_xv = G_random[x][v][k_xv]["label"]

            if label_uy != label_xv:
                continue  # Skip swap if labels don't match

        elif not uy_exists and not xv_exists:
            # If neither (u, y) nor (x, v) exist, treat them as the same "non-edge"
            label_uy = "none"
            label_xv = "none"

        else:
            continue  # Skip if only one of them exists

        # Perform the swap
        if uy_exists:
            G_random[u][y][k_uy]["label"], G_random[u][v][k1]["label"] = label_uv, label_uy
        else:
            G_random.add_edge(u, y, key=k1, label=label_uv)
            G_random.remove_edge(u, v, key=k1)

        if xv_exists:
            G_random[x][v][k_xv]["label"], G_random[x][y][k2]["label"] = label_xy, label_xv
        else:
            G_random.add_edge(x, v, key=k2, label=label_xy)
            G_random.remove_edge(x, y, key=k2)

        swap_count += 1
        print(f"Swapped edges: ({u}, {v}) <-> ({u}, {y}) and ({x}, {y}) <-> ({x}, {v})")

    print(f"Total swaps performed: {swap_count}")
    print(f"Final Graph - Nodes: {len(G_random.nodes)}, Edges: {len(G_random.edges)}")

    return G_random

# def swap_edges(G_prime, num_swaps=100):
#     """
#     Swap edge labels between randomly chosen compatible edges in the graph.

#     Parameters:
#     G_prime (networkx.MultiDiGraph): The labeled graph.
#     num_swaps (int): The number of swap attempts.

#     Returns:
#     G_random (networkx.MultiDiGraph): A new graph with shuffled edges.
#     """
#     G_random = nx.MultiDiGraph()
#     G_random.update(G_prime)

#     edges = list(G_random.edges(keys=True, data=True))  # List of all edges with keys
#     edge_dict = {(u, v, k): d["label"] for u, v, k, d in edges}  # Store labels

#     for _ in range(num_swaps):
#         # Randomly pick two different edges
#         (u, v, k1, data1), (x, y, k2, data2) = random.sample(edges, 2)

#         label_uv = data1["label"]
#         label_xy = data2["label"]

#         # Ensure we can swap
#         if (u, y, k1) in edge_dict and (x, v, k2) in edge_dict:
#             label_uy = edge_dict[(u, y, k1)]
#             label_xv = edge_dict[(x, v, k2)]

#             if label_uv == label_xy and label_uy == label_xv:
#                 # Perform label swap in the actual graph
#                 G_random[u][v][k1]["label"], G_random[u][y][k1]["label"] = label_uy, label_uv
#                 G_random[x][y][k2]["label"], G_random[x][v][k2]["label"] = label_xv, label_xy

#                 # Update the edge dictionary
#                 edge_dict[(u, v, k1)], edge_dict[(u, y, k1)] = label_uy, label_uv
#                 edge_dict[(x, y, k2)], edge_dict[(x, v, k2)] = label_xv, label_xy

#                 print(f"Swapped edges: ({u}, {v}) <-> ({u}, {y}) and ({x}, {y}) <-> ({x}, {v})")

#     return G_random

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

    num_swaps = 5

    for txid in taxon_ids:
        ppi_path = f"data/oxidative_stress/{txid}/stress_ppi.csv"
        reg_path = f"data/oxidative_stress/{txid}/stress_reg.csv"
        # output_dir = f"data/oxidative_stress/{txid}/randomized_networks/"
        # out_ppi_path = f"{output_dir}stress_ppi{num_swaps}.csv"
        # out_reg_path = f"{output_dir}stress_reg{num_swaps}.csv"

        G = read_csv(
            ppi_path,
            reg_path
        )

        print(f"Original graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

        # Example: Relabel edges
        G_prime = label_edges(G)

        # Compare expected number of edges vs actual
        unique_node_pairs = set()
        for u, v in G.edges():
            unique_node_pairs.add(tuple(sorted([u, v])))  # Sorting ensures (A, B) == (B, A)

        expected_labeled_edges = len(unique_node_pairs)
        actual_labeled_edges = len(G_prime.edges())
        if actual_labeled_edges == expected_labeled_edges:
            print("Edge count matches expectation!")
        else:
            print("Edge count mismatch! Check relabeling logic.")

        G_random = swap_edges(G_prime, num_swaps)
        
        # # Example: Print basic graph stats
        # print(f"Original graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
        
        print(f"Labeled subnetwork: {G_prime.number_of_nodes()} nodes, {G_prime.number_of_edges()} edges")
        print(f"Shuffled subnetwork: {G_random.number_of_nodes()} nodes, {G_random.number_of_edges()} edges")

        # Compare edge label distributions
        original_label_counts = {label: 0 for label in set(nx.get_edge_attributes(G_prime, "label").values())}
        shuffled_label_counts = {label: 0 for label in set(nx.get_edge_attributes(G_prime, "label").values())}

        for _, _, d in G_prime.edges(data=True):
            original_label_counts[d["label"]] += 1

        for _, _, d in G_random.edges(data=True):
            shuffled_label_counts[d["label"]] += 1

        print("Original label counts:", original_label_counts)
        print("Shuffled label counts:", shuffled_label_counts)

        # # Check output graphs
        # print("First 5 edges in G_prime:")
        # for edge in list(G_prime.edges(data=True))[:5]:
        #     print(edge)

        # print("\nFirst 5 edges in G_random:")
        # for edge in list(G_random.edges(data=True))[:5]:
        #     print(edge)

        # Check if edges have changed after shuffling
        edges_prime = set((u, v, tuple(sorted(d.items()))) for u, v, d in G_prime.edges(data=True))
        edges_random = set((u, v, tuple(sorted(d.items()))) for u, v, d in G_random.edges(data=True))

        if edges_prime == edges_random:
            print("No changes in edges! Shuffling may not be working.")
        else:
            print("Edges have changed after shuffling.")

        only_in_prime = edges_prime - edges_random
        only_in_random = edges_random - edges_prime

        print(f"Edges unique to G_prime: {len(only_in_prime)}")
        print(f"Edges unique to G_random: {len(only_in_random)}")

        # Print a few examples of swapped edges
        if only_in_prime and only_in_random:
            print("Example of changed edges:")
            print("Before shuffle:", list(only_in_prime)[:5])
            print("After shuffle:", list(only_in_random)[:5])

        # Example: Draw the graphs
        # nx.draw_networkx(G_prime, with_labels=True, font_size=10)
        # plt.show()

        # nx.draw_networkx(G_random, with_labels=True, font_size=10)
        # plt.show()

        # split_to_csv(G_random, out_ppi_path, out_reg_path)

if __name__ == "__main__":
    main()