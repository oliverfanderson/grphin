import networkx as nx
import random
import csv
import os
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
    # print(f"Removed self-loops: {len(self_loops)} edges")
    # print(self_loops)

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

    # Step 1: Track PPI and directed Reg edges separately
    edge_info = {}
    for u, v, d in G.edges(data=True):
        edge_label = d["label"]
        key = tuple(sorted((u, v)))  # Ensure undirected edges are processed consistently

        if key not in edge_info:
            edge_info[key] = {"ppi": 0, "reg_uv": False, "reg_vu": False}

        if edge_label == "ppi":
            edge_info[key]["ppi"] += 1  # Count PPI edges
        elif edge_label == "reg":
            if (u, v) == key:  # Check if the edge follows the key's order
                edge_info[key]["reg_uv"] = True  # Mark directed reg edge u → v
            else:
                edge_info[key]["reg_vu"] = True  # Mark directed reg edge v → u

    # Step 2: Assign new labels based on rules
    processed_pairs = set()
    for (u, v), counts in edge_info.items():
        if (u, v) in processed_pairs:
            continue  # Skip if already processed

        num_ppi = counts["ppi"]
        has_reg_uv = counts["reg_uv"]
        has_reg_vu = counts["reg_vu"]
        has_coreg = has_reg_uv and has_reg_vu  # Check reciprocal regulation

        # Determine new edge label
        if num_ppi > 0 and not has_reg_uv and not has_reg_vu:
            new_label = "ppi"  # Only PPI
        elif num_ppi == 0 and has_reg_uv and not has_reg_vu:
            new_label = "reg"  # Only one Reg (u → v)
        elif num_ppi == 0 and has_reg_vu and not has_reg_uv:
            new_label = "reg"  # Only one Reg (v → u)
        elif num_ppi > 0 and (has_reg_uv or has_reg_vu) and not has_coreg:
            new_label = "mix"  # One Reg + One PPI
        elif has_coreg and num_ppi == 0:
            new_label = "coreg"  # Reciprocal regulation (u → v and v → u)
        elif has_coreg and num_ppi > 0:
            new_label = "coreg_ppi"  # Reciprocal regulation + PPI
        else:
            continue  # Shouldn't happen

        # Add the relabeled edge once
        G_prime.add_edge(u, v, label=new_label)
        processed_pairs.add((u, v))  # Mark as processed

    return G_prime

def swap_edges(G_prime, num_swaps):
    """Performs constrained edge swaps in a MultiDiGraph while preserving connectivity.
    
    Parameters:
        G (nx.MultiDiGraph): The input graph to shuffle.
        num_swaps (int): The number of swaps to attempt.
        
    Returns:
        nx.MultiDiGraph: A shuffled version of G.
    """
    G_random = nx.MultiDiGraph()
    G_random.update(G_prime)
    edges = list(G_random.edges(keys=True, data=True))  # (u, v, key, data)
    swaps = 0
    
    while swaps < num_swaps:
        # Select two random edges (ensuring distinct nodes)
        (u, v, key1, data1), (x, y, key2, data2) = random.sample(edges, 2)

        if len({u, v, x, y}) < 4:
            continue  # Skip if nodes are not unique
        
        # Ensure the edges have the same label
        uv_type = data1.get("label")
        xy_type = data2.get("label")
        if uv_type != xy_type or uv_type is None:
            continue

        # Ensure the edges have the same label
        uy_type = G_random[u][y][key1]["label"] if G_random.has_edge(u, y, key1) else None
        xv_type = G_random[x][v][key2]["label"] if G_random.has_edge(x, v, key2) else None
        if uy_type != xv_type:
            continue
        # print(f"DEBUG: Attempting swap {swaps + 1}:")
        # print(f"DEBUG: U: {u}, V: {v}, Key1: {key1}, Data1: {data1}")
        # print(f"DEBUG: X: {x}, Y: {y}, Key2: {key2}, Data2: {data2}")
        # print(f"UV: {uv_type}, XY: {xy_type}, UY: {uy_type}, XV: {xv_type}")

        uv_edge = G_random[u][v][key1]
        xy_edge = G_random[x][y][key2]
        if G_random.has_edge(u, y, key1):
            uy_edge = G_random[u][y][key1]
        else:
            uy_edge = None
        if G_random.has_edge(x, v, key2):
            xv_edge = G_random[x][v][key2]
        else:
            xv_edge = None
        # print(f"DEBUG: Swapping edges: {uv_edge}, {xy_edge}, {uy_edge}, {xv_edge}")
        
        # Perform the swap: (u, v) ↔ (u, y) and (x, y) ↔ (x, v)
        if uy_edge is None:
            G_random.remove_edge(u, v, key1)
            G_random.remove_edge(x, y, key2)
            G_random.add_edge(u, y, key=key1, **uv_edge)
            G_random.add_edge(x, v, key=key2, **xy_edge)
        else:
            # Look if we can just update labels
            G_random.remove_edge(u, v, key1)
            G_random.remove_edge(x, y, key2)
            G_random.remove_edge(u, y, key1)
            G_random.remove_edge(x, v, key2)
            G_random.add_edge(u, y, key=key1, **uv_edge)
            G_random.add_edge(x, v, key=key2, **xy_edge)
            G_random.add_edge(u, v, key=key1, **uy_edge)
            G_random.add_edge(x, y, key=key2, **xv_edge)
            # print(f"DEBUG: Swapped edges: {uv_edge}, {xy_edge}, {uy_edge}, {xv_edge}")
        
        # Update the edges list
        edges = list(G_random.edges(keys=True, data=True))
        swaps += 1
    
    return G_random

def split_to_csv(G_random, out_ppi_path, out_reg_path):
    """
    Writes the graph to CSV files based on edge labels.

    Parameters:
    G_random (networkx.Graph): The randomized graph.
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
        for u, v, key, data in G_random.edges(data=True, keys=True):
            # print(f"Edge ({u}, {v}, {key}): {data}")
            label = data.get("label", None)
            if label == "ppi":
                ppi_writer.writerow([u, v])
                ppi_writer.writerow([v, u])
            elif label == "reg":
                reg_writer.writerow([u, v])
            elif label == "mix":
                ppi_writer.writerow([u, v])
                ppi_writer.writerow([v, u])
                reg_writer.writerow([u, v])
            elif label == "coreg":
                reg_writer.writerow([u, v])
                reg_writer.writerow([v, u])
            elif label == "coreg_ppi":
                ppi_writer.writerow([u, v])
                ppi_writer.writerow([v, u])
                reg_writer.writerow([u, v])
                reg_writer.writerow([v, u])
 

    print(f"PPI edges written to: {out_ppi_path}")
    print(f"Reg edges written to: {out_reg_path}")


def main():
   # List of taxon IDs to process
    # taxon_ids = ["txid6239", "txid7227", "txid7955", "txid224308", "txid559292"]
    taxon_ids = ["txid559292"]

    num_swaps = 500

    for txid in taxon_ids:
        ppi_path = f"data/oxidative_stress/{txid}/stress_ppi.csv"
        reg_path = f"data/oxidative_stress/{txid}/stress_reg.csv"
        output_dir = f"data/oxidative_stress/{txid}/randomized_networks"
        out_ppi_path = f"{output_dir}/stress_ppi{num_swaps}.csv"
        out_reg_path = f"{output_dir}/stress_reg{num_swaps}.csv"

        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

        G = read_csv(
            ppi_path,
            reg_path
        )

        print(f"Original graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

        # Relabel edges with 5 two-node graphlet types
        G_prime = label_edges(G)

        # Compare expected number of edges vs actual in G_prime
        unique_node_pairs = set()
        for u, v in G.edges():
            unique_node_pairs.add(tuple(sorted([u, v])))  # Sorting ensures (A, B) == (B, A)

        # Break loop if expected number of labeled edges does not match actual in G_prime
        expected_labeled_edges = len(unique_node_pairs)
        actual_labeled_edges = len(G_prime.edges())
        if actual_labeled_edges != expected_labeled_edges:
            print(f"Expected {expected_labeled_edges} labeled edges, but found {actual_labeled_edges}!")
            break

        # Randomize the graph
        G_random = swap_edges(G_prime, num_swaps)
        
        # Validate the random graph
        original_nodes = G_prime.number_of_nodes()
        randomized_nodes = G_random.number_of_nodes()
        if original_nodes != randomized_nodes:
            print("Number of nodes does not match after shuffling!")
            break

        original_edges = G_prime.number_of_edges()
        randomized_edges = G_random.number_of_edges()
        if original_edges != randomized_edges:
            print("Number of edges does not match after shuffling!")
            print(f"Original edges: {original_edges}, Randomized edges: {randomized_edges}!")
            break

        # Compare edge label distributions
        original_label_counts = {label: 0 for label in set(nx.get_edge_attributes(G_prime, "label").values())}
        shuffled_label_counts = {label: 0 for label in set(nx.get_edge_attributes(G_prime, "label").values())}

        for _, _, d in G_prime.edges(data=True):
            original_label_counts[d["label"]] += 1

        for _, _, d in G_random.edges(data=True):
            shuffled_label_counts[d["label"]] += 1

        print("Original label counts:", original_label_counts)
        print("Shuffled label counts:", shuffled_label_counts)

        if not all(original_label_counts[label] == shuffled_label_counts[label] for label in original_label_counts):
            print("Edge label distributions do not match after shuffling!")
            break
        
        # Check degree sequence
        original_degree = dict(G_prime.degree())
        random_degree = dict(G_random.degree())
        if original_degree != random_degree:
            print("Degree sequence does not match after shuffling!")
            break

        # Check if edges have changed after shuffling
        edges_prime = set((u, v, tuple(sorted(d.items()))) for u, v, d in G_prime.edges(data=True))
        edges_random = set((u, v, tuple(sorted(d.items()))) for u, v, d in G_random.edges(data=True))

        # Get unique edges in G_prime and G_random
        only_in_prime = edges_prime - edges_random
        only_in_random = edges_random - edges_prime

        if edges_prime == edges_random:
            print("No changes in edges! Shuffling may not be working.")
        else:
            print(f"{len(only_in_prime)} edges have changed after shuffling.")

        if len(only_in_prime) != len(only_in_random):
            print("Mismatch between swapped edges in graphs.")
            print(f"Edges unique to G_prime: {len(only_in_prime)}")
            print(f"Edges unique to G_random: {len(only_in_random)}")
            break
        

        # Print a few examples of swapped edges
        # if only_in_prime and only_in_random:
        #     print("Before shuffle:", list(only_in_prime)[:5])
        #     print("After shuffle:", list(only_in_random)[:5])

        # Troubleshooting: Draw the graphs
        # nx.draw_networkx(G_prime, with_labels=True, font_size=10)
        # plt.show()

        # nx.draw_networkx(G_random, with_labels=True, font_size=10)
        # plt.show()

        split_to_csv(G_random, out_ppi_path, out_reg_path)

if __name__ == "__main__":
    main()