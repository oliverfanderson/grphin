import networkx as nx
import random
import csv
import os
from matplotlib import pyplot as plt
import argparse
from collections import defaultdict
import copy

parser = argparse.ArgumentParser()
parser.add_argument(
    "-i", "--iterations", required=True, help="Number of iterations to run"
)
parser.add_argument(
    "-s", "--swaps", required=True, help="Number of edges to swap per iteration"
)
args = parser.parse_args()

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
    """
    Helper function to process edges and add them to the graph.
    """

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
    print(f"Removed self-loops: {len(self_loops)} edges")
    # print(self_loops)

    return G

def label_edges(G):
    """
    Assigns new edge labels based on the combination of PPI and Reg edges between node pairs.
    Ensures that each node pair is processed only once.

    Parameters:
        G (networkx.MultiDiGraph): Input graph with edges labeled as 'ppi' or 'reg'.

    Returns:
        G_prime (networkx.DiGraph): A new graph with relabeled edges based on 2-node graphlets: 'ppi', 'reg', 'mix', 'coreg', 'coreg_ppi'.
    """

    G_prime = nx.DiGraph()
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
        elif num_ppi > 0 and (has_reg_uv and not has_reg_vu) and not has_coreg:
            new_label = "mix"  # One Reg + One PPI (u → v)
        elif num_ppi > 0 and (has_reg_vu and not has_reg_uv) and not has_coreg:
            new_label = "mix"  # One Reg + One PPI (v → u)
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
    """
    Performs constrained edge swaps in a MultiDiGraph while preserving connectivity by ensuring that swapped edges maintain the same label.
    
    Parameters:
        G_prime (nx.DiGraph): The input graph to randomize with 2-node graphlet edge types labeled.
        num_swaps (int): The number of swaps to attempt.
        
    Returns:
        G_random (nx.DiGraph): A randomized version of G_prime.
    """

    G_random = nx.DiGraph()
    G_random.update(G_prime)
    edges = list(G_random.edges(data=True))  # (u, v, data)
    swaps = 0
    # random.seed(42)
    
    while swaps < num_swaps:
        # Select two random edges (ensuring distinct nodes)
        (u, v, data1), (x, y, data2) = random.sample(edges, 2)

        if len({u, v, x, y}) < 4:
            continue  # Skip if nodes are not unique
        
        # Ensure the edges have the same label
        uv_type = data1.get("label")
        xy_type = data2.get("label")
        if uv_type != xy_type or uv_type is None or xy_type is None:
            continue

        # Ensure the edges have the same label
        uy_type = G_random[u][y]["label"] if G_random.has_edge(u, y) else None
        xv_type = G_random[x][v]["label"] if G_random.has_edge(x, v) else None
        if uy_type != xv_type:
            continue
        # print(f"Swapping edges: {u}->{v} and {x}->{y}")
        # print(f"UV: {uv_type}, XY: {xy_type}, UY: {uy_type}, XV: {xv_type}")

        uv_edge = G_random[u][v]
        xy_edge = G_random[x][y]
        if G_random.has_edge(v, u):
            vu_edge = G_random[v][u]
        else:
            vu_edge = None
        if G_random.has_edge(y, x):
            yx_edge = G_random[y][x]
        else:
            yx_edge = None
        if G_random.has_edge(u, y):
            uy_edge = G_random[u][y]
        else:
            uy_edge = None
        if G_random.has_edge(y, u):
            yu_edge = G_random[y][u]
        else:
            yu_edge = None
        if G_random.has_edge(x, v):
            xv_edge = G_random[x][v]
        else:
            xv_edge = None
        if G_random.has_edge(v, x):
            vx_edge = G_random[v][x]
        else:
            vx_edge = None
        # print(f"Swapping ({u}, {v}) -> ({u}, {y}) and ({x}, {y}) -> ({x}, {v})")
        # print(f"UV: {uv_edge}, XY: {xy_edge}, UY: {uy_edge}, XV: {xv_edge}")

        ## Condition 1: No u,y and x,v edges
        # Perform the swap: (u, v) ↔ (u, y) and (x, y) ↔ (x, v)
        if vu_edge is None and yx_edge is None and yu_edge is None and vx_edge is None:
            if uy_edge is None and xv_edge is None:
                G_random.remove_edge(u, v)
                G_random.remove_edge(x, y)
                G_random.add_edge(u, y, **uv_edge)
                G_random.add_edge(x, v, **xy_edge)
            else:
                # IDEA: Maybe we can just swap label instead of removing and adding new edges
                G_random.remove_edge(u, v)
                G_random.remove_edge(x, y)
                G_random.remove_edge(u, y)
                G_random.remove_edge(x, v)
                G_random.add_edge(u, y, **uv_edge)
                G_random.add_edge(x, v, **xy_edge)
                G_random.add_edge(u, v, **uy_edge)
                G_random.add_edge(x, y, **xv_edge)
        
        # Update the edges list
        edges = list(G_random.edges(data=True))
        swaps += 1
    
    return G_random

# def swap_edges(G_prime, num_swaps):
#     G_random = nx.DiGraph()
#     G_random.update(G_prime)
#     edges = list(G_random.edges(data=True))  # (u, v, data)
#     swaps = 0
#     random.seed(42)
    
#     while swaps < num_swaps:
#         # Select two random edges (ensuring distinct nodes)
#         (u, v, data1), (x, y, data2) = random.sample(edges, 2)

#         if len({u, v, x, y}) < 4:
#             continue  # Skip if nodes are not unique
        
#         # Ensure the edges have the same label
#         uv_type = data1.get("label")
#         xy_type = data2.get("label")
#         if uv_type != xy_type or not uv_type:
#             continue

#         # Ensure the edges have the same label for the swap
#         uy_type = G_random.get_edge_data(u, y, default={}).get("label")
#         xv_type = G_random.get_edge_data(x, v, default={}).get("label")
#         if uy_type != xv_type:
#             continue

#         # Get the edges and their data
#         uv_edge = data1
#         xy_edge = data2
#         uy_edge = G_random.get_edge_data(u, y, default={})
#         xv_edge = G_random.get_edge_data(x, v, default={})

#         # Perform the swap: (u, v) ↔ (u, y) and (x, y) ↔ (x, v)
#         if not uy_edge and not xv_edge:
#             G_random.remove_edge(u, v)
#             G_random.remove_edge(x, y)
#             G_random.add_edge(u, y, **uv_edge)
#             G_random.add_edge(x, v, **xy_edge)
#         else:
#             # Swap the edges
#             G_random.remove_edge(u, v)
#             G_random.remove_edge(x, y)
#             G_random.remove_edge(u, y)
#             G_random.remove_edge(x, v)
#             G_random.add_edge(u, y, **uv_edge)
#             G_random.add_edge(x, v, **xy_edge)
#             if uy_edge:
#                 G_random.add_edge(u, v, **uy_edge)
#             if xv_edge:
#                 G_random.add_edge(x, y, **xv_edge)
        
#         # Update the edges list
#         edges = list(G_random.edges(data=True))
#         swaps += 1
    
#     return G_random

def split_to_csv(G_random, out_ppi_path, out_reg_path):
    """
    Writes the randomized graph to CSV files based on 2-node graphlet edge labels.

    Parameters:
        G_random (networkx.DiGraph): The randomized graph.
        out_ppi_path (string): A filepath to write the set of PPI edges (tuples of (u, v)).
        out_reg_path (string): A filepath to write the set of Reg edges (tuples of (u, v)).

    Returns:
        out_ppi_path (CSV): The randomized PPI edges CSV file.
        out_reg_path (CSV): The randomized Reg edges CSV file.
    """

    # # Keep track of written edges
    # ppi_written = set()  # Set to track edges written to the PPI file
    # reg_written = set()  # Set to track edges written to the Reg file

    # Write edges to CSV files
    with open(out_ppi_path, "w", newline="") as ppi_out, open(out_reg_path, "w", newline="") as reg_out:
        ppi_writer = csv.writer(ppi_out, quotechar='"', quoting=csv.QUOTE_ALL)
        reg_writer = csv.writer(reg_out, quotechar='"', quoting=csv.QUOTE_ALL)
        
        # Write CSV headers
        ppi_writer.writerow(["id1", "id2"])
        reg_writer.writerow(["id1", "id2"])
        
        # Iterate over edges
        for u, v, data in G_random.edges(data=True):
            label = data.get("label", None)
            if label == "ppi":
                # print(f"Label: PPI. Adding PPI edge: {u} -> {v}")
                ppi_writer.writerow([u, v])
                # print(f"Label: PPI. Adding PPI edge: {v} -> {u}")
                ppi_writer.writerow([v, u])
            elif label == "reg":
                # print(f"Label: Reg. Adding Reg edge: {u} -> {v}")
                reg_writer.writerow([u, v])
            elif label == "mix":
                # print(f"Label: Mix. Adding PPI edge: {u} -> {v}")
                ppi_writer.writerow([u, v])
                # print(f"Label: Mix. Adding PPI edge: {v} -> {u}")
                ppi_writer.writerow([v, u])
                # print(f"Label: Mix. Adding Reg edge: {u} -> {v}")
                reg_writer.writerow([u, v])
            elif label == "coreg":
                # print(f"Label: Coreg. Adding Reg edge: {u} -> {v}")
                reg_writer.writerow([u, v])
                # print(f"Label: Coreg. Adding Reg edge: {v} -> {u}")
                reg_writer.writerow([v, u])
            elif label == "coreg_ppi":
                # print(f"Label: Coreg_PPI. Adding PPI edge: {u} -> {v}")
                ppi_writer.writerow([u, v])
                # print(f"Label: Coreg_PPI. Adding PPI edge: {v} -> {u}")
                ppi_writer.writerow([v, u])
                # print(f"Label: Coreg_PPI. Adding Reg edge: {u} -> {v}")
                reg_writer.writerow([u, v])
                # print(f"Label: Coreg_PPI. Adding Reg edge: {v} -> {u}")
                reg_writer.writerow([v, u])

    print(f"PPI edges written to: {out_ppi_path}")
    print(f"Reg edges written to: {out_reg_path}")


def main():
    """
    A function to generate randomized networks based on 2-node graphlet edge labels.
    
    Parameters:
        -s / --swaps: Command-line argument for number of swaps.
        -i / --iterations: Command-line argument for number of iterations.

    Returns:
        Randomized PPI and Reg interaction CSV files for each taxon ID with (-s) swaps performed on each (-i) iteration.

    Example:
        python3 enrichment.py --swaps 1000 --iterations 10

        This will generate 10 randomized networks with 1000 edge swaps for each taxon ID (txid6239, txid7227, txid7955, txid224308, txid559292).
    """

   # List of taxon IDs to process
    # taxon_ids = ["txid6239", "txid7227", "txid7955", "txid224308", "txid559292"]
    taxon_ids = ["txid6239"]

    num_swaps = int(args.swaps)
    num_iterations = int(args.iterations)
    for iteration in range(num_iterations):
        for txid in taxon_ids:
            ppi_path = f"data/oxidative_stress/{txid}/stress_ppi.csv"
            reg_path = f"data/oxidative_stress/{txid}/stress_reg.csv"
            output_dir = f"data/oxidative_stress/{txid}/randomized_networks"
            out_ppi_path = f"{output_dir}/stress_ppi{iteration}.csv"
            out_reg_path = f"{output_dir}/stress_reg{iteration}.csv"

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
            # Dictionary to store edges and their labels
            edge_labels = defaultdict(set)

            # Populate the dictionary with edge labels
            for u, v, data in G_random.edges(data=True):
                label = data.get("label", None)
                edge_labels[tuple(sorted((u, v)))].add(label)

            # Print edges that have more than one label
            for edge, labels in edge_labels.items():
                if len(labels) > 1:
                    print(f"Edge {edge} has multiple labels: {labels}")
                    break

            # Check nodes and edges match after randomization
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
            shuffled_label_counts = {label: 0 for label in set(nx.get_edge_attributes(G_random, "label").values())}

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