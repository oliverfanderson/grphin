import argparse
import ast
from collections import defaultdict
from functools import lru_cache
import re
import sys
import time
from matplotlib import pyplot as plt
import networkx as nx
import csv
from pathlib import Path
import numpy as np
from collections import defaultdict


def count_two_node_graphlet(G, output_dir):
    """
    wrapper function that initializes variables, counts graphlet and orbits in the network, and generates output files for two node graphlets and orbits
    Args:
        G (nx.MultiDiGraph) : stores the multi directed graph representation of our interactome.
        output_dir (str) : string path of output folder
    Returns:
        two_node_graphlet_count (dict)
        two_node_graphlet_id (dict)
        two_node_orbit_dict (dict)
    """
    two_node_graphlet_count, two_node_graphlet_id = initialize_two_node_graphlet_data()

    two_node_graphlet_count, two_node_orbit_dict = get_two_node_graphlets_and_orbits(
        G, two_node_graphlet_count
    )

    generate_two_node_output_files(
        two_node_graphlet_id,
        two_node_graphlet_count,
        two_node_orbit_dict,
        output_dir,
    )

    return two_node_graphlet_count, two_node_graphlet_id, two_node_orbit_dict


def generate_two_node_output_files(
    two_node_graphlet_id,
    two_node_graphlet_count,
    two_node_orbit_dict,
    output_dir,
):
    # print("\nTwo node graphlet counts")
    # for key in two_node_graphlet_id:
    #     print(f"G_{two_node_graphlet_id[key]} = {two_node_graphlet_count[key]}")

    # print("\nTwo node graphlet orbit counts")
    # for key in two_node_orbit_dict:
    #     print(f"{key} = {len(two_node_orbit_dict[key])}")

    with open(f"{output_dir}/two_node_graphlet_counts.csv", "w") as f:
        for key in two_node_graphlet_id:
            f.write(f"G_{two_node_graphlet_id[key]}, {two_node_graphlet_count[key]}\n")
    f.close()

    with open(f"{output_dir}/two_node_orbit_counts.csv", "w") as f:
        for key in two_node_orbit_dict:
            f.write(f"{key}, {len(two_node_orbit_dict[key])}\n")
        f.close()

    return None


def initialize_graphlet_data(network_ppi_path, network_reg_path):
    protein_id = get_protein_id_dict(network_ppi_path, network_reg_path)
    G = read_csv(
        network_ppi_path,
        network_reg_path,
        protein_id,
        node_size_limit=999999999,
        edge_size_limit=999999999,
    )
    G_prime = simplify_graph_to_undirected(G)
    graphlet_config = load_graphlet_config("graphlet_config.csv")

    return protein_id, G, G_prime, graphlet_config


def read_csv(
    ppi_path,
    reg_path,
    protein_id_dict,
    node_size_limit=float("inf"),
    edge_size_limit=float("inf"),
):
    """Reads CSV files and constructs a graph with edges labeled as 'ppi' or 'reg'."""
    G = nx.MultiDiGraph()
    visited_nodes = set()

    process_edges(
        ppi_path,
        G,
        protein_id_dict,
        visited_nodes,
        "ppi",
        node_size_limit,
        edge_size_limit,
    )
    process_edges(
        reg_path,
        G,
        protein_id_dict,
        visited_nodes,
        "reg",
        node_size_limit,
        edge_size_limit,
    )
    print()

    return G


def process_edges(
    file_path, G, protein_id_dict, visited_nodes, label, node_limit, edge_limit
):
    """Helper function to process edges and add them to the graph."""
    # print(f"Processing {label} edges")
    with open(file_path, "r") as file:
        csv_reader = csv.reader(file)
        next(csv_reader)  # Skip header
        node_count = len(visited_nodes)
        edge_count = G.number_of_edges()

        for row in csv_reader:
            if node_count > node_limit or edge_count > edge_limit:
                break

            id1 = protein_id_dict[row[0]]
            id2 = protein_id_dict[row[1]]

            if id1 not in visited_nodes:
                visited_nodes.add(id1)
                node_count += 1

            if id2 not in visited_nodes:
                visited_nodes.add(id2)
                node_count += 1

            G.add_edge(id1, id2, label=label)
            edge_count += 1


def initialize_two_node_graphlet_data():
    """
    Get all 2-node graphlet binary edge vectors
    see diagram for reference
    """

    # 2-node graphlet binary edge vectors
    # See 2-node graphlet figure in GRPhIN paper
    a_1 = (1, 0, 0)
    a_2 = (1, 0, 0)
    a_hash = hash(a_1 + a_2)

    b_1 = (0, 1, 0)
    b_2 = (0, 0, 1)
    b_hash = hash(b_1 + b_2)

    c_1 = (1, 1, 0)
    c_2 = (1, 0, 1)
    c_hash = hash(c_1 + c_2)

    d_1 = (1, 1, 1)
    d_2 = (1, 1, 1)
    d_hash = hash(d_1 + d_2)

    e_1 = (0, 1, 1)
    e_2 = (0, 1, 1)
    e_hash = hash(e_1 + e_2)

    two_node_graphlet_counts = {a_hash: 0, b_hash: 0, c_hash: 0, d_hash: 0, e_hash: 0}
    two_node_graphlet_id = {a_hash: 1, b_hash: 2, c_hash: 3, d_hash: 4, e_hash: 5}
    return two_node_graphlet_counts, two_node_graphlet_id


def get_protein_id_dict(ppi_path, reg_path):
    """Creates a dictionary mapping protein IDs to unique integers."""

    res_dict = {}
    start_index = 0
    start_index = update_protein_id_dict(ppi_path, res_dict, start_index)
    update_protein_id_dict(reg_path, res_dict, start_index)
    return res_dict


def update_protein_id_dict(file_path, res_dict, start_index):
    """Helper function to update the protein ID dictionary from a file."""

    with open(file_path, "r") as file:
        csv_reader = csv.reader(file)
        next(csv_reader)
        for row in csv_reader:
            for protein_id in row[:2]:
                res_dict.setdefault(protein_id, start_index)
                if res_dict[protein_id] == start_index:
                    start_index += 1
    return start_index


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



def load_graphlet_config(file_path):
    """Load graphlet lookup table from a CSV file. Creates a list."""

    graphlet_config = []
    with open(file_path, mode="r") as file:
        reader = csv.DictReader(file)
        for row in reader:
            graphlet_config.append(
                {
                    "key": ast.literal_eval(row["key"]),
                    "a_expected": ast.literal_eval(row["a_expected"]),
                    "b_expected": ast.literal_eval(row["b_expected"]),
                    "c_expected": ast.literal_eval(row["c_expected"]),
                    "orbits": (
                        int(row["orbit1"]),
                        int(row["orbit2"]),
                        int(row.get("orbit3", -1)),
                    ),  # Handle missing orbit3
                }
            )
    return graphlet_config


def load_orbit_config(file_path):
    """Load orbit lookup table from a CSV file. Creates a dictionary."""

    graphlet_config = {}
    with open(file_path, mode="r") as file:
        reader = csv.DictReader(file)
        for row in reader:
            graphlet_config[ast.literal_eval(row["key"])] = {
                "a_expected": ast.literal_eval(row["a_expected"]),
                "b_expected": ast.literal_eval(row["b_expected"]),
                "c_expected": ast.literal_eval(row["c_expected"]),
                "orbits": (
                    int(row["orbit1"]),
                    int(row["orbit2"]),
                    int(row.get("orbit3", -1)),
                ),  # Handle missing orbit3}
            }
    return graphlet_config


def get_two_node_graphlets_and_orbits(G, two_node_graphlet_dict):
    """
    Get the counts of 2-node graphlets in the graph and their orbit positions.

    Parameters:
        G (nx.MultiDiGraph): Stores the MultiDiGraph representation of the mixed interaction network.
        two_node_graphlet_dict (dict)
    Returns:
        two_node_graphlet_dict (dict)
        two_node_orbit_dict (dict)
    """
    G_adj_list = get_two_node_adjacency_list(G)
    two_node_orbit_dict = {1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: []}
    for neighbors in G_adj_list:
        i = neighbors[0]
        for j in neighbors[1]:
            vectors = G_adj_list[i][1][j] + G_adj_list[j][1][i]
            two_node_orbit_dict = get_two_node_orbit_position(
                i, j, G_adj_list[i][1][j], G_adj_list[j][1][i], two_node_orbit_dict
            )
            if hash(tuple(vectors)) in two_node_graphlet_dict:
                two_node_graphlet_dict[hash(tuple(vectors))] += 1
    return two_node_graphlet_dict, two_node_orbit_dict


def get_two_node_adjacency_list(G):
    """Get the adjacency list for a MultiDiGraph"""

    print("\nGetting 2-node adjacency list...")

    adj_list_vector = [{} for _ in range(len(G.nodes()))]

    for i, j, data in G.edges(data=True):
        label = data.get("label")
        if label == "ppi":
            if j not in adj_list_vector[i]:
                adj_list_vector[i][j] = [1, 0, 0]
            if i not in adj_list_vector[j]:
                adj_list_vector[j][i] = [1, 0, 0]
        elif label == "reg":
            if j not in adj_list_vector[i]:
                adj_list_vector[i][j] = [0, 1, 0]
            else:
                adj_list_vector[i][j][1] += 1
            if i not in adj_list_vector[j]:
                adj_list_vector[j][i] = [0, 0, 1]
            else:
                adj_list_vector[j][i][2] += 1

    final_adj_list_vector = [
        (i, neighbors) for i, neighbors in enumerate(adj_list_vector)
    ]
    return final_adj_list_vector


def get_two_node_orbit_position(a, b, a_vec, b_vec, orbit_dict):
    """
    check orbit positioning via vectors
    """
    if a_vec == [1, 0, 0] and b_vec == [1, 0, 0]:
        if a not in orbit_dict[1]:
            orbit_dict[1] += [a]
        if b not in orbit_dict[1]:
            orbit_dict[1] += [b]
    elif (
        a_vec == [0, 1, 0]
        and b_vec == [0, 0, 1]
        or a_vec == [0, 0, 1]
        and b_vec == [0, 1, 0]
    ):
        if a_vec == [0, 1, 0] and a not in orbit_dict[2]:
            orbit_dict[2] += [a]
            orbit_dict[3] += [b]
        elif b not in orbit_dict[2]:
            orbit_dict[2] += [b]
            orbit_dict[3] += [a]

    elif (
        a_vec == [1, 1, 0]
        and b_vec == [1, 0, 1]
        or a_vec == [1, 0, 1]
        and b_vec == [1, 1, 0]
    ):
        if a_vec == [1, 1, 0] and a not in orbit_dict[4]:
            orbit_dict[4] += [a]
            orbit_dict[5] += [b]
        elif b not in orbit_dict[4]:
            orbit_dict[4] += [b]
            orbit_dict[5] += [a]

    elif a_vec == [1, 1, 1] and b_vec == [1, 1, 1]:
        if a not in orbit_dict[6]:
            orbit_dict[6] += [a]
        if b not in orbit_dict[6]:
            orbit_dict[6] += [b]

    elif a_vec == [0, 1, 1] and b_vec == [0, 1, 1]:
        if a not in orbit_dict[7]:
            orbit_dict[7] += [a]
        if b not in orbit_dict[6]:
            orbit_dict[7] += [b]

    return orbit_dict


def hash_tuple(xy):
    """Hash a tuple helper function."""

    return hash(tuple(xy))


# Precompute the direct mapping from tuple to index
tuple_to_index = {
    (0, 0, 0): 0,
    (1, 0, 0): 1,
    (0, 1, 0): 2,
    (0, 0, 1): 3,
    (1, 1, 0): 4,
    (1, 0, 1): 5,
    (0, 1, 1): 6,
    (1, 1, 1): 7,
}

default_edge = [0, 0, 0]


def get_three_node_graphlet_dict(ab, ac, ba, bc, ca, cb):
    """
    Direct lookup with precomputed tuple-to-index mapping
    """
    return (
        tuple_to_index[tuple(ab)],
        tuple_to_index[tuple(ac)],
        tuple_to_index[tuple(ba)],
        tuple_to_index[tuple(bc)],
        tuple_to_index[tuple(ca)],
        tuple_to_index[tuple(cb)],
    )


def grphin_algorithm(
    G: nx.MultiDiGraph,
    G_prime: nx.Graph,
    three_node_graphlet_count,
    three_node_orbit_protein_data,
):
    """
    GRPhIN algorithm to count the of graphlets and orbits for a given network

    Args:
        G (nx.MultiDiGraph): Stores the MultiDiGraph representation of our interactome from undirected and directed input files.
        G_prime (nx.Graph): The simplified graph.
        three_node_graphlet_count (dict): Keys are the hash of the graphlet-orbit form: (((0, 1), (0, 1), (1, 1)), 0). Values are the ordered ID values associated with each graphlet.
        three_node_orbit_protein_data (dict): Keys are the hash of the graphlet-orbit form: (((0, 1), (0, 1), (1, 1)), 0). Values are a list of proteins that appear in this orbit.
    Returns:
        three_node_graphlet_count (dict)
        three_node_orbit_protein_data (dict)
        algorithm_run_time (float): Float that stores how long the grphin_algorithm took to execute
    """

    graphlet_config = load_orbit_config("graphlet_config.csv")
    start_time = time.time()

    # create all the binary edge vectors
    adj_list_vector = {node: {} for node in G.nodes()}
    for i, j, data in G.edges(data=True):
        label = data.get("label")
        adj_list_vector[i].setdefault(j, [0, 0, 0])
        adj_list_vector[j].setdefault(i, [0, 0, 0])

        if label == "ppi":
            adj_list_vector[i][j][0] = 1
            adj_list_vector[j][i][0] = 1
        elif label == "reg":
            adj_list_vector[i][j][1] += 1
            adj_list_vector[j][i][2] += 1

    # find all combinations of potential 3 node graphlets
    # pick an edge between A and B
    # for each edge pair, find the union of neighbors between A and B
    three_node_combination = set()  # Use a set for fast triplet lookups

    # Preprocess neighbors for each node once
    neighbors_dict = {i: set(G_prime.neighbors(i)) for i in G_prime.nodes()}
    completed_i = set()
    run_time_data = []

    triple_counter = 0
    count = 0
    for i in G_prime.nodes():
        node_start_time = time.time()
        print(f"Node: {count+1}/{len(G_prime.nodes)}", end="\r")
        for j in neighbors_dict[i]:
            for k in neighbors_dict[j].difference(completed_i):
                if (i != k) and (i != j) and (j != k):

                    triplet = tuple(sorted([i, j, k]))
                    if triplet not in three_node_combination:
                        triple_counter += 1
                        three_node_combination.add(triplet)

                        ab = adj_list_vector[i].get(j, default_edge)
                        ac = adj_list_vector[i].get(k, default_edge)
                        ba = adj_list_vector[j].get(i, default_edge)
                        bc = adj_list_vector[j].get(k, default_edge)
                        ca = adj_list_vector[k].get(i, default_edge)
                        cb = adj_list_vector[k].get(j, default_edge)

                        a_b, a_c, b_a, b_c, c_a, c_b = get_three_node_graphlet_dict(
                            ab, ac, ba, bc, ca, cb
                        )

                        # Order A, B, C edge values internally using min() and max() instead of sorted()
                        a_edges = (min(a_b, a_c), max(a_b, a_c))
                        b_edges = (min(b_a, b_c), max(b_a, b_c))
                        c_edges = (min(c_a, c_b), max(c_a, c_b))

                        # Sort the tuples efficiently
                        sorted_tuples = tuple(sorted([a_edges, b_edges, c_edges]))
                        # catch missing graphlets in config
                        if hash(sorted_tuples) not in three_node_graphlet_count:
                            print("Warning! Missing graphlet in config.")
                            print(sorted_tuples)

                        three_node_graphlet_count[hash(sorted_tuples)] += 1

                        # update orbit_dict
                        orbit_change = get_orbit_position_change(
                            a_edges,
                            b_edges,
                            c_edges,
                            graphlet_config[sorted_tuples]["a_expected"],
                            graphlet_config[sorted_tuples]["b_expected"],
                            graphlet_config[sorted_tuples]["c_expected"],
                            i,
                            j,
                            k,
                        )
                        for idx, orbit in enumerate(
                            graphlet_config[sorted_tuples]["orbits"]
                        ):
                            if orbit == -1:  # Skip missing orbits
                                continue
                            graphlet_key = (sorted_tuples, orbit)

                            # catch missing orbits in config
                            if hash(graphlet_key) not in three_node_orbit_protein_data:
                                print("Warning! Missing orbit in config.")

                            three_node_orbit_protein_data[hash(graphlet_key)] += [
                                orbit_change[idx]
                            ]
        run_time_data.append(time.time() - node_start_time)
        # Once we're done processing i, mark it as completed
        completed_i.add(i)
        count += 1

    # print("triple counter", triple_counter)

    algorithm_run_time = time.time() - start_time
    print("\nRun time: %.3f seconds" % algorithm_run_time)

    return three_node_graphlet_count, three_node_orbit_protein_data, algorithm_run_time


def get_orbit_per_graphlet(
    orbit_dict,
    sorted_tuples,
    a_edges,
    b_edges,
    c_edges,
    i,
    j,
    k,
    graphlet_config,
):
    """Assign orbit information based on the graphlet configuration."""
    for graphlet in graphlet_config:
        if sorted_tuples == graphlet["key"]:
            orbit_change = get_orbit_position_change(
                a_edges,
                b_edges,
                c_edges,
                graphlet["a_expected"],
                graphlet["b_expected"],
                graphlet["c_expected"],
                i,
                j,
                k,
            )
            for idx, orbit in enumerate(graphlet["orbits"]):
                if orbit == -1:  # Skip missing orbits
                    continue
                graphlet_key = (sorted_tuples, orbit)

                # catch missing orbits in config
                if hash(graphlet_key) not in orbit_dict:
                    print("Warning! Missing orbit in config.")

                orbit_dict[hash(graphlet_key)] += [orbit_change[idx]]

    return orbit_dict


def get_orbit_position_change(
    a_edges, b_edges, c_edges, a_expected, b_expected, c_expected, i, j, k
):
    """
    Given the which i, j, k proteins were chosen, orient them compared to the expected a, b, c
    """
    # a, b, c
    if a_edges == a_expected and b_edges == b_expected and c_edges == c_expected:
        return [i, j, k]
    # a, c, b
    if a_edges == a_expected and b_edges == c_expected and c_edges == b_expected:
        return [i, k, j]
    # b, a, c
    if a_edges == b_expected and b_edges == a_expected and c_edges == c_expected:
        return [j, i, k]
    # b, c, a
    if a_edges == b_expected and b_edges == c_expected and c_edges == a_expected:
        return [j, k, i]
    # c, a, b
    if a_edges == c_expected and b_edges == a_expected and c_edges == b_expected:
        return [k, i, j]
    # c, b, a
    if a_edges == c_expected and b_edges == b_expected and c_edges == a_expected:
        return [k, j, i]
    raise ValueError(
        f"Unexpected edges configuration: a_edges={a_edges}, b_edges={b_edges}, c_edges={c_edges}"
    )


def convert_orbit_protein_dict_to_np_matrix(
    G,
    three_node_orbit_id,
    three_node_orbit_protein_data,
    three_node_orbit_namespace,
    output_dir,
):
    """
    Convert three_node_orbit_protein_data into a numpy matrix where the columns are the orbit IDs, the rows are the protein IDs, and the entries are the counts of the proteins in those orbits.

    Args:
        G (nx.MultiDiGraph): Stores the MultiDiGraph representation of our interactome.
        three_node_orbit_id (dict): Keys are the hash of the graphlet-orbit form :(((0, 1), (0, 1), (1, 1)), 0). Values are the ordered ID values associated with each graphlet.
        three_node_orbit_protein_data (dict): Keys are the hash of the graphlet-orbit form: (((0, 1), (0, 1), (1, 1)), 0). Values are a list of proteins that appear in this orbit.
        **this variable may be redundant** three_node_orbit_namespace (dict): Keys are the hash of the graphlet-orbit form: (((0, 1), (0, 1), (1, 1)), 0). Values are the string name for the graphlet-orbit form: (((0, 5), (0, 7), (4, 7)), 2), Orbit 2.
        output_dir (str): Path to desired output folder
    Returns:
        node_orbit_count_matrix (np.arr)
    """
    node_orbit_out = f'{output_dir}/node_orbit.csv'

    node_orbit_count_matrix = np.zeros(
        (len(G.nodes), len(three_node_orbit_id)), dtype=int
    )

    # reassign orbit_dict where each orbit contains a list of (protein, protein_counts)
    for orbit in three_node_orbit_protein_data:
        protein_list = three_node_orbit_protein_data[orbit]
        protein_count_dict = defaultdict(int)
        for protein in protein_list:
            protein_count_dict[protein] += 1

        protein_count_list = []
        for protein in set(protein_list):
            protein_count_list += [(protein, protein_count_dict[protein])]

        three_node_orbit_protein_data[orbit] = protein_count_list

    for orbit in three_node_orbit_protein_data:
        match = re.match(
            r"(\(\(\(.*?\)\), \d+\)), Orbit \d+", three_node_orbit_namespace[orbit]
        )
        if match:
            result = match.group(1)
            try:
                result_tuple = ast.literal_eval(result)
            except (SyntaxError, ValueError) as e:
                print(f"Error converting result to tuple: {e}")
                result_tuple = None
            if result_tuple is not None:
                for node_count_set in three_node_orbit_protein_data[orbit]:
                    if hash(result_tuple) in three_node_orbit_id:
                        orbit_index = three_node_orbit_id[hash(result_tuple)]
                        node_orbit_count_matrix[node_count_set[0]][int(orbit_index)] = (
                            node_count_set[1]
                        )
    print(f"\nNode orbits file saved to: {node_orbit_out}")

    np.savetxt(
        node_orbit_out,
        node_orbit_count_matrix,
        delimiter=",",
        fmt="%d",
    )

    return node_orbit_count_matrix


def initialize_three_node_graphlet_data(
    graphlet_config,
    three_node_graphlet_id,
    three_node_orbit_id,
    three_node_graphlet_count,
    three_node_orbit_protein_data,
    three_node_graphlet_namespace,
    three_node_orbit_namespace,
):
    """
    Using the graphlet_config information, initialize dictionary variables with the correct keys and default values.

    Parameters:
        graphlet_config (dict): Keys are graphlet in the form of ((0, 1), (0, 1), (1, 1)) for example. Values are "key": ast.literal_eval(row["key"]),"a_expected": (tuple),"b_expected": (tuple),"c_expected": (tuple),"orbits": (tuple).
        three_node_graphlet_id (dict): Keys are hash of the graphlet form: ((0, 1), (0, 1), (1, 1)). Values are the ordered ID values associated with each node orbit in the graphlet.
        three_node_orbit_id (dict): Keys are hash of the graphlet form: ((0, 1), (0, 1), (1, 1)). Values are the ordered ID values associated with each graphlet.
        three_node_graphlet_count (dict): Keys are the hash of the graphlet-orbit form: (((0, 1), (0, 1), (1, 1)), 0). Values the ordered ID values associated with each graphlet.
        three_node_orbit_protein_data (dict): Keys are the hash of the graphlet-orbit form: (((0, 1), (0, 1), (1, 1)), 0). Values are a list of proteins that appear in this orbit.
        **this variable may be redundant** three_node_graphlet_namespace (dict): Keys are the hash of the graphlet form: ((0, 1), (0, 1), (1, 1)). Values are the string name for the graphlet form: ((0, 1), (0, 1), (1, 1)).
        **this variable may be redundant** three_node_orbit_namespace (dict): Keys are the hash of the graphlet-orbit form: (((0, 1), (0, 1), (1, 1)), 0). Values are the string name for the graphlet-orbit form: (((0, 5), (0, 7), (4, 7)), 2), Orbit 2.
    Returns:
        three_node_graphlet_id (dict)
        three_node_orbit_id (dict)
        three_node_graphlet_count (dict)
        three_node_orbit_protein_data (dict)
        three_node_graphlet_namespace (dict)
        three_node_orbit_namespace (dict)
    """

    orbit_count = 0
    graphlet_count = 0
    for graphlet in graphlet_config:
        three_node_graphlet_id[hash(graphlet["key"])] = graphlet_count
        three_node_graphlet_count[hash(graphlet["key"])] = 0
        three_node_graphlet_namespace[hash(graphlet["key"])] = graphlet["key"]
        for i in set(graphlet["orbits"]):
            three_node_orbit_id[hash((graphlet["key"], i))] = orbit_count
            three_node_orbit_namespace[hash((graphlet["key"], i))] = str(
                f"{(graphlet['key'], i)}, Orbit {i}"
            )
            three_node_orbit_protein_data[hash((graphlet["key"], i))] = []
            orbit_count += 1
        graphlet_count += 1

    return (
        three_node_graphlet_id,
        three_node_orbit_id,
        three_node_graphlet_count,
        three_node_orbit_protein_data,
        three_node_graphlet_namespace,
        three_node_orbit_namespace,
    )


def write_files(
    G,
    G_prime,
    run_time,
    three_node_graphlet_count,
    three_node_graphlet_namespace,
    three_node_orbit_protein_data,
    three_node_orbit_namespace,
    three_node_orbit_id,
    three_node_graphlet_id,
    protein_id,
    output_dir,
):
    """
    A function to write GRPhIN output files.

    Parameters:
        graphlet_config (dict): keys are graphlet in the form of ((0, 1), (0, 1), (1, 1)) for example. values are "key": ast.literal_eval(row["key"]),"a_expected": (tuple),"b_expected": (tuple),"c_expected": (tuple),"orbits": (tuple)
        G (nx.MultiDiGraph): Stores the graph representation of our mixed interaction network.
        G_prime (nx.Graph): The simplified undirected graph.
        run_time (float): Float that stores how long the grphin_algorithm took to execute.
        three_node_graphlet_count (dict): Keys are the hash of the graphlet-orbit form (((0, 1), (0, 1), (1, 1)), 0). Values are the ordered ID values associated with each graphlet.
        **this variable may be redundant** three_node_graphlet_namespace (dict): Keys are the hash of the graphlet form: ((0, 1), (0, 1), (1, 1)). Values are the string name for the graphlet form: ((0, 1), (0, 1), (1, 1)).
        three_node_orbit_protein_data (dict): Keys are the hash of the graphlet-orbit form: (((0, 1), (0, 1), (1, 1)), 0). Values are a list of proteins that appear in this orbit.
        **this variable may be redundant** three_node_orbit_namespace (dict): Keys are the hash of the graphlet-orbit form: (((0, 1), (0, 1), (1, 1)), 0). Values are the string name for the graphlet-orbit form: (((0, 5), (0, 7), (4, 7)), 2), Orbit 2.
        protein_id (dict): Keys are protein names and values are ids (int).
        output_dir (str): Path to output folder.
    Returns:
        - stats.csv
        - protein_id_mapper.csv
        - orbit_id_mapper.csv
        - orbit_hash_mapper.csv
        - graphlet_id_mapper.csv
        - graphlet_hash_mapper.csv
        - graphlet_counts.csv
    """

    print("Getting statitics...")
    with open(f"{output_dir}/stats.csv", "w+") as f:
        f.write(f"Number of nodes: {len(G.nodes())}\n")
        f.write(f"Number of edges: {len(G.edges())}\n")
        f.write(f"Simplified Graph:\n")
        f.write(f"Number of nodes: {len(G_prime.nodes())}\n")
        f.write(f"Number of edges: {len(G_prime.edges())}\n")
        f.write(f"Run time : %.3f seconds\n" % run_time)
        f.write("Three node graphlet counts\n")
        count = 0
        for key in three_node_graphlet_count:
            f.write(
                f"{three_node_graphlet_namespace[key]} = {three_node_graphlet_count[key]}\n"
            )
            count += three_node_graphlet_count[key]
        f.write(f"Total graphlets found: {count}\n")
        f.write(f"Unique graphlet counts : {len(three_node_graphlet_count)}\n")
        f.write(f"Three node orbit counts\n")
        total_orbit_count = 0

        for orbit in three_node_orbit_protein_data:
            orbit_count = 0
            for protein_set in three_node_orbit_protein_data[orbit]:
                orbit_count += protein_set[1]
            total_orbit_count += orbit_count
            f.write(f"{three_node_orbit_namespace[orbit]} : {orbit_count}\n")
        f.write(f"Total orbits found: {total_orbit_count}\n")
        f.write(f"Unique orbits counts : {len(three_node_orbit_protein_data)}\n")
    f.close()

    with open(f"{output_dir}/protein_id_mapper.csv", "w+") as f:
        for protein in protein_id:
            f.write(f"{protein}, {protein_id[protein]}\n")
    f.close()

    with open(f"{output_dir}/orbit_id_mapper.csv", "w+") as f:
        for orbit in three_node_orbit_id:
            if orbit in three_node_orbit_namespace:
                f.write(
                    f"{three_node_orbit_namespace[orbit]} : {three_node_orbit_id[orbit]}\n"
                )

    with open(f"{output_dir}/orbit_hash_mapper.csv", "w+") as f:
        for orbit in three_node_orbit_id:
            if orbit in three_node_orbit_namespace:
                f.write(f"{orbit} : {three_node_orbit_id[orbit]}\n")

    with open(f"{output_dir}/graphlet_id_mapper.csv", "w+") as f:
        for graphlet in three_node_graphlet_id:
            if graphlet in three_node_graphlet_namespace:
                f.write(
                    f"{three_node_graphlet_namespace[graphlet]} : {three_node_graphlet_id[graphlet]}\n"
                )

    with open(f"{output_dir}/graphlet_hash_mapper.csv", "w+") as f:
        for graphlet in three_node_graphlet_count:
            if graphlet in three_node_graphlet_namespace:
                f.write(f"{graphlet} : {three_node_graphlet_namespace[graphlet]}\n")

    with open(f"{output_dir}/graphlet_counts.csv", "w+") as f:
        for graphlet in three_node_graphlet_id:
            if graphlet in three_node_graphlet_namespace:
                f.write(
                    f"{three_node_graphlet_namespace[graphlet]} : {three_node_graphlet_count[graphlet]}\n"
                )

def write_graphlet_counts(
        three_node_graphlet_id,
        three_node_graphlet_namespace,
        three_node_graphlet_count,
        output_dir,
        num_random_nets
        ):
    """
    For use in "graphlets only" mode. Writes graphlet counts to a CSV file. For use in generating graphlet counts from randomized networks.
    
    Parameters:
        graphlet_config (dict): keys are graphlet in the form of ((0, 1), (0, 1), (1, 1)) for example. values are "key": ast.literal_eval(row["key"]),"a_expected": (tuple),"b_expected": (tuple),"c_expected": (tuple),"orbits": (tuple)
        three_node_graphlet_count (dict): Keys are the hash of the graphlet-orbit form (((0, 1), (0, 1), (1, 1)), 0). Values are the ordered ID values associated with each graphlet.
        three_node_graphlet_namespace (dict): Keys are the hash of the graphlet form: ((0, 1), (0, 1), (1, 1)). Values are the string name for the graphlet form: ((0, 1), (0, 1), (1, 1)).
        output_dir (str): Path to output folder.
        num_random_nets (int): The number of random randomized networks to run GRPhIN on.
    Returns:
        - graphlet_counts{i}.csv
    """

    for i in range(num_random_nets):
        with open(f"{output_dir}/graphlet_counts{i}.csv", "w+") as f:
            for graphlet in three_node_graphlet_id:
                if graphlet in three_node_graphlet_namespace:
                    f.write(
                        f"{three_node_graphlet_namespace[graphlet]} : {three_node_graphlet_count[graphlet]}\n"
                    )


def count_three_node_graphlets(graphlet_config, protein_id, G, G_prime, output_dir, graphlets_only):
    """
    Wrapper function for variable initialization, running GRPhIN algorithm, and outputing result files for 3-node graphlets.

    Args:
        graphlet_config (dict): Keys are graphlet in the form of: ((0, 1), (0, 1), (1, 1)). Values are "key": ast.literal_eval(row["key"]),"a_expected": (tuple),"b_expected": (tuple),"c_expected": (tuple),"orbits": (tuple).
        protein_id (dict): Keys are protein names and values are IDs (int).
        G (nx.MultiDiGraph): Stores the graph representation of our mixed interaction network.
        G_prime (nx.Graph): The simplified undirected graph.
        output_dir (str): String to output folder.
        graphlets_only (bool): Boolean indicating whether to run GRPhIN in graphlets only mode.

    Returns:
        None
    """

    # initialize variables required for counting and storing graphlet and orbit information
    three_node_graphlet_count = {}
    three_node_graphlet_namespace = {}
    three_node_orbit_protein_data = {}
    three_node_orbit_namespace = {}
    three_node_graphlet_id = {}
    three_node_orbit_id = {}
    run_time = 0

    # using the graphlet_config to initialize dictionaries with information
    (
        three_node_graphlet_id,
        three_node_orbit_id,
        three_node_graphlet_count,
        three_node_orbit_protein_data,
        three_node_graphlet_namespace,
        three_node_orbit_namespace,
    ) = initialize_three_node_graphlet_data(
        graphlet_config,
        three_node_graphlet_id,
        three_node_orbit_id,
        three_node_graphlet_count,
        three_node_orbit_protein_data,
        three_node_graphlet_namespace,
        three_node_orbit_namespace,
    )

    # use the grphin algorithm to count graphlets and orbits
    (three_node_graphlet_count, three_node_orbit_protein_data, run_time) = (
        grphin_algorithm(
            G, G_prime, three_node_graphlet_count, three_node_orbit_protein_data
        )
    )

    node_orbit_count_matrix = convert_orbit_protein_dict_to_np_matrix(
        G,
        three_node_orbit_id,
        three_node_orbit_protein_data,
        three_node_orbit_namespace,
        output_dir,
    )

    # Output necessary files
    if graphlets_only == True:
        write_graphlet_counts(
            three_node_graphlet_id,
            three_node_graphlet_namespace,
            three_node_graphlet_count,
            output_dir,
            1,
        )
    else: 
        write_files(
            G,
            G_prime,
            run_time,
            three_node_graphlet_count,
            three_node_graphlet_namespace,
            three_node_orbit_protein_data,
            three_node_orbit_namespace,
            three_node_orbit_id,
            three_node_graphlet_id,
            protein_id,
            output_dir,
        )

def plot_runtime_stats():
    """
    Plots runtime statistics for internal benchmarking and optimization.
    """

    species = ["bsub", "cerevisiae", "drerio", "elegans", "fly"]
    runtime_full_algorithm_data = [6.132, 4507.445, 407.541, 303.814, 340.379]
    runtime_triplet_iteration_data = [0.461, 194.703, 51.803, 24.740, 29.512]

    stacked_bar_data = {"below" : runtime_full_algorithm_data, "above" : runtime_triplet_iteration_data}

    fig, ax = plt.subplots()
    width = 0.5
    bottom = np.zeros(len(species))

    for stack, data in stacked_bar_data.items():
        p = ax.bar(species, data, width=width, label=stack, bottom=bottom)
        bottom += data

    plt.show()


def main(input_ppi, input_reg, output_dir, graphlets_only=False):
    """
    A function that generates node orbit counts, graphlet counts and summary statistics for a mixed network from a directed and an undirected edge file.
    
    Parameters:
        -u / --undirected: Command-line argument for undirected edges input file (Example: PPI edges).
        -d / --directed: Command-line argument for directed edges input file (Example: Regulatory edges).
        -o / --output_dir: Command-line argument for output directory.
        -g / --graphlets_only: Command-line argument to run GRPhIN in graphlets only mode.

    Returns:
        Counts of graphlets and node orbit positions within a mixed graph.

    Example:
        python3 grphin.py -u data/bsub_ppi.csv -d data/bsub_reg.csv -o output/bsub

        This will count the occurences of all 2 and 3-node mixed graphlets in the mixed network generated from the two input files. It will also output node orbit positions and various statistics to the specified output directory.
    """

    # Initialize graphlet data with config file
    protein_id, G, G_prime, graphlet_config = initialize_graphlet_data(
        input_ppi, input_reg
    )

    # Print graph information
    print(f"Complete Graph:")
    print(f"Number of nodes: {len(G.nodes())}")
    print(f"Number of edges: {len(G.edges())}")
    print(f"\nSimplified Graph:")
    print(f"Number of nodes: {len(G_prime.nodes())}")
    print(f"Number of edges: {len(G_prime.edges())}")

    if graphlets_only:
        print("Graphlets only mode enabled.")
        # Count three-node graphlets
        count_three_node_graphlets(graphlet_config, protein_id, G, G_prime, output_dir, graphlets_only)
    else:
        # Count two-node graphlets
        count_two_node_graphlet(G, output_dir)

        # Count three-node graphlets
        count_three_node_graphlets(graphlet_config, protein_id, G, G_prime, output_dir, graphlets_only)

    print("GRPhIN Algorithm finished successfully.")

    ## LEAVE THIS LINE UNCOMMENTED TO PLOT RUNTIME STATISTICS
    # plot_runtime_stats()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="GRPhIN Algorithm"
    )
    parser.add_argument(
        "-u",
        "--undirected", 
        type=str, 
        help="Path to the undirected/PPI edges input file.",
        required=True
    )
    parser.add_argument(
        "-d",
        "--directed",
        type=str, 
        help="Path to the directed/Reg edges input file.",
        required=True
    )
    parser.add_argument(
        "-o",
        "--output_dir", 
        type=str, 
        help="Path to the output directory.",
        required=True
    )
    parser.add_argument(
        "-g",
        "--graphlets_only", 
        type=bool,
        help="Run GRPhIN in graphlets only mode."
    )

    args = parser.parse_args()

    main(args.undirected, args.directed, args.output_dir, args.graphlets_only)