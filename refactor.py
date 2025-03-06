import ast
from collections import defaultdict
import re
import sys
import time
from matplotlib import pyplot as plt
import networkx as nx
import csv
from pathlib import Path
import numpy as np
from collections import defaultdict


def count_two_node_graphlet(G, output_dir, species):

    two_node_graphlet_count, two_node_graphlet_id = initialize_two_node_graphlet_data()

    two_node_graphlet_count, two_node_orbit_dict = get_two_node_graphlet_stats(
        G, two_node_graphlet_count
    )

    generate_two_node_output_files(
        two_node_graphlet_id,
        two_node_graphlet_count,
        two_node_orbit_dict,
        output_dir,
        species,
    )

    return two_node_graphlet_count, two_node_graphlet_id, two_node_orbit_dict


def generate_two_node_output_files(
    two_node_graphlet_id,
    two_node_graphlet_count,
    two_node_orbit_dict,
    output_dir,
    species,
):
    print("\ntwo node graphlet counts")
    for key in two_node_graphlet_id:
        print(f"G_{two_node_graphlet_id[key]} = {two_node_graphlet_count[key]}")

    print("\ntwo node graphlet orbit counts")
    for key in two_node_orbit_dict:
        print(f"{key} = {len(two_node_orbit_dict[key])}")

    with open(f"{output_dir}/{species}/two_node_graphlet_counts.csv", "w") as f:
        for key in two_node_graphlet_id:
            f.write(f"G_{two_node_graphlet_id[key]}, {two_node_graphlet_count[key]}\n")
    f.close()

    with open(f"{output_dir}/{species}/two_node_orbit_counts.csv", "w") as f:
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
    print(f"processing {label} edges")
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
    # see diagram for reference
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
    """Load graphlet lookup table from a CSV file."""
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


def get_two_node_graphlet_stats(G, two_node_graphlet_dict):
    G_adj_list = get_two_node_adjacency_list(G)
    orbits_dict = {1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: []}
    for neighbors in G_adj_list:
        i = neighbors[0]
        for j in neighbors[1]:
            vectors = G_adj_list[i][1][j] + G_adj_list[j][1][i]
            orbits_dict = get_two_node_orbit_position(
                i, j, G_adj_list[i][1][j], G_adj_list[j][1][i], orbits_dict
            )
            if hash(tuple(vectors)) in two_node_graphlet_dict:
                two_node_graphlet_dict[hash(tuple(vectors))] += 1
    return two_node_graphlet_dict, orbits_dict


def get_two_node_adjacency_list(G):
    """Get the adjacency list for a MultiDiGraph"""
    print("\ngetting adjacency list")

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

    # check orbit positioning via vectors
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
    return hash(tuple(xy))


def get_three_node_graphlet_dict(ab, ac, ba, bc, ca, cb):
    # Process into a hashed tuple
    ab = hash_tuple(ab)
    ac = hash_tuple(ac)
    ba = hash_tuple(ba)
    bc = hash_tuple(bc)
    ca = hash_tuple(ca)
    cb = hash_tuple(cb)

    # Define the dictionary
    my_dict = {
        hash((0, 0, 0)): 0,
        hash((1, 0, 0)): 1,
        hash((0, 1, 0)): 2,
        hash((0, 0, 1)): 3,
        hash((1, 1, 0)): 4,
        hash((1, 0, 1)): 5,
        hash((0, 1, 1)): 6,
        hash((1, 1, 1)): 7,
    }

    return my_dict[ab], my_dict[ac], my_dict[ba], my_dict[bc], my_dict[ca], my_dict[cb]


def plot_run_time_data(run_time_data):
    x_list = [i for i in range(len(run_time_data))]

    norm_run_time = [
        (x - min(run_time_data)) / (max(run_time_data) - min(run_time_data))
        for x in run_time_data
    ]

    plt.figure(figsize=(14, 8))  # Correcting figure size
    plt.plot(x_list, norm_run_time, "o", markersize=2)
    plt.xlabel("Index")  # Adding labels for clarity
    plt.ylabel("Run Time")
    plt.title("Run Time Data Plot")
    plt.show()


def compare_csv_files(file_path1, file_path2):
    """
    Compares two CSV files and returns a list of differences.
    """
    differences = []
    with open(file_path1, "r") as file1, open(file_path2, "r") as file2:
        reader1 = csv.reader(file1)
        reader2 = csv.reader(file2)

        for row_num, (row1, row2) in enumerate(zip(reader1, reader2), start=1):
            if row1 != row2:
                differences.append(f"Row {row_num}: File1 - {row1}, File2 - {row2}")

    if len(differences) == 0:
        print("\033[42;37mPASSED TEST\033[0m")
    else:
        print("\033[41;37mFAILED TEST\033[0m")

    return None


def grphin_algorithm(
    G: nx.MultiDiGraph, G_prime: nx.Graph, three_node_graphlet_dict, orbit_dict
):
    graphlet_config = load_graphlet_config("graphlet_config.csv")
    start_time = time.time()

    # create all the binary edge vectors
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

    # find all combinations of potential 3 node graphlets
    # pick an edge between A and B
    # for each edge pair, find the union of neighbors between A and B
    three_node_combination = set()  # Use a set for fast triplet lookups

    # Preprocess neighbors for each node once
    neighbors_dict = {i: set(G_prime.neighbors(i)) for i in G_prime.nodes()}
    completed_i = set()
    run_time_data = []

    node_list = [node for node, _ in sorted(G_prime.degree(), key=lambda x: x[1], reverse=False)]

    print("len of g prime deg", len(node_list))
    print("len of g prime nodes", len(G_prime.nodes()))

    triple_counter = 0
    count = 0
    # for i in node_list:
    for i in G_prime.nodes():
        node_start_time = time.time()
        print(f"Node: {count}/{len(G_prime.nodes)}", end="\r")
        for j in neighbors_dict[i]:
            for k in neighbors_dict[j]:
                if (
                    (i != k) and (i != j) and (j != k)
                ):
                    
                    triplet = tuple(sorted([i, j, k]))
                    if triplet not in three_node_combination:
                        triple_counter+=1
                        three_node_combination.add(triplet)

                        a = i
                        b = j
                        c = k

                        ab = ac = ba = bc = ca = cb = 0
                        if b in adj_list_vector[a]:
                            ab = adj_list_vector[a][b]
                        else:
                            ab = [0, 0, 0]
                        if c in adj_list_vector[a]:
                            ac = adj_list_vector[a][c]
                        else:
                            ac = [0, 0, 0]
                        if a in adj_list_vector[b]:
                            ba = adj_list_vector[b][a]
                        else:
                            ba = [0, 0, 0]
                        if c in adj_list_vector[b]:
                            bc = adj_list_vector[b][c]
                        else:
                            bc = [0, 0, 0]
                        if a in adj_list_vector[c]:
                            ca = adj_list_vector[c][a]
                        else:
                            ca = [0, 0, 0]
                        if b in adj_list_vector[c]:
                            cb = adj_list_vector[c][b]
                        else:
                            cb = [0, 0, 0]
                        a_b, a_c, b_a, b_c, c_a, c_b = get_three_node_graphlet_dict(
                            ab, ac, ba, bc, ca, cb
                        )

                        # order A, B, C edge values internally
                        a_edges = tuple(sorted([a_b, a_c]))
                        b_edges = tuple(sorted([b_a, b_c]))
                        c_edges = tuple(sorted([c_a, c_b]))

                        # Create a list of tuples in order [A, B, C]
                        tuples_list = [a_edges, b_edges, c_edges]

                        # Sort the tuples first by the first index, then by the second index
                        sorted_tuples = tuple(
                            sorted(tuples_list, key=lambda x: (x[0], x[1]))
                        )

                        # catch missing graphlets in config
                        if hash(sorted_tuples) not in three_node_graphlet_dict:
                            print("MISSING GRAPHLET IN CONFIG")
                            print(sorted_tuples)

                        three_node_graphlet_dict[hash(sorted_tuples)] += 1
                        orbit_dict = get_orbit_per_graphlet(
                            orbit_dict,
                            sorted_tuples,
                            a_edges,
                            b_edges,
                            c_edges,
                            i,
                            j,
                            k,
                            graphlet_config,
                        )
                    var = 0
        run_time_data.append(time.time() - node_start_time)
        # Once we're done processing i, mark it as completed
        completed_i.add(i)
        count+=1

    print("triple counter", triple_counter)

    algorithm_run_time = time.time() - start_time
    print("run time : %.3f seconds" % algorithm_run_time)
    # plot_run_time_data(run_time_data)

    return three_node_graphlet_dict, orbit_dict, algorithm_run_time


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
                    print("MISSING ORBIT IN CONFIG")

                orbit_dict[hash(graphlet_key)] += [orbit_change[idx]]

    return orbit_dict


def get_orbit_position_change(
    a_edges, b_edges, c_edges, a_expected, b_expected, c_expected, i, j, k
):
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


def get_node_orbit_matrix(
    G,
    three_node_orbit_id,
    three_node_orbit_protein_data,
    three_node_orbit_namespace,
    output_dir,
    species,
):

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

    np.savetxt(
        f"{output_dir}/{species}/node_orbit.csv",
        node_orbit_count_matrix,
        delimiter=",",
        fmt="%d",
    )

    compare_csv_files(
        Path("final_output/bsub/node_orbit.csv"),
        Path("output_refactor/bsub/node_orbit.csv"),
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


def write_stats(
    G,
    G_prime,
    run_time,
    three_node_graphlet_count,
    three_node_graphlet_namespace,
    three_node_orbit_protein_data,
    three_node_orbit_namespace,
    output_dir,
):
    print("getting stats")
    with open(f"{output_dir}/bsub/stats.csv", "w") as f:
        f.write(f"Number of nodes: {len(G.nodes())}\n")
        f.write(f"Number of edges: {len(G.edges())}\n")
        f.write(f"Simplified Graph:\n")
        f.write(f"Number of nodes: {len(G_prime.nodes())}\n")
        f.write(f"Number of edges: {len(G_prime.edges())}\n")
        f.write(f"run time : %.3f seconds\n" % run_time)
        f.write("three node graphlet counts\n")
        count = 0
        for key in three_node_graphlet_count:
            f.write(
                f"{three_node_graphlet_namespace[key]} = {three_node_graphlet_count[key]}\n"
            )
            count += three_node_graphlet_count[key]
        f.write(f"Total graphlets found: {count}\n")
        f.write(f"unique graphlet counts : {len(three_node_graphlet_count)}\n")
        f.write(f"three node orbit counts\n")
        total_orbit_count = 0

        for orbit in three_node_orbit_protein_data:
            orbit_count = 0
            for protein_set in three_node_orbit_protein_data[orbit]:
                orbit_count += protein_set[1]
            total_orbit_count += orbit_count
            f.write(f"{three_node_orbit_namespace[orbit]} : {orbit_count}\n")
        f.write(f"Total orbits found: {total_orbit_count}\n")
        f.write(f"unique orbits counts : {len(three_node_orbit_protein_data)}\n")
    f.close()


def count_three_node_graphlets(graphlet_config, G, G_prime, output_dir, species):

    three_node_graphlet_count = {}
    three_node_graphlet_namespace = {}
    three_node_orbit_protein_data = {}
    three_node_orbit_namespace = {}
    three_node_graphlet_id = {}
    three_node_orbit_id = {}
    run_time = 0
    node_orbit_arr = np.zeros((len(G.nodes), len(three_node_orbit_id)), dtype=int)

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

    (
        three_node_graphlet_count,
        three_node_orbit_protein_data,
        run_time,
    ) = grphin_algorithm(
        G, G_prime, three_node_graphlet_count, three_node_orbit_protein_data
    )

    node_orbit_count_matrix = get_node_orbit_matrix(
        G,
        three_node_orbit_id,
        three_node_orbit_protein_data,
        three_node_orbit_namespace,
        output_dir,
        species,
    )

    # print("GRAPHLET COUNTS")
    # for graphlet in three_node_graphlet_count:
    #     print(three_node_graphlet_namespace[graphlet], three_node_graphlet_count[graphlet])

    # print("\nORBIT COUNTS")
    # for orbit in three_node_orbit_protein_data:
    #     print(three_node_orbit_namespace[orbit], len(three_node_orbit_protein_data[orbit]))

    write_stats(
        G,
        G_prime,
        run_time,
        three_node_graphlet_count,
        three_node_graphlet_namespace,
        three_node_orbit_protein_data,
        three_node_orbit_namespace,
        output_dir,
    )


def main():
    network_ppi_path = Path("data/bsub_ppi.csv")
    network_reg_path = Path("data/bsub_reg.csv")
    # network_ppi_path = Path("data/cerevisiae_ppi.csv")
    # network_reg_path = Path("data/cerevisiae_reg.csv")
    # network_ppi_path = Path("data/fly_ppi.csv")
    # network_reg_path = Path("data/fly_reg.csv")
    # network_ppi_path = Path("data/elegans_ppi.csv")
    # network_reg_path = Path("data/elegans_reg.csv")
    # network_ppi_path = Path("data/drerio_ppi.csv")
    # network_reg_path = Path("data/drerio_reg.csv")
    stress_proteins_path = Path(
        "data/oxidative_stress/txid224308/txid224308-stress-proteins.csv"
    )
    counting_algorithm = True
    output_dir = Path("output_refactor")
    species = "bsub"

    # initialize graphlet data
    protein_id, G, G_prime, graphlet_config = initialize_graphlet_data(
        network_ppi_path, network_reg_path
    )

    # print graph information
    print(f"Complete Graph:")
    print(f"Number of nodes: {len(G.nodes())}")
    print(f"Number of edges: {len(G.edges())}")
    print(f"\nSimplified Graph:")
    print(f"Number of nodes: {len(G_prime.nodes())}")
    print(f"Number of edges: {len(G_prime.edges())}")

    # two_node_graphlet_count, two_node_graphlet_id, two_node_orbit_dict = (
    #     count_two_node_graphlet(G, output_dir, species)
    # )

    count_three_node_graphlets(graphlet_config, G, G_prime, output_dir, species)


if __name__ == "__main__":
    main()
