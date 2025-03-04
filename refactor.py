import ast
import networkx as nx
import csv
from pathlib import Path

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

def main():
    network_ppi_path = Path("data/bsub_ppi.csv")
    network_reg_path = Path("data/bsub_reg.csv")
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

    two_node_graphlet_count, two_node_graphlet_id, two_node_orbit_dict = (
        count_two_node_graphlet(G, output_dir, species)
    )

if __name__ == "__main__":
    main()
