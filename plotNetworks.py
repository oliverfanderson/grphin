import argparse
import networkx as nx
import csv
from pyvis.network import Network

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
            geneName1 = row[2]
            geneName2 = row[3]

            if geneName1 not in visited_nodes:
                visited_nodes.add(geneName1)
                node_count += 1

            if geneName2 not in visited_nodes:
                visited_nodes.add(geneName2)
                node_count += 1

            G.add_edge(geneName1, geneName2, label=label)
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

    return G

def main(ppi_path, reg_path, output_dir):
    G = read_csv(
            ppi_path,
            reg_path
        )
    
    nt = Network('700px', '840px', directed=True)
    # nt = Network(directed=True)
    G.remove_edges_from(nx.selfloop_edges(G))

    # Add nodes to the network graph
    for node in G.nodes():
        nt.add_node(node, label=node)  # Add node label
        nt.get_node(node)['color'] = '#cfe2f3'
        nt.get_node(node)['shape'] = 'ellipse'

    if 'txid224308' in args.output_dir:
        nt.get_node('sigA')['color'] = '#fad6a5'
        nt.get_node('sigB')['color'] = '#fad6a5'
        nt.get_node('katX')['color'] = '#fad6a5'
        nt.get_node('katE')['color'] = '#fad6a5'
        nt.get_node('ahpC')['color'] = '#fad6a5'
        nt.get_node('ahpF')['color'] = '#fad6a5'
        nt.get_node('ccpA')['color'] = '#fad6a5'
        nt.get_node('perR')['color'] = '#fad6a5'
    if 'txid7955' in args.output_dir:
        nt.get_node('nanog')['color'] = '#fad6a5'
        nt.get_node('sod1')['color'] = '#fad6a5'
        nt.get_node('sod2')['color'] = '#fad6a5'
        nt.get_node('prdx1')['color'] = '#fad6a5'
        nt.get_node('erp44')['color'] = '#fad6a5'
        nt.get_node('keap1a')['color'] = '#fad6a5'
        nt.get_node('keap1b')['color'] = '#fad6a5'
        nt.get_node('smad2')['color'] = '#fad6a5'
        nt.get_node('gata1a')['color'] = '#fad6a5'
        nt.get_node('sall4')['color'] = '#fad6a5'
        nt.get_node('foxh1')['color'] = '#fad6a5'  # Skip this subnetwork if not in the specified list

    # Track added edges to prevent duplication
    seen_ppi_edges = set()

    # Iterate over edges
    for u, v, data in G.edges(data=True):
        label = data.get("label", None)  # Get label from edge data
        if label == "reg":
            nt.add_edge(u, v, color="red", width = 3, arrows="to", arrowStrikethrough=True)  # Directed red edge
        elif label == "ppi":
            if (u, v) not in seen_ppi_edges and (v, u) not in seen_ppi_edges:
                nt.add_edge(u, v, color="black", width = 4, arrows="")  # Undirected black edge
                seen_ppi_edges.add((u, v))  # Mark the pair as seen
                seen_ppi_edges.add((v, u))  # Ensure (v, u) isn't added separately

    nt.write_html(f'{output_dir}nx.html')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot Stress Subnetworks")
    parser.add_argument(
        "-u",
        "--undirected",
        type=str,
        help="Path to the undirected/PPI edges input file.",
        required=True,
    )
    parser.add_argument(
        "-d",
        "--directed",
        type=str,
        help="Path to the directed/Reg edges input file.",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output_dir",
        type=str,
        help="Path to the output directory.",
        required=True,
    )

    args = parser.parse_args()

    main(args.undirected, args.directed, args.output_dir)