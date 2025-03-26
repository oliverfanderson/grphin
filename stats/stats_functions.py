import ast
import csv
import os
import re
import sys
from matplotlib import pyplot as plt
import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from refactor import load_graphlet_config, initialize_three_node_graphlet_data


def plot_runtime_stats():
    species = ["bsub", "cerevisiae", "drerio", "elegans", "fly"]
    runtime_full_algorithm_data = [6.132, 4507.445, 407.541, 303.814, 340.379]
    runtime_triplet_iteration_data = [0.461, 194.703, 51.803, 24.740, 29.512]

    stacked_bar_data = {
        "below": runtime_full_algorithm_data,
        "above": runtime_triplet_iteration_data,
    }

    fig, ax = plt.subplots()
    width = 0.5
    bottom = np.zeros(len(species))

    for stack, data in stacked_bar_data.items():
        p = ax.bar(species, data, width=width, label=stack, bottom=bottom)
        bottom += data

    plt.show()


def plot_three_node_graphlet_distribution(
    graphlet_dict, graphlet_mapper, indexed_graphlet_dict, selected_network, output_dir
):
    print("plotting graphlet distribution")
    hist_data = []
    x_label = [*range(0, len(indexed_graphlet_dict), 1)]
    for graphlet in indexed_graphlet_dict:
        if graphlet in graphlet_dict:
            hist_data.append(graphlet_dict[graphlet])
        else:
            hist_data.append(0)
    hist_data = [value if value > 0 else 0.1 for value in hist_data]

    fig = plt.figure(figsize=(14, 6))
    plt.bar(x_label, hist_data, color="skyblue", edgecolor="black")
    plt.yscale("log")
    plt.title(f"{selected_network} Graphlet Count Distribution", fontsize=16)
    plt.xlabel("Graphlet Index", fontsize=14)
    plt.ylabel("Count (log scale)", fontsize=14)
    plt.xticks(x_label[::2], fontsize=8)
    plt.grid(axis="y", linestyle="--", alpha=0.7)
    plt.tight_layout()
    plt.savefig(f"{output_dir}/three_node_graphlet_dist.pdf")
    plt.show()

    sorted_graphlet_dict = {
        key: value
        for key, value in sorted(
            graphlet_dict.items(), key=lambda item: item[1], reverse=True
        )
    }

    with open(f"{output_dir}/top_graphlet_counts.csv", "w") as f:
        f.write(f"graphlet\tid\tcount\n")
        for graphlet in sorted_graphlet_dict:
            f.write(
                f"{graphlet_mapper[graphlet]}\t{indexed_graphlet_dict[graphlet]}\t{sorted_graphlet_dict[graphlet]}\n"
            )

    return None


def read_output_files(output_dir):
    three_node_graphlet_count = {}
    three_node_graphlet_namespace = {}
    three_node_orbit_protein_data = {}
    three_node_orbit_namespace = {}
    three_node_graphlet_id = {}
    three_node_orbit_id = {}
    run_time = 0

    graphlet_config = load_graphlet_config("graphlet_config.csv")

    # for graphlet in graphlet_config:
    #     print(graphlet)

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

    print("post-processing only")
    with open(f"{output_dir}/graphlet_hash_mapper.csv", "r") as file:
        csv_reader = csv.reader(file, delimiter=":")
        for row in csv_reader:
            graphlet_hash = int(row[0])
            graphlet = ast.literal_eval(row[1])
            three_node_graphlet_namespace[graphlet_hash] = graphlet
    print("post processing step 1/6")

    # for key in three_node_graphlet_namespace:
    #     print(key, three_node_graphlet_namespace[key])

    with open(f"{output_dir}/graphlet_counts.csv", "r") as file:
        csv_reader = csv.reader(file, delimiter=":")
        for row in csv_reader:
            graphlet = ast.literal_eval(row[0])
            count = int(row[1].strip())
            three_node_graphlet_count[hash(graphlet)] = count
    print("post processing step 2/6")

    # for key in three_node_graphlet_count:
    #     print(key, three_node_graphlet_count[key])

    with open(f"{output_dir}/orbit_hash_mapper.csv", "r") as file:
        csv_reader = csv.reader(file, delimiter=":")
        for row in csv_reader:
            orbit_hash = int(row[0])
            three_node_orbit_namespace[orbit_hash] = row[1]
    print("post processing step 3/6")

    # for key in three_node_orbit_namespace:
    #     print(key, three_node_orbit_namespace[key])

    with open(f"{output_dir}/graphlet_id_mapper.csv", "r") as file:
        csv_reader = csv.reader(file, delimiter=":")
        for row in csv_reader:
            graphlet = ast.literal_eval(row[0])
            id = row[1].strip()
            three_node_graphlet_id[hash(graphlet)] = id
    print("post processing step 4/6")

    # for key in three_node_graphlet_id:
    #     print(key, three_node_graphlet_id[key])

    orbit_id_dict = {}
    with open(f"{output_dir}/orbit_id_mapper.csv", "r") as file:
        csv_reader = csv.reader(file, delimiter=":")
        for row in csv_reader:
            match = re.search(r"\(\(\(.*?\)\), \d+\)", row[0])
            if match:
                result = ast.literal_eval(match.group())
                id = row[1].strip()
                three_node_orbit_id[hash(result)] = id
                orbit_id_dict[id] = result
    print("post processing step 5/6")

    # for key in three_node_orbit_id:
    #     print(key, three_node_orbit_id[key])

    # for key in orbit_id_dict:
    #     print(key, orbit_id_dict[key])

    node_orbit_arr = np.loadtxt(
        f"{output_dir}/node_orbit.csv", delimiter=",", dtype=int
    )
    rows, cols = node_orbit_arr.shape

    count = 0
    for orbit in orbit_id_dict:
        for protein in range(0, rows, 1):
            orbit_hash = hash(orbit_id_dict[orbit])
            if orbit_hash not in three_node_orbit_protein_data:
                three_node_orbit_protein_data[orbit_hash] = []
            protein_count = node_orbit_arr[protein][int(orbit)]
            three_node_orbit_protein_data[orbit_hash] += [(protein, protein_count)]
            print(count, "/", rows * cols, end="\r")
            count += 1

    # for key in three_node_orbit_protein_data:
    #     print(key, three_node_orbit_protein_data[key][1])

    return (
        three_node_graphlet_count,
        three_node_graphlet_namespace,
        three_node_orbit_protein_data,
        three_node_orbit_namespace,
        three_node_graphlet_id,
        three_node_orbit_id,
        graphlet_config,
    )


def main():
    print("running stats")
    read_output_files("output_refactor/bsub")


if __name__ == "__main__":
    main()
