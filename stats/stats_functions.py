import ast
from collections import defaultdict
import csv
import os
from pathlib import Path
import random
import re
import sys
from matplotlib import pyplot as plt
import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from refactor import (
    initialize_graphlet_data,
    load_graphlet_config,
    initialize_three_node_graphlet_data,
)


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
    three_node_graphlet_count,
    three_node_graphlet_namespace,
    three_node_graphlet_id,
    species,
    output_dir,
):
    print("plotting graphlet distribution")
    hist_data = []
    x_label = [*range(0, len(three_node_graphlet_id), 1)]
    for graphlet in three_node_graphlet_id:
        if graphlet in three_node_graphlet_count:
            hist_data.append(three_node_graphlet_count[graphlet])
        else:
            hist_data.append(0)
    hist_data = [value if value > 0 else 0.1 for value in hist_data]

    fig = plt.figure(figsize=(14, 6))
    plt.bar(x_label, hist_data, color="skyblue", edgecolor="black")
    plt.yscale("log")
    plt.title(f"{species} Graphlet Count Distribution", fontsize=16)
    plt.xlabel("Graphlet Index", fontsize=14)
    plt.ylabel("Count (log scale)", fontsize=14)
    plt.xticks(x_label[::2], fontsize=10)
    plt.grid(axis="y", linestyle="--", alpha=0.7)
    plt.tight_layout()
    plt.savefig(f"{output_dir}/three_node_graphlet_dist.pdf")
    # plt.show()

    sorted_graphlet_dict = {
        key: value
        for key, value in sorted(
            three_node_graphlet_count.items(), key=lambda item: item[1], reverse=True
        )
    }

    with open(f"{output_dir}/top_graphlet_counts.csv", "w") as f:
        f.write(f"graphlet\tid\tcount\n")
        for graphlet in sorted_graphlet_dict:
            f.write(
                f"{three_node_graphlet_namespace[graphlet]}\t{three_node_graphlet_id[graphlet]}\t{sorted_graphlet_dict[graphlet]}\n"
            )

    return None


def read_output_files(output_dir, species):
    three_node_graphlet_count = {}
    three_node_graphlet_namespace = {}
    three_node_orbit_protein_data = {}
    three_node_orbit_namespace = {}
    three_node_graphlet_id = {}
    three_node_orbit_id = {}
    run_time = 0

    graphlet_config = load_graphlet_config("graphlet_config.csv")

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

    with open(f"{output_dir}/graphlet_counts.csv", "r") as file:
        csv_reader = csv.reader(file, delimiter=":")
        for row in csv_reader:
            graphlet = ast.literal_eval(row[0])
            count = int(row[1].strip())
            three_node_graphlet_count[hash(graphlet)] = count
    print("post processing step 2/6")

    with open(f"{output_dir}/orbit_hash_mapper.csv", "r") as file:
        csv_reader = csv.reader(file, delimiter=":")
        for row in csv_reader:
            orbit_hash = int(row[0])
            three_node_orbit_namespace[orbit_hash] = row[1]
    print("post processing step 3/6")

    with open(f"{output_dir}/graphlet_id_mapper.csv", "r") as file:
        csv_reader = csv.reader(file, delimiter=":")
        for row in csv_reader:
            graphlet = ast.literal_eval(row[0])
            id = row[1].strip()
            three_node_graphlet_id[hash(graphlet)] = id
    print("post processing step 4/6")

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

    compare_files(
        Path(f"output_refactor/{species}/node_orbit.csv"),
        Path(f"final_output/{species}/node_orbit.csv"),
    )

    return (
        three_node_graphlet_count,
        three_node_graphlet_namespace,
        three_node_orbit_protein_data,
        three_node_orbit_namespace,
        three_node_graphlet_id,
        three_node_orbit_id,
        graphlet_config,
        node_orbit_arr,
    )


def get_stress_proteins(protein_id_dict, stress_path, delimiter):
    stress_proteins = []

    with open(stress_path, "r") as file:
        next(file)
        reader = csv.reader(file, delimiter=delimiter)
        for line in reader:
            if line[0] in protein_id_dict:
                stress_proteins.append(protein_id_dict[line[0]])
    return stress_proteins


def plot_stress_orbit_distribution(
    three_node_orbit_protein_data,
    three_node_orbit_id,
    stress_proteins_list,
    protein_id,
    species,
    output_dir,
    node_orbit_arr,
):

    observed_median_count = {}
    # print("Stress proteins : ", stress_proteins_list)
    stress_protein_orbit_dict = {}
    print("Counting stress proteins\n")
    for orbit in three_node_orbit_protein_data:
        # print(f"orbit {orbit}", end="\r")
        stress_protein_orbit_dict[orbit] = []
        for protein in stress_proteins_list:
            protein_orbit_count = node_orbit_arr[protein][
                int(three_node_orbit_id[orbit])
            ]
            stress_protein_orbit_dict[orbit] += [protein_orbit_count]

    print("Counting stress medians\n")
    for orbit in stress_protein_orbit_dict:
        # print(f"orbit {orbit}", end="\r")
        orbit_list = []
        if len(stress_protein_orbit_dict[orbit]) == 0:
            orbit_list = [0]
        else:
            orbit_list = stress_protein_orbit_dict[orbit]
        sorted_list = np.sort(orbit_list)
        median = np.median(sorted_list)
        observed_median_count[orbit] = median

    sample = 1000
    sample_size = len(stress_proteins_list)
    print("getting non stress proteins \n")
    non_stress_proteins = []
    for protein in protein_id:
        # print(f"protein {protein}", end="\r")
        if protein_id[protein] not in stress_proteins_list:
            non_stress_proteins.append(protein_id[protein])

    sample_results = {}

    for i in range(0, sample):
        print("sample :", i, end="\r")
        non_stress_sample = random.sample(non_stress_proteins, sample_size)
        array_stack = []
        for protein in non_stress_sample:
            array_stack.append(node_orbit_arr[protein])
        stacked = np.vstack(array_stack)
        orbit_medians_list = np.median(stacked, axis=0)

        for orbit in three_node_orbit_protein_data:
            orbit_index = three_node_orbit_id[orbit]
            if orbit not in sample_results:
                sample_results[orbit] = []
            sample_results[orbit] += [orbit_medians_list[int(orbit_index)]]

    significance = {}

    # Calculate if an orbit is significant
    for orbit in three_node_orbit_protein_data:
        count = 0
        # from the vector of all the random medians at a given orbit
        for random_median in sample_results[orbit]:
            if observed_median_count[orbit] > random_median:
                count += 1

        if count >= float(sample) * 0.99:
            significance[orbit] = 1
            print(three_node_orbit_id[orbit], ": significant")
        else:
            significance[orbit] = 0
            # print(three_node_orbit_id[orbit], ": NOT significant")

    with open(f"{output_dir}/stress_orbit_significance.csv", "w") as f:
        f.write(f"orbit_id\tsignificant?\tobserved_median\tvector_random_medians\n")
        for orbit in significance:
            f.write(
                f"{three_node_orbit_id[orbit]}\t{significance[orbit]}\t{observed_median_count[orbit]}\t{sample_results[orbit]}\n"
            )

    i = 0
    for orbit in three_node_orbit_protein_data:
        if significance[orbit] == 1:
            hist_data = sample_results[orbit]

            fig = plt.figure(figsize=(14, 6))
            plt.hist(hist_data)
            plt.axvline(x=int(observed_median_count[orbit]), color="r", linestyle="--")
            plt.title(
                f"{species} Random Median Samples at Orbit {int(three_node_orbit_id[orbit])} Distribution",
                fontsize=16,
            )
            plt.savefig(f"{output_dir}/sig_orbit/orbit{i}.pdf")
            plt.close()
        i += 1
    return None


def compare_files(file1, file2):
    with open(file1, "rb") as f1, open(file2, "rb") as f2:
        content1 = f1.read()
        content2 = f2.read()
        if content1 == content2:
            print("node orbit files are the same")
        else:
            print("DIFFERENT node orbit files")
        return content1 == content2


def species_wide_3_node_plots(n):
    species_list = ["bsub", "fly", "cerevisiae", "drerio", "elegans"]
    graphlet_counts = defaultdict(int)
    species_graphlet_counts = {}
    for species in species_list:
        with open(f"output_refactor/{species}/top_graphlet_counts.csv", "r") as file:
            csv_reader = csv.reader(file, delimiter="\t")
            next(csv_reader)
            for row in csv_reader:
                graphlet_id = row[1]
                count = int(row[2])
                graphlet_counts[graphlet_id] += count
                if species not in species_graphlet_counts:
                    species_graphlet_counts[species] = {}
                if graphlet_id not in species_graphlet_counts[species]:
                    species_graphlet_counts[species][graphlet_id] = count

    top_graphlets = sorted(graphlet_counts.items(), key=lambda x: x[1], reverse=True)[
        :n
    ]

    with open(
        f"output_refactor/species_wide/global_top_three_node_graphlet_counts.csv", "w"
    ) as f:
        f.write(f"graphlet_id\ttotal_count\tbsub\tdrerio\telegans\tfly\n")
        for graphlet, count in top_graphlets:
            f.write(
                f"{graphlet[0]}\t{species_graphlet_counts['bsub'][graphlet]}\t{species_graphlet_counts['drerio'][graphlet]}\t{species_graphlet_counts['elegans'][graphlet]}\t{species_graphlet_counts['fly'][graphlet]}\n"
            )
    return None


def species_wide_two_node_plots():
    species_list = ["bsub", "fly", "cerevisiae", "drerio", "elegans"]
    graphlet_ids = ["G_1", "G_2", "G_3", "G_4", "G_5"]
    species_graphlet_data = {}

    for species in species_list:
        species_graphlet_data[species] = {}
        with open(
            f"output_refactor/{species}/two_node_graphlet_counts.csv", "r"
        ) as file:
            csv_reader = csv.reader(file, delimiter=",")
            for row in csv_reader:
                species_graphlet_data[species][row[0]] = int(row[1].strip())

    graphlet_counts = {graphlet: [] for graphlet in graphlet_ids}

    for species_name in species_list:
        species_counts = species_graphlet_data[species_name]
        total_count = sum(species_counts[graphlet] for graphlet in graphlet_ids)
        for graphlet in graphlet_ids:
            percentage = (
                (species_counts[graphlet] / total_count) * 100 if total_count > 0 else 0
            )
            graphlet_counts[graphlet].append(percentage)

    index = np.arange(len(species_list))
    fig, ax = plt.subplots(figsize=(10, 6))
    bar_width = 0.5
    bottom_values = np.zeros(len(species_list))
    colors = ["blue", "green", "red", "yellow", "pink"]

    for i, graphlet in enumerate(graphlet_ids):
        bar = ax.bar(
            index,
            graphlet_counts[graphlet],
            bar_width,
            bottom=bottom_values,
            label=graphlet,
            color=colors[i],
        )
        bottom_values += np.array(graphlet_counts[graphlet])

    ax.set_xlabel("Species")
    ax.set_ylabel("Percentage of Graphlet Counts")
    ax.set_title("Stacked Bar Chart of Graphlet Percentages per Species")
    ax.set_xticks(index)
    ax.set_xticklabels(species_list)
    ax.legend()

    plt.tight_layout()
    plt.savefig(
        "output_refactor/species_wide/species_two_node_graphlet_dist_percentage.pdf"
    )
    # plt.show()
    plt.close()

    orbit_ids = ["1", "2", "3", "4", "5", "6", "7"]
    species_orbit_data = {}

    for species in species_list:
        species_orbit_data[species] = {}
        with open(f"output_refactor/{species}/two_node_orbit_counts.csv", "r") as file:
            csv_reader = csv.reader(file, delimiter=",")
            for row in csv_reader:
                species_orbit_data[species][row[0]] = int(row[1].strip())

    orbit1_count = []
    orbit2_count = []
    orbit3_count = []
    orbit4_count = []
    orbit5_count = []
    orbit6_count = []
    orbit7_count = []

    for species_name in species_list:
        species_counts = species_orbit_data[species_name]
        total_count = sum(species_counts.values())
        orbit1_count.append(species_counts["1"] / total_count * 100)
        orbit2_count.append(species_counts["2"] / total_count * 100)
        orbit3_count.append(species_counts["3"] / total_count * 100)
        orbit4_count.append(species_counts["4"] / total_count * 100)
        orbit5_count.append(species_counts["5"] / total_count * 100)
        orbit6_count.append(species_counts["6"] / total_count * 100)
        orbit7_count.append(species_counts["7"] / total_count * 100)

    orbit1_count = np.array(orbit1_count)
    orbit2_count = np.array(orbit2_count)
    orbit3_count = np.array(orbit3_count)
    orbit4_count = np.array(orbit4_count)
    orbit5_count = np.array(orbit5_count)
    orbit6_count = np.array(orbit6_count)
    orbit7_count = np.array(orbit7_count)

    index = np.arange(len(species_list))

    fig, ax = plt.subplots(figsize=(10, 6))
    bar_width = 0.5
    bar1 = ax.bar(index, orbit1_count, bar_width, label="1", color="blue")
    bar2 = ax.bar(
        index, orbit2_count, bar_width, bottom=orbit1_count, label="2", color="green"
    )
    bar3 = ax.bar(
        index,
        orbit3_count,
        bar_width,
        bottom=orbit1_count + orbit2_count,
        label="3",
        color="red",
    )
    bar4 = ax.bar(
        index,
        orbit4_count,
        bar_width,
        bottom=orbit1_count + orbit2_count + orbit3_count,
        label="4",
        color="yellow",
    )
    bar5 = ax.bar(
        index,
        orbit5_count,
        bar_width,
        bottom=orbit1_count + orbit2_count + orbit3_count + orbit4_count,
        label="5",
        color="pink",
    )
    bar6 = ax.bar(
        index,
        orbit6_count,
        bar_width,
        bottom=orbit1_count + orbit2_count + orbit3_count + orbit4_count + orbit5_count,
        label="6",
        color="orange",
    )
    bar7 = ax.bar(
        index,
        orbit7_count,
        bar_width,
        bottom=orbit1_count
        + orbit2_count
        + orbit3_count
        + orbit4_count
        + orbit5_count
        + orbit6_count,
        label="7",
        color="purple",
    )

    ax.set_xlabel("Species")
    ax.set_ylabel("Percentage of orbit ids")
    ax.set_title("Stacked Bar Chart of Orbit Percentages per Species")
    ax.set_xticks(index)
    ax.set_xticklabels(species_list)
    ax.legend()
    plt.tight_layout()
    plt.savefig(f"output_refactor/species_wide/species_two_node_orbit_dist.pdf")
    # plt.show()
    plt.close()

    return None


def main():
    print("running stats")

    species = "fly"
    output_dir = f"output_refactor/{species}"
    input_ppi = f"data/{species}_ppi.csv"
    input_reg = f"data/{species}_reg.csv"

    if species == "bsub":
        stress_dir = "data/oxidative_stress/txid224308/txid224308-stress-proteins.csv"
    if species == "drerio":
        stress_dir = "data/oxidative_stress/txid7955/txid7955-stress-proteins.csv"
    if species == "fly":
        stress_dir = "data/oxidative_stress/txid7227/txid7227-stress-proteins.csv"
    if species == "elegans":
        stress_dir = "data/oxidative_stress/txid6239/txid6239-stress-proteins.csv"
    if species == "cerevisiae":
        stress_dir = "data/oxidative_stress/txid559292/txid559292-stress-proteins.csv"

    protein_id, G, G_prime, graphlet_config = initialize_graphlet_data(
        input_ppi, input_reg
    )

    (
        three_node_graphlet_count,
        three_node_graphlet_namespace,
        three_node_orbit_protein_data,
        three_node_orbit_namespace,
        three_node_graphlet_id,
        three_node_orbit_id,
        graphlet_config,
        node_orbit_arr,
    ) = read_output_files(output_dir, species)

    plot_three_node_graphlet_distribution(
        three_node_graphlet_count,
        three_node_graphlet_namespace,
        three_node_graphlet_id,
        species,
        output_dir,
    )

    # significance orbit stats

    # stress_proteins_list = get_stress_proteins(protein_id, stress_dir, "\t")

    # significance = plot_stress_orbit_distribution(
    #     three_node_orbit_protein_data,
    #     three_node_orbit_id,
    #     stress_proteins_list,
    #     protein_id,
    #     species,
    #     output_dir,
    #     node_orbit_arr,
    # )

    # species wide stats
    species_wide_3_node_plots(10)

    # two node graphlet stats
    species_wide_two_node_plots()


if __name__ == "__main__":
    main()
