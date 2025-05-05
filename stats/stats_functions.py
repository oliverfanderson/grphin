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
import requests
import seaborn as sns
import pandas as pd
from PyPDF2 import PdfMerger
from pyvis.network import Network
import networkx as nx
import webbrowser


sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from grphin import (
    initialize_graphlet_data,
    load_graphlet_config,
    initialize_three_node_graphlet_data,
)

PANTHER_URL = "https://pantherdb.org/services/oai/pantherdb/enrich/overrep"
species_txid = {
    "fly": "7227",
    "bsub": "224308",
    "drerio": "7955",
    "cerevisiae": "559292",
    "elegans": "6239",
}


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

    # plt.show()


def plot_three_node_graphlet_distribution(
    three_node_graphlet_count,
    three_node_graphlet_namespace,
    three_node_graphlet_id,
    species,
    output_dir,
):
    print("plotting graphlet distribution")
    hist_data = []
    x_label = [*range(1, len(three_node_graphlet_id) + 1, 1)]

    for graphlet in three_node_graphlet_id:
        if graphlet in three_node_graphlet_count:
            hist_data.append(three_node_graphlet_count[graphlet])
        else:
            hist_data.append(0)
    hist_data = [value if value > 0 else 0.1 for value in hist_data]

    n = 10
    top_n_graphlets = sorted(
        three_node_graphlet_count.items(), key=lambda item: item[1], reverse=True
    )[:n]
    top_n_graphlet_set = set([g[0] for g in top_n_graphlets])

    bar_colors = [
        "red" if graphlet in top_n_graphlet_set else "skyblue"
        for graphlet in three_node_graphlet_id
    ]

    fig = plt.figure(figsize=(14, 6))
    plt.bar(x_label, hist_data, color=bar_colors, edgecolor="black")
    plt.yscale("log")
    plt.title(f"{species} Graphlet Count Distribution", fontsize=16)
    plt.xlabel("Graphlet Index", fontsize=14)
    plt.ylabel("Count (log scale)", fontsize=14)
    plt.xticks(x_label[::2], fontsize=12)
    plt.grid(axis="y", linestyle="--", alpha=0.7)
    plt.tight_layout()
    plt.savefig(f"{output_dir}/{species}/three_node_graphlet_dist.pdf")
    plt.close()

    # Save top graphlet counts
    sorted_graphlet_dict = {
        key: value
        for key, value in sorted(
            three_node_graphlet_count.items(), key=lambda item: item[1], reverse=True
        )
    }

    with open(f"{output_dir}/{species}/top_graphlet_counts.csv", "w") as f:
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
    with open(f"{output_dir}/{species}/graphlet_hash_mapper.csv", "r") as file:
        csv_reader = csv.reader(file, delimiter=":")
        for row in csv_reader:
            graphlet_hash = int(row[0])
            graphlet = ast.literal_eval(row[1])
            three_node_graphlet_namespace[graphlet_hash] = graphlet
    print("post processing step 1/6")

    with open(f"{output_dir}/{species}/graphlet_counts.csv", "r") as file:
        csv_reader = csv.reader(file, delimiter=":")
        for row in csv_reader:
            graphlet = ast.literal_eval(row[0])
            count = int(row[1].strip())
            three_node_graphlet_count[hash(graphlet)] = count
    print("post processing step 2/6")

    with open(f"{output_dir}/{species}/orbit_hash_mapper.csv", "r") as file:
        csv_reader = csv.reader(file, delimiter=":")
        for row in csv_reader:
            orbit_hash = int(row[0])
            three_node_orbit_namespace[orbit_hash] = row[1]
    print("post processing step 3/6")

    with open(f"{output_dir}/{species}/graphlet_id_mapper.csv", "r") as file:
        csv_reader = csv.reader(file, delimiter=":")
        for row in csv_reader:
            graphlet = ast.literal_eval(row[0])
            id = row[1].strip()
            three_node_graphlet_id[hash(graphlet)] = id
    print("post processing step 4/6")

    orbit_id_dict = {}
    with open(f"{output_dir}/{species}/orbit_id_mapper.csv", "r") as file:
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
        f"{output_dir}/{species}/node_orbit.csv", delimiter=",", dtype=int
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
        Path(f"{output_dir}/{species}/node_orbit.csv"),
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


def analyze_stress_proteins(
    three_node_orbit_protein_data,
    three_node_orbit_id,
    stress_proteins_list,
    protein_id,
    species,
    output_dir,
    node_orbit_arr,
):

    id_to_protein = {}
    for protein, id in protein_id.items():
        id_to_protein[id] = protein

    orbit_stress_median_dict = {}
    orbit_stress_protein_dict = (
        {}
    )  # for each orbit, store counts of stress proteins i.e. orbit-1 = [12,4,5,0,3,2,6]

    print("Counting stress proteins\n")
    for orbit in three_node_orbit_protein_data:
        orbit_stress_protein_dict[orbit] = []
        for protein in stress_proteins_list:
            protein_orbit_count = node_orbit_arr[protein][
                int(three_node_orbit_id[orbit])
            ]
            orbit_stress_protein_dict[orbit] += [protein_orbit_count]

    print("Counting stress medians\n")
    for orbit in orbit_stress_protein_dict:
        orbit_list = []
        if len(orbit_stress_protein_dict[orbit]) == 0:
            orbit_list = [0]
        else:
            orbit_list = orbit_stress_protein_dict[orbit]
        sorted_list = np.sort(
            orbit_list
        )  # loses the order of which protein in the list
        median = np.median(sorted_list)
        orbit_stress_median_dict[orbit] = (
            median  # for each orbit, store the median counts for all stress proteins
        )
        # TODO: why do we do median instead of mean?

    sample = 1000
    sample_size = len(
        stress_proteins_list
    )  # TODO: the size will effect the significance right? more vs less stress proteins

    print("getting non stress proteins \n")
    non_stress_proteins = []
    for protein in protein_id:
        if protein_id[protein] not in stress_proteins_list:
            non_stress_proteins.append(protein_id[protein])

    sample_results = {}
    # dict where the keys are orbits
    # the values are lists of all sampled non stress protein median counts at each orbit

    for i in range(0, sample):
        print("sample :", i, end="\r")
        non_stress_sample = random.sample(
            non_stress_proteins, sample_size
        )  # this can affect our analysis right?

        array_stack = []
        # for each non-stress protein, get all the orbit counts as a list. then store this in array_stack
        # non-stress-protein-1 = [0,21,2,....2]
        # non-stress-protein-2 = [1,2,242,....21]
        # non-stress-protein-3 = [3,1,25,....231]
        # array_stack = [non-stress-protein-1, non-stress-protein-4, non-stress-protein-3]

        for protein in non_stress_sample:
            array_stack.append(node_orbit_arr[protein])

        # since all elements in array_stack have the same length,
        # we can stack them and find the median for each index (representing an orbit)
        stacked = np.vstack(array_stack)
        orbit_medians_list = np.median(stacked, axis=0)

        # append each sampled orbit medians to sample results
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
            if orbit_stress_median_dict[orbit] > random_median:
                count += 1
        if count >= float(sample) * 0.99:
            significance[orbit] = 1
            print(three_node_orbit_id[orbit], ": significant")
        else:
            significance[orbit] = 0

    with open(f"{output_dir}/{species}/stress_orbit_significance.csv", "w") as f:
        f.write(f"orbit_id\tsignificant?\tobserved_median\tvector_random_medians\n")
        for orbit in significance:
            f.write(
                f"{int(three_node_orbit_id[orbit]) + 1}\t{significance[orbit]}\t{orbit_stress_median_dict[orbit]}\t{sample_results[orbit]}\n"
            )

    # plot each significant orbits' p value with observed stress proteins' median count at the orbit
    i = 0
    for orbit in three_node_orbit_protein_data:
        if significance[orbit] == 1:
            hist_data = sample_results[orbit]

            fig = plt.figure(figsize=(14, 6))
            plt.hist(hist_data)
            plt.axvline(
                x=int(orbit_stress_median_dict[orbit]), color="r", linestyle="--"
            )
            plt.title(
                f"{species} Random Median Samples at Orbit {int(three_node_orbit_id[orbit])} Distribution",
                fontsize=16,
            )
            plt.savefig(f"{output_dir}/{species}/sig_orbit/per_orbit/orbit{i + 1}.pdf")
            plt.close()
        i += 1

    # find all the stress proteins that are significant orbits
    significant_orbit_list = [
        orbit for orbit in significance if significance[orbit] == 1
    ]

    sig_orbit_stress_protein_dict = {}

    for protein in stress_proteins_list:
        for orbit in significant_orbit_list:
            if node_orbit_arr[protein][int(three_node_orbit_id[orbit])] > 0:
                if orbit not in sig_orbit_stress_protein_dict:
                    sig_orbit_stress_protein_dict[orbit] = []
                sig_orbit_stress_protein_dict[orbit] += [
                    (protein, node_orbit_arr[protein][int(three_node_orbit_id[orbit])])
                ]

    # do go enrichment analysis
    # go_class = "GO:0008150" # biological process
    # go_class = "GO:0003674" # molecular function
    # go_class = "GO:0005575" # Cellular Component

    fdr_dist_data = []
    go_class_list = ["GO:0008150", "GO:0003674", "GO:0005575"]
    # go_class_list = ["GO:0008150"]
    # stress_protein_data = []
    protein_size_data = []
    orbit_data = []
    go_data = []
    for go_class in go_class_list:
        with open(
            f"{output_dir}/{species}/sig_orbit/go_enrichment_{go_class}.txt", "w+"
        ) as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerows(
                [
                    [
                        "orbit",
                        "term_id",
                        "term_label",
                        "pValue",
                        "fdr",
                        "fold_enrichment",
                        "expected",
                        "number_in_list",
                        "number_in_reference",
                        "plus_minus",
                    ]
                ]
            )
            print(go_class)
            for orbit in sig_orbit_stress_protein_dict:
                print()
                gene_list = []

                # filter the number of stress proteins
                # default
                for protein_count_tuple in sig_orbit_stress_protein_dict[orbit]:
                    protein = id_to_protein[protein_count_tuple[0]]
                    # print(protein)
                    count = protein_count_tuple[1]
                    gene_list.append(protein)

                # # filter based on only top 25th percentile of counts
                # sorted_protein_count_list = sorted(sig_orbit_stress_protein_dict[orbit], key=lambda x:[1])
                # sorted_count_list = [x[1] for x in sorted_protein_count_list]
                # percentile = np.percentile(sorted_count_list, 75)
                # for protein_count_tuple in sorted_protein_count_list:
                #     protein = id_to_protein[protein_count_tuple[0]]
                #     count = protein_count_tuple[1]
                #     if count >= percentile:
                #         gene_list.append(protein)
                # if go_class == "GO:0008150":
                #     protein_size_data.append(len(gene_list))
                #     orbit_data.append(f"{int(three_node_orbit_id[orbit])}")
                # # stress_protein_data.append({"gene_size": len(gene_list), "orbit": int(three_node_orbit_id[orbit]), "go_class": go_class, "species" :species})
                # gene_list = format_gene_list(gene_list)
                # print(f"{three_node_orbit_id[orbit]} - {len(gene_list)}")

                gene_list = (",".join(gene_list),)
                query_params = build_query_params(
                    gene_list,
                    species_txid[species],
                    go_class,
                    "FISHER",
                    "FDR",
                    "NONE",
                )
                results = call_panther_api(query_params)
                parsed_results = parse_results(results)

                for entry in parsed_results:
                    fdr_dist_data.append(
                        {"fdr": entry["fdr"], "species": species, "go_class": go_class}
                    )

                top_hits = filter_and_sort_results(parsed_results, fdr_threshold=0.01)
                print(
                    f"orbit-{three_node_orbit_id[orbit]}significant p-value hits = {len(top_hits)}"
                )

                for hit in top_hits:
                    writer.writerows(
                        [
                            [
                                str(int(three_node_orbit_id[orbit]) + 1),
                                hit["term_id"],
                                hit["term_label"],
                                hit["pValue"],
                                hit["fdr"],
                                hit["fold_enrichment"],
                                hit["expected"],
                                hit["number_in_list"],
                                hit["number_in_reference"],
                                hit["plus_minus"],
                            ]
                        ]
                    )
            f.close()

    with open(f"{output_dir}/{species}/sig_orbit/summary.txt", "w+") as f:
        headers = ["orbit", "protein", "count"]
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(headers)
        for orbit in sig_orbit_stress_protein_dict:
            for protein_count_tuple in sig_orbit_stress_protein_dict[orbit]:
                protein = protein_count_tuple[0]
                count = protein_count_tuple[1]
                row = [three_node_orbit_id[orbit], id_to_protein[protein], count]
                writer.writerow(
                    [int(three_node_orbit_id[orbit]) + 1, id_to_protein[protein], count]
                )
        f.close()

    # print(stress_protein_data["orbit"])
    # print(type(stress_protein_data["orbit"]))

    plt.figure(figsize=(12, 8))
    plt.title(f"{species} GO enrichment gene set size per orbit")
    plt.ylabel("Gene set size")
    plt.xlabel("Orbit")
    plt.bar(orbit_data, protein_size_data, color="skyblue")
    plt.savefig(f"{output_dir}/{species}/sig_orbit/{species}_orbit_gene_set_plt.pdf")
    # plt.show()
    plt.close()
    return fdr_dist_data


def compare_files(file1, file2):
    with open(file1, "rb") as f1, open(file2, "rb") as f2:
        content1 = f1.read()
        content2 = f2.read()
        if content1 == content2:
            print("node orbit files are the same")
        else:
            print("DIFFERENT node orbit files")
        return content1 == content2


def species_wide_3_node_plots(n, output_dir):
    species_list = ["bsub", "fly", "cerevisiae", "drerio", "elegans"]
    graphlet_counts = defaultdict(int)
    species_graphlet_counts = {}
    for species in species_list:
        with open(f"{output_dir}/{species}/top_graphlet_counts.csv", "r") as file:
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
        f"{output_dir}/species_wide/global_top_three_node_graphlet_counts.csv", "w"
    ) as f:
        f.write(f"graphlet_id\ttotal_count\tbsub\tdrerio\telegans\tfly\n")
        for graphlet, count in top_graphlets:
            f.write(
                f"{graphlet[0]}\t{species_graphlet_counts['bsub'][graphlet]}\t{species_graphlet_counts['drerio'][graphlet]}\t{species_graphlet_counts['elegans'][graphlet]}\t{species_graphlet_counts['fly'][graphlet]}\n"
            )
    return None


def species_wide_two_node_plots(output_dir):
    species_list = ["bsub", "fly", "cerevisiae", "drerio", "elegans"]
    graphlet_ids = ["G_1", "G_2", "G_3", "G_4", "G_5"]
    species_graphlet_data = {}

    for species in species_list:
        species_graphlet_data[species] = {}
        with open(f"{output_dir}/{species}/two_node_graphlet_counts.csv", "r") as file:
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
        f"{output_dir}/species_wide/species_two_node_graphlet_dist_percentage.pdf"
    )
    # plt.show()
    plt.close()

    orbit_ids = ["1", "2", "3", "4", "5", "6", "7"]
    species_orbit_data = {}

    for species in species_list:
        species_orbit_data[species] = {}
        with open(f"{output_dir}/{species}/two_node_orbit_counts.csv", "r") as file:
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
    plt.savefig(f"{output_dir}/species_wide/species_two_node_orbit_dist.pdf")
    # plt.show()
    plt.close()

    return None


def format_gene_list(gene_list):
    return ",".join(gene_list)


def build_query_params(
    gene_list, txid, go_class, enrichment_test, correction, mapped_info
):

    return {
        "geneInputList": format_gene_list(gene_list),
        "organism": txid,
        "annotDataSet": go_class,
        "enrichmentTestType": enrichment_test,
        "correction": correction,
        "mappedInfo": mapped_info,
    }


def call_panther_api(query_params):
    response = requests.get(PANTHER_URL, params=query_params)
    if response.ok:
        return response.json()
    else:
        raise RuntimeError(f"API failed {response.status_code} - {response.text}")


def parse_results(data):
    results = data.get("results", {}).get("result", [])
    parsed = []

    for entry in results:
        result_entry = {
            "term_id": entry.get("term", {}).get("id", "NA"),
            "term_label": entry.get("term", {}).get("label", "NA"),
            "pValue": entry.get("pValue", "NA"),
            "fdr": entry.get("fdr", "NA"),
            "fold_enrichment": entry.get("fold_enrichment", "NA"),
            "expected": entry.get("expected", "NA"),
            "number_in_list": entry.get("number_in_list", "NA"),
            "number_in_reference": entry.get("number_in_reference", "NA"),
            "plus_minus": entry.get("plus_minus", "NA"),
        }
        parsed.append(result_entry)
    return parsed


def print_results(parsed_results):
    for res in parsed_results:
        print(
            f"{res['term']} - p-value: {res['p-value']}, FDR: {res['fdr']}, "
            f"Fold Enrichment: {res['fold_enrichment']}, Count: {res['number_in_list']}"
        )


def filter_and_sort_results(results, fdr_threshold=0.05):
    filtered = []
    for r in results:
        try:
            fdr = float(r["fdr"])
            if fdr <= fdr_threshold:
                filtered.append(r)
                # print(r)
        except (ValueError, TypeError):
            continue
    # Sort by FDR ascending
    sorted_filtered = sorted(filtered, key=lambda x: float(x["fdr"]))
    return sorted_filtered


def merge_pdfs(output_dir, species_list):
    merger = PdfMerger()
    for species in species_list:
        pdf_path = f"{output_dir}/{species}/sig_orbit/{species}_orbit_gene_set_plt.pdf"
        if os.path.exists(pdf_path):
            merger.append(pdf_path)
    merger.write(f"{output_dir}/species_wide/orbit_gene_set_plt.pdf")

    merger = PdfMerger()
    for species in species_list:
        pdf_path = f"{output_dir}/{species}/three_node_graphlet_dist.pdf"
        if os.path.exists(pdf_path):
            merger.append(pdf_path)
    merger.write(f"{output_dir}/species_wide/three_node_graphlet_dist.pdf")


def get_go_network(path):
    G = nx.DiGraph()
    with open(path, "r") as f:
        csvreader = csv.reader(f, delimiter="\t")
        next(csvreader)
        i = 0
        for row in csvreader:
            g1 = row[0]
            g2 = row[1]

            if not G.has_node(g1):
                G.add_node(g1)
            if not G.has_node(g2):
                G.add_node(g2)
            if not G.has_edge(g1, g2):
                G.add_edge(g2, g1)
            # if i == 10:
            #     break
            i += 1

    return G


def analyze_go_enrichment(species_list):
    mixed_orbits_list = get_mixed_orbit_list()

    go_enrichment_files_list = [
        "go_enrichment_GO:0003674.txt",
        "go_enrichment_GO:0005575.txt",
        "go_enrichment_GO:0008150.txt",
    ]

    output_dir = "stats/output"
    go_path = Path("data/go_hierarchy/is_a_import_2024-07-17.tsv")
    G = get_go_network(go_path)
    go_id_to_name = {}
    results = {}

    mixed_orbits_summary = {}
    node_attributes_dict = {}

    for species in species_list:
        for file in go_enrichment_files_list:
            go_type = file.split("_")[-1].split(".")[0]
            file_path = Path(f"{output_dir}/{species}/sig_orbit/{file}")
            node_list = []
            with open(file_path, "r") as f:
                csvreader = csv.reader(f, delimiter="\t")
                next(csvreader)
                for row in csvreader:
                    orbit = row[0]
                    if int(orbit) in mixed_orbits_list:
                        go_id = row[1]
                        go_name = row[2]

                        if go_id not in go_id_to_name:
                            go_id_to_name[go_id] = go_name

                        key = f"{species}_{go_type}_{orbit}"
                        if key not in results:
                            results[key] = []
                            node_attributes_dict[key] = {}
                        results[key].append(go_id)
                        node_attributes_dict[key][go_id] = row

                        if f"{species}_{go_type}" not in mixed_orbits_summary:
                            mixed_orbits_summary[f"{species}_{go_type}"] = []
                        mixed_orbits_summary[f"{species}_{go_type}"].append(row)

    # for key in mixed_orbits_summary:
    #     # print(key)
    #     species = key.split("_")[0]
    #     go_type = key.split("_")[1]
    #     with open(f"{output_dir}/{species}/sig_orbit/mixed_{go_type}", "w+") as f:
    #         writer = csv.writer(f, delimiter="\t")
    #         writer.writerows(
    #             [
    #                 [
    #                     "orbit",
    #                     "term_id",
    #                     "term_label",
    #                     "pValue",
    #                     "fdr",
    #                     "fold_enrichment",
    #                     "expected",
    #                     "number_in_list",
    #                     "number_in_reference",
    #                     "plus_minus",
    #                 ]
    #             ]
    #         )
    #         for row in mixed_orbits_summary[key]:
    #             # print(row)
    #             writer.writerows([row])

    for key in results:
        filename = "_".join(key.split("_")[1:3])
        species = key.split("_")[0]
        go_type = filename.split("_")[0]
        orbit = key.split("_")[-1]
        output_path = f"{output_dir}/{species}/sig_orbit/vis/mixed_{filename}.html"
        H = nx.induced_subgraph(G, results[key])

        nt = Network("1000px", "1000px")

        # for edge in H.edges():
        #     nt.add_edge(edge[0], edge[1], arrows="to")

        go_attributes = node_attributes_dict[f"{species}_{go_type}_{orbit}"]
        nt.from_nx(H)

        write_type = ""
        if os.path.exists(f"{output_dir}/{species}/sig_orbit/mixed_{go_type}"):
            write_type = "a"
        else:
            write_type = "w"
        with open(f"{output_dir}/{species}/sig_orbit/mixed_{go_type}", write_type) as f:
            writer = csv.writer(f, delimiter="\t")
            for node in nt.nodes:
                go_id = node["id"]
                go_name = go_id_to_name.get(go_id, go_id)
                node["label"] = go_name
                node["orbit"] = go_attributes[node["id"]][0]
                node["pValue"] = go_attributes[node["id"]][3]
                node["fdr"] = go_attributes[node["id"]][4]
                node["fold_enrichment"] = go_attributes[node["id"]][5]
                node["expected"] = go_attributes[node["id"]][6]
                node["number_in_list"] = go_attributes[node["id"]][7]
                node["number_in_reference"] = go_attributes[node["id"]][8]
                node["plus_minus"] = go_attributes[node["id"]][9]

                writer.writerows(
                    [
                        [
                            int(node["orbit"]),
                            go_id,
                            node["label"],
                            node["pValue"],
                            node["fdr"],
                            node["fold_enrichment"],
                            node["expected"],
                            node["number_in_list"],
                            node["number_in_reference"],
                            node["plus_minus"],
                        ]
                    ]
                )

        for edge in nt.edges:
            edge["arrows"] = "from"

        nt.show(output_path, notebook=False)
        abs_path = os.path.abspath(output_path)
        # print(f"Visualization saved to: {abs_path}")
        # # webbrowser.open(f"file://{abs_path}")
        print()


def get_mixed_orbit_list():
    result = []
    result.extend(list(range(2, 62)))
    result.extend(list(range(79, 113)))
    result.extend(list(range(113, 187)))
    result.extend(list(range(187, 244)))

    return result


def main():
    print("running stats")
    species_list = ["bsub", "drerio", "fly", "elegans", "cerevisiae"]
    # species_list = ["fly"]

    output_dir = "stats/output"

    go_enrichment_data = False

    if not go_enrichment_data:
        fdr_dist_data = []

        for species in species_list:
            print(species)
            input_ppi = f"data/{species}_ppi.csv"
            input_reg = f"data/{species}_reg.csv"
            stress_dir = f"data/oxidative_stress/txid{species_txid[species]}/txid{species_txid[species]}-stress-proteins.csv"
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

            stress_proteins_list = get_stress_proteins(protein_id, stress_dir, "\t")

            fdr_dist_data.extend(
                analyze_stress_proteins(
                    three_node_orbit_protein_data,
                    three_node_orbit_id,
                    stress_proteins_list,
                    protein_id,
                    species,
                    output_dir,
                    node_orbit_arr,
                )
            )

        merge_pdfs(output_dir, species_list)

    analyze_go_enrichment(species_list)

    # species wide stats
    species_wide_3_node_plots(10, output_dir)

    # two node graphlet stats
    species_wide_two_node_plots(output_dir)


if __name__ == "__main__":
    main()
