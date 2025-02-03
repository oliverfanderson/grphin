import csv
from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt


def main():
    print("")
    species = ["bsub", "fly", "cerevisiae", "drerio", "elegans"]
    name = ["B. subtilis", "D. melanogaster", "S. cerevisiae", "D. rerio", "C. elegans"]
    txid = ["224308", "7227", "559292", "7955", "6239"]
    output_dir = Path(f"final_output/{species[2]}")
    stress_proteins_path = Path(
        f"data/oxidative_stress/txid{txid[2]}/txid{txid[2]}-stress-proteins.csv"
    )
    protein_id_mapper_path = Path(f"final_output/{species[2]}/protein_id_mapper.csv")
    go_enrichment_path = Path(f"data/go_enrichment/{txid[2]}_mol.txt")
    node_orbit_arr = np.loadtxt(
        f"{output_dir}/node_orbit.csv", delimiter=",", dtype=int
    )
    orbit_id = 28

    protein_id_dict = {}
    id_protein_dict = {}
    with open(protein_id_mapper_path, "r") as f:
        csv_reader = csv.reader(f)
        for line in csv_reader:
            protein_id_dict[line[0]] = int(line[1].strip())
            id_protein_dict[int(line[1].strip())] = line[0]

    stress_protein_list = []
    with open(stress_proteins_path, "r") as f:
        csv_reader = csv.reader(f)
        next(f)
        for line in csv_reader:
            if protein_id_dict[line[0]] not in stress_protein_list:
                stress_protein_list.append(protein_id_dict[line[0]])

    print(stress_protein_list)
    rows, cols = node_orbit_arr.shape

    # for protein in stress_protein_list:
    #     print(protein)
    stress_protein_orbit_list = []
    for i in range(rows):
        if i - 1 in stress_protein_list and node_orbit_arr[i - 1][orbit_id] > 0:
            stress_protein_orbit_list.append(id_protein_dict[i])
    for protein in stress_protein_list:
        print(id_protein_dict[protein])
    print()
    print(",".join(stress_protein_orbit_list))

    labels = []
    counts = []
    total = 0
    with open(go_enrichment_path, "r") as f:
        csv_reader = csv.reader(f, delimiter="\t")
        for line in csv_reader:
            labels.append(line[1])
            counts.append(float(line[2]))
            total += int(line[2])
    sizes = [x / total for x in counts]

    # fig = plt.figure(figsize=(14, 6))
    # plt.pie(counts, labels=labels)
    # plt.legend(bbox_to_anchor=(0, 0.5), loc="center right")
    # plt.show()

    fig1, ax1 = plt.subplots()
    ax1.pie(sizes, labels=labels, autopct="%1.1f%%", startangle=140, pctdistance=0.8)
    ax1.axis("equal")
    plt.title(f"{name[2]} Stress Proteins at orbit-{orbit_id} GO Molecular Function Enrichment")
    plt.show()


if __name__ == "__main__":
    main()
