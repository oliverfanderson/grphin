import csv
from pathlib import Path
import numpy as np


def main():
    print("")
    species = ["bsub", "fly", "cerevisiae", "drerio", "elegans"]
    txid = ["224308", "7227", "559292", "7955", "6239"]
    output_dir = Path(f"final_output/{species[2]}")
    stress_proteins_path = Path(f"data/oxidative_stress/txid{txid[2]}/txid{txid[2]}-stress-proteins.csv")
    protein_id_mapper_path = Path(f"final_output/{species[2]}/protein_id_mapper.csv")
    node_orbit_arr = np.loadtxt(
        f"{output_dir}/node_orbit.csv", delimiter=",", dtype=int
    )

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

    for protein in stress_protein_list:
        print(id_protein_dict[protein])

if __name__ == "__main__":
    main()