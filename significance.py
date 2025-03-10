import math
import pandas as pd

species_dict = {
    "txid6239" : "elegans",
    "txid7227" : "fly",
    "txid7955" : "drerio",
    "txid224308" : "bsub",
    "txid559292" : "cerevisiae"
}

def main():
    # List of taxon IDs to process
    taxon_ids = ["txid6239", "txid7227", "txid7955", "txid224308", "txid559292"]

    for txid in taxon_ids:
        stress_graphlet_counts = f"final_output/{species_dict[txid]}_stress/graphlet_counts.csv"

        # Open the stress file and parse the data into a dictionary
        data_dict = {}
        with open(stress_graphlet_counts, 'r') as file:
            for line in file:
                key, value = line.split(':')
                key = eval(key.strip())  # Convert the key into a Python tuple
                value = int(value.strip())  # Convert the value into an integer
                data_dict[key] = value

        # Convert the stress dictionary to a DataFrame
        stress = pd.DataFrame(list(data_dict.items()), columns=['Graphlet', 'Count'])
        print(stress.head())

if __name__ == "__main__":
    main()