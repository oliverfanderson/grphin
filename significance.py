import math
import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument(
    "-i", "--iterations", required=True, help="Number of iterations to run"
)
args = parser.parse_args()

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
    # taxon_ids = ["txid224308"]

    iterations = int(args.iterations)  # Number of randomized networks to compare

    for txid in taxon_ids:
        stress_graphlet_counts = f"final_output/{species_dict[txid]}_stress/graphlet_counts.csv"
        output_file = f"data/oxidative_stress/{txid}/graphlet_significance.csv"

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
        # print(stress.head())

        # Initialize the significance DataFrame with the same columns as the stress DataFrame, and set the tally column to zero
        significance = stress.copy()
        significance['Tally'] = 0

        for i in range(iterations):

            random_graphlet_counts = f"data/oxidative_stress/{txid}/randomized_networks/graphlet_counts{i}.csv"

            # Open the stress file and parse the data into a dictionary
            data_dict = {}
            with open(random_graphlet_counts, 'r') as file:
                for line in file:
                    key, value = line.split(':')
                    key = eval(key.strip())  # Convert the key into a Python tuple
                    value = int(value.strip())  # Convert the value into an integer
                    data_dict[key] = value

            # Convert the stress dictionary to a DataFrame
            random = pd.DataFrame(list(data_dict.items()), columns=['Graphlet', 'Count'])
            # print(random.head())

            # Count the number of times the graphlets appear in the randomized networks more than the stress network
            for _, stress_row in stress.iterrows():
                for _, random_row in random.iterrows():
                    if stress_row['Graphlet'] == random_row['Graphlet']:
                        if stress_row['Count'] <= random_row['Count']:
                            significance.loc[significance['Graphlet'] == stress_row['Graphlet'], 'Tally'] += 1

        # Calculate the significance by dividing the tally by the number of iterations
        significance['Significance'] = significance['Tally'] / iterations
        # significance.drop('Tally', axis=1, inplace=True)

        # Save the significance DataFrame to a CSV file
        significance.to_csv(output_file, sep='\t', index=False)
        # print(significance.head())
        print(f"Graphlet significance written to {output_file}")

if __name__ == "__main__":
    main()