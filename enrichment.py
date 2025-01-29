import math
import pandas as pd
from scipy.stats import norm

species_list = ['bsub', 'cerevisiae', 'drerio', 'elegans', 'fly']

d_global_dict = {
    "bsub": 24150,
    "cerevisiae": 803494,
    "drerio": 141926,
    "elegans": 184276,
    "fly": 501168,
}

d_sub_dict = {
    "bsub": 112,
    "cerevisiae": 13910,
    "drerio": 198,
    "elegans": 5824,
    "fly": 1870,
}

def create_graphlet_presence_matrix():
    """
    Creates a binary presence/absence matrix of graphlets across species.

    Parameters:
    - species_list: List of species names
    - base_dir: Directory where the enrichment results are stored

    Returns:
    - A DataFrame where rows are graphlets, columns are species, and values are True/False.
    """
    graphlet_dict = {}

    # Load data and collect graphlets for each species
    for species in species_list:
        file_path = f"final_output/{species}_stress/enrichment_results.csv"
        
        try:
            df = pd.read_csv(file_path)
            graphlets = set(df['Graphlet'])  # Extract graphlets
            graphlet_dict[species] = graphlets
        except FileNotFoundError:
            print(f"Warning: File not found for {species}. Skipping...")
            graphlet_dict[species] = set()  # Empty set if file missing

    # Get all unique graphlets
    all_graphlets = sorted(set().union(*graphlet_dict.values()))

    # Create a DataFrame with species as columns and graphlets as index
    presence_matrix = pd.DataFrame(index=all_graphlets, columns=species_list, dtype=bool)

    # Populate the matrix with True/False values
    for species in species_list:
        presence_matrix[species] = presence_matrix.index.isin(graphlet_dict[species])

    return presence_matrix
def main():
    for species in species_list:
        stress_dir = f'final_output/{species}_stress/'
        stress_graphlet_counts = f'{stress_dir}graphlet_counts.csv'

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

        main_dir = f'final_output/{species}/'
        total_graphlet_counts = f'{main_dir}graphlet_counts.csv'

        # Open the total file and parse the data into a dictionary
        data_dict = {}
        with open(total_graphlet_counts, 'r') as file:
            for line in file:
                key, value = line.split(':')
                key = eval(key.strip())  # Convert the key into a Python tuple
                value = int(value.strip())  # Convert the value into an integer
                data_dict[key] = value

        # Convert the total dictionary to a DataFrame
        total = pd.DataFrame(list(data_dict.items()), columns=['Graphlet', 'Count'])

        # Prepare global and sub-network degree sums
        d_global = d_global_dict.get(species)  # Sum of degree of total network
        d_sub = d_sub_dict.get(species)  # Sum of degree of stress network

        # Create a list to store results
        results = []

        # Loop through each graphlet in the stress set
        for _, row in stress.iterrows():
            graphlet = row['Graphlet']
            m_sub = row['Count']  # Count in stress network

            # Get the equivalent count in the total set
            m_global = total.loc[total['Graphlet'] == graphlet, 'Count'].values
            if len(m_global) > 0:  # Check if the graphlet exists in the total set
                m_global = m_global[0]

                # Calculate enrichment and Z-score
                mu = m_global * (d_sub / d_global)
                enrichment = m_sub / mu if mu != 0 else float('inf')  # Avoid division by zero
                Z = (m_sub - mu) / math.sqrt(mu) if mu > 0 else 0  # Avoid negative or zero sqrt

                # Append the results
                results.append({
                    'Graphlet': graphlet,
                    # 'm_sub': m_sub,
                    # 'm_global': m_global,
                    # 'mu': mu,
                    'Enrichment': enrichment,
                    'Z': Z
                })
            else:
                # If graphlet is missing in the total set, handle appropriately
                results.append({
                    'Graphlet': graphlet,
                    # 'm_sub': m_sub,
                    # 'm_global': None,
                    # 'mu': None,
                    'Enrichment': None,
                    'Z': None
                })

        # Convert results to a DataFrame for easier analysis
        results_df = pd.DataFrame(results)
        results_df = results_df.sort_values(by='Enrichment', ascending=False)
        results_df['p_value'] = 2 * norm.sf(abs(results_df['Z']))

        sig_results_df = results_df[(results_df['p_value'] < 0.05) & (results_df['Enrichment'] > 1)]

        # Save the results to a CSV file
        output_file = f'{stress_dir}enrichment_results.csv'
        results_df.to_csv(output_file, index=False)

        # Print the number of results
        print(f"Number of results for {species}: {len(results_df)}")

        # Save the significant results to a CSV file
        sig_output_file = f'{stress_dir}enrichment_results_sig.csv'
        sig_results_df.to_csv(sig_output_file, index=False)

        # Print the number of significant results
        print(f"Number of significant results for {species}: {len(sig_results_df)}")

        # Print the first few rows of the results
        print(f"Results for {species}:")
        print(results_df.head())

    # Find graphlets across species
    graphlet_matrix = create_graphlet_presence_matrix()
    print(f'Graphlet presence matrix:\n{graphlet_matrix}')
    graphlet_matrix.to_csv("data/oxidative_stress/graphlet_presence_matrix.csv")

    # Get graphlets that are present in all species (all values in the row are True)
    common_graphlets = graphlet_matrix[graphlet_matrix.all(axis=1)]
    common_graphlets.to_csv("data/oxidative_stress/common_graphlets.csv")

    # Display the result
    print(f'Graphlets in all species:\n{common_graphlets}')


if __name__ == "__main__":
    main()