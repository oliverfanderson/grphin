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

    # Save the results to a CSV file
    output_file = f'{stress_dir}enrichment_results.csv'
    results_df.to_csv(output_file, index=False)

    # Print the first few rows of the results
    print(f"Results for {species}:")
    print(results_df.head())