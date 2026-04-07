import pandas as pd
import numpy as np
import sys

# Load the existing data file
input_filename = sys.argv[1] #'simulated_gwas_data.tsv'
df = pd.read_csv(input_filename, sep='\t')

# Function to generate new, realistic floating-point values
def simulate_new_values(df):
    """
    Simulates new random values for the float columns in a DataFrame.
    """
    num_rows = len(df)

    # Simulate 'beta' values from a normal distribution centered around 0
    # GWAS betas are typically small, so we use a small standard deviation
    df['beta'] = np.random.normal(loc=0, scale=0.05, size=num_rows)

    # Simulate 'standard_error' values from a uniform distribution
    # SEs are typically small and positive
    df['standard_error'] = np.random.uniform(low=0.01, high=0.1, size=num_rows)

    # Simulate 'effect_allele_frequency' (EAF) from a uniform distribution
    # EAF ranges from 0 to 1
    df['effect_allele_frequency'] = np.random.uniform(low=0.01, high=0.99, size=num_rows)

    # Simulate 'p_value' values from a uniform distribution
    # P-values are uniformly distributed under the null hypothesis (no association)
    df['p_value'] = np.random.uniform(low=0, high=1, size=num_rows)

    # Ensure a few variants have very low p-values to simulate true signals
    num_signals = int(num_rows * 0.05)  # e.g., 5% of variants are significant
    if num_signals > 0:
        significant_indices = np.random.choice(df.index, size=num_signals, replace=False)
        df.loc[significant_indices, 'p_value'] = np.random.uniform(low=1e-10, high=5e-8, size=num_signals)

    return df

# Simulate new values for the DataFrame
df_simulated = simulate_new_values(df.copy())

# Save the updated DataFrame to a new tab-separated file
output_filename = sys.argv[2] #'simulated_gwas_data_new_floats.tsv'
df_simulated.to_csv(output_filename, sep='\t', index=False)

print(f"New simulated dataset with updated floats created: {output_filename}")
print("\nFirst 5 rows of the new dataset:")
print(df_simulated.head())