import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Load the standard data
standard_file = 'standard.csv'
standard_data = np.loadtxt(standard_file)

# Define the cutoff values
cutoff_values = [500, 1000, 2000, 5000, 10000]
output_files = [f'output_parallel_{cutoff}.csv' for cutoff in cutoff_values]

# Create a DataFrame to store relative differences
df_rel_diffs = pd.DataFrame()

for cutoff, file in zip(cutoff_values, output_files):
    output_data = np.loadtxt(file)
    # Compute relative differences
    relative_diff = (output_data - standard_data) / np.abs(standard_data)  # relative difference
    df_rel_diffs[f'rel_diff_cutoff_{cutoff}'] = relative_diff

# Remove outliers using IQR
def remove_outliers(df):
    Q1 = df.quantile(0.25)
    Q3 = df.quantile(0.75)
    IQR = Q3 - Q1
    return df[~((df < (Q1 - 1.5 * IQR)) | (df > (Q3 + 1.5 * IQR))).any(axis=1)]

# Remove outliers in the DataFrame
df_cleaned = remove_outliers(df_rel_diffs)

# Sample size for area plot
sample_size = 20000
x = np.arange(sample_size)

# Set colors for deeper tones for larger cutoffs
colors = ['lightblue', 'skyblue', 'cornflowerblue', 'royalblue', 'darkblue']

# Area plot of relative differences for the first 100 particles (after removing outliers)
plt.figure(figsize=(12, 8))
for i, cutoff in enumerate(cutoff_values):
    plt.fill_between(x, df_cleaned[f'rel_diff_cutoff_{cutoff}'][:sample_size], color=colors[i], alpha=0.6, label=f'Cutoff Radius {cutoff} (1e-10 m)')
plt.title('Relative Difference from Standard Net Force')
plt.xlabel('Particle Index')
plt.ylabel('Difference for Different Cutoff Radius')
plt.legend()
plt.grid(True)
plt.show()
