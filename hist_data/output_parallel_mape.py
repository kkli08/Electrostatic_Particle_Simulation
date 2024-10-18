import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Load the standard data
standard_file = 'standard.csv'
standard_data = np.loadtxt(standard_file)

# Define the cutoff values
cutoff_values = [500, 1000, 2000, 5000, 10000, 15000, 20000, 30000, 40000, 50000]
output_files = [f'output_parallel_{cutoff}.csv' for cutoff in cutoff_values]

# Create a list to store the percentage errors
percent_errors = []

# Calculate the percentage error for each cutoff radius
for cutoff, file in zip(cutoff_values, output_files):
    output_data = np.loadtxt(file)

    # Avoid division by zero: replace zeros in standard_data with a small value (e.g., 1e-10)
    standard_data_no_zero = np.where(standard_data == 0, 1e-10, standard_data)

    # Compute the absolute percentage error
    abs_percentage_error = np.abs((output_data - standard_data) / np.abs(standard_data_no_zero)) * 100

    # Calculate the mean percentage error, ignoring NaN values
    mean_percentage_error = np.nanmean(abs_percentage_error)
    percent_errors.append(mean_percentage_error)

# Print the calculated mean percentage errors
print("Cutoff values:", cutoff_values)
print("Mean percentage errors:", percent_errors)

# Plot the percentage errors against the cutoff values
plt.figure(figsize=(10, 6))
plt.plot(np.arange(len(cutoff_values)), percent_errors, color='skyblue', marker='o', linestyle='-', linewidth=2, markersize=8)

# Set x-axis labels to show the cutoff values equally spaced
plt.xticks(np.arange(len(cutoff_values)), labels=cutoff_values)

# Annotate each point with its percentage error value, slightly shifted to the right
for i, (x, y) in enumerate(zip(np.arange(len(cutoff_values)), percent_errors)):
    plt.text(x + 0.1, y + 0.5, f'{y:.2f}%', ha='left', fontsize=10)  # Shift slightly right with x + 0.1

plt.title('Mean Percentage Error vs Cutoff Radius')
plt.xlabel('Cutoff Radius (1e-10 m)')
plt.ylabel('Mean Percentage Error (%)')
plt.grid(True)

# Adjust the y-axis limits for better visibility
plt.ylim(min(percent_errors) * 0.9, max(percent_errors) * 1.1)  # Add some padding to y-axis limits

plt.show()
