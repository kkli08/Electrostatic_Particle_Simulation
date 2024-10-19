import os
import re
import numpy as np
import matplotlib.pyplot as plt

# Directory containing the timing txt files (update this if necessary)
directory = '.'  # Set this to the appropriate directory path

# Regex pattern to match the total time and per-thread time
total_time_pattern = re.compile(r"Total computation time: ([\d.]+) seconds")
thread_time_pattern = re.compile(r"Thread handling particles from (\d+) to (\d+) took ([\d.]+) seconds")

# Store results
threads = []
total_times = []
average_particle_times = []

# Check if the directory exists
if not os.path.exists(directory):
    print(f"Directory {directory} does not exist.")
else:
    # Iterate over all files in the directory
    for filename in os.listdir(directory):
        if filename.startswith("timing_cutoff_") and filename.endswith(".txt"):
            print(f"Processing file: {filename}")  # Debugging: print each file that is being processed

            # Extract the number of threads from the filename
            num_threads_match = re.search(r'threads_(\d+)', filename)
            if num_threads_match:
                num_threads = int(num_threads_match.group(1))
                threads.append(num_threads)

                # Initialize variables for total and thread times
                total_time = 0
                thread_particle_times = []

                # Open and read the file
                with open(os.path.join(directory, filename), 'r') as file:
                    for line in file:
                        # Check if the line matches the total computation time
                        total_time_match = total_time_pattern.search(line)
                        if total_time_match:
                            total_time = float(total_time_match.group(1))

                        # Check if the line matches the per-thread time
                        thread_time_match = thread_time_pattern.search(line)
                        if thread_time_match:
                            start_idx = int(thread_time_match.group(1))
                            end_idx = int(thread_time_match.group(2))
                            thread_time = float(thread_time_match.group(3))

                            # Calculate number of particles handled by the thread
                            num_particles_handled = end_idx - start_idx

                            # Calculate average time per particle for this thread
                            if num_particles_handled > 0:
                                avg_time_per_particle = thread_time / num_particles_handled
                                thread_particle_times.append(avg_time_per_particle)

                # Calculate the overall average time per particle across all threads
                average_particle_time = np.mean(thread_particle_times) if thread_particle_times else 0

                # Store the results
                total_times.append(total_time)
                average_particle_times.append(average_particle_time)

    # Check if we found any matching files
    if not threads:
        print("No matching txt files found.")
    else:
        # Sort the results by the number of threads (in case they are unordered)
        sorted_indices = np.argsort(threads)
        threads = np.array(threads)[sorted_indices]
        total_times = np.array(total_times)[sorted_indices]
        average_particle_times = np.array(average_particle_times)[sorted_indices]

        # Plot 1: Total Time Consumed vs Number of Threads
        plt.figure(figsize=(10, 6))
        plt.plot(threads, total_times, color='skyblue', marker='o', linestyle='-', linewidth=2, markersize=8)
        plt.title('Total Time Consumed vs Number of Threads')
        plt.xlabel('Number of Threads')
        plt.ylabel('Total Time Consumed (seconds)')
        plt.grid(True)
        plt.xticks(threads)
        plt.show()

        # Plot 2: Average Time per Particle vs Number of Threads
        plt.figure(figsize=(10, 6))
        plt.plot(threads, average_particle_times, color='skyblue', marker='o', linestyle='-', linewidth=2, markersize=8)
        plt.title('Average Time per Particle vs Number of Threads')
        plt.xlabel('Number of Threads')
        plt.ylabel('Average Time per Particle (seconds)')
        plt.grid(True)
        plt.xticks(threads)
        plt.show()
