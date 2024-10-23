import os
import re
import matplotlib.pyplot as plt

# Function to read average leader computation time from timing_summary.txt files
def read_average_leader_times():
    leader_numbers = []
    average_times = []

    for i in range(1, 8):  # Leaders from 1 to 7
        folder = f"{i}"
        summary_file = os.path.join(folder, "timing_summary.txt")

        if os.path.exists(summary_file):
            with open(summary_file, 'r') as f:
                content = f.read()
                match = re.search(r"Average leader computation time: ([\d\.]+) seconds\.", content)
                if match:
                    avg_time = float(match.group(1))
                    leader_numbers.append(i)
                    average_times.append(avg_time)
                    print(f"Read average time {avg_time} seconds for {i} leader(s).")
                else:
                    print(f"Average leader computation time not found in {summary_file}.")
        else:
            print(f"{summary_file} does not exist.")

    return leader_numbers, average_times

# Function to read total computation times from mode 2 files
def read_mode2_times():
    mode2_folder = "../mode_2_threads_compare"
    thread_numbers = []
    total_times = []

    for threads in [10, 20, 30, 40, 50, 60, 70]:
        filename = f"timing_cutoff_0.000005_threads_{threads}.txt"
        filepath = os.path.join(mode2_folder, filename)

        if os.path.exists(filepath):
            with open(filepath, 'r') as f:
                content = f.read()
                match = re.search(r"Total computation time: ([\d\.]+) seconds\.", content)
                if match:
                    total_time = float(match.group(1))
                    thread_numbers.append(threads)
                    total_times.append(total_time)
                    print(f"Read total time {total_time} seconds for {threads} threads in Mode 2.")
                else:
                    print(f"Total computation time not found in {filepath}.")
        else:
            print(f"{filepath} does not exist.")

    return thread_numbers, total_times

# Function to read total computation times from mode 3 files
def read_mode3_times():
    thread_numbers = []
    total_times = []

    for leaders in range(1, 8):  # Leaders from 1 to 7
        folder = f"{leaders}"
        leader_times = []
        num_threads_per_leader = 0

        for leader_id in range(leaders):
            filename = f"timing_cutoff_0.000005_threads_10_leader_{leader_id}.txt"
            filepath = os.path.join(folder, filename)

            if os.path.exists(filepath):
                with open(filepath, 'r') as f:
                    content = f.read()

                    # Get total computation time for this leader
                    match = re.search(r"Leader \d+ total computation time: ([\d\.]+) seconds\.", content)
                    if match:
                        leader_time = float(match.group(1))
                        leader_times.append(leader_time)
                        print(f"Read leader {leader_id} time {leader_time} seconds in folder {folder}.")
                    else:
                        print(f"Leader total computation time not found in {filepath}.")

                    # Count threads
                    num_threads = len(re.findall(r"Thread 0x[\da-f]+ total computation time:", content))
                    num_threads_per_leader = num_threads  # Assuming same for all leaders

            else:
                print(f"{filepath} does not exist.")

        if leader_times:
            # Total computation time is the maximum of leader times (since they run in parallel)
            total_time = max(leader_times)
            total_times.append(total_time)

            # Total number of threads is number of leaders * threads per leader
            total_threads = leaders * num_threads_per_leader
            thread_numbers.append(total_threads)

            print(f"Calculated total time {total_time} seconds for {total_threads} threads in Mode 3.")

    return thread_numbers, total_times

# Plotting functions
def plot_average_leader_times(leader_numbers, average_times):
    plt.figure(figsize=(10, 6))
    plt.plot(leader_numbers, average_times, marker='o', color='skyblue')
    plt.title('Average Leader Computation Time vs. Number of Leaders')
    plt.xlabel('Number of Leaders')
    plt.ylabel('Average Leader Computation Time (seconds)')
    plt.grid(True)
    # ax = plt.gca()  # Get current axis
    # ax.spines['top'].set_visible(False)
    # ax.spines['right'].set_visible(False)
    plt.savefig('average_leader_computation_time.png')
    plt.show()
    print("Saved plot as average_leader_computation_time.png")

def plot_mode_comparison(mode2_threads, mode2_times, mode3_threads, mode3_times):
    plt.figure(figsize=(10, 6))
    plt.plot(mode2_threads, mode2_times, marker='o', color='red', label='Mode 2 (No Load Balancing)')
    plt.plot(mode3_threads, mode3_times, marker='o', color='skyblue', label='Mode 3 (With Load Balancing)')
    plt.title('Total Computation Time vs. Number of Threads')
    plt.xlabel('Number of Threads')
    plt.ylabel('Total Computation Time (seconds)')
    plt.legend()
    plt.grid(True)
    ax = plt.gca()  # Get current axis
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.savefig('mode_comparison.png')
    plt.show()
    print("Saved plot as mode_comparison.png")

def main():
    # Generate the first graph
    leader_numbers, average_times = read_average_leader_times()
    plot_average_leader_times(leader_numbers, average_times)

    # Generate the second graph
    mode2_threads, mode2_times = read_mode2_times()
    mode3_threads, mode3_times = read_mode3_times()
    plot_mode_comparison(mode2_threads, mode2_times, mode3_threads, mode3_times)

if __name__ == "__main__":
    main()
