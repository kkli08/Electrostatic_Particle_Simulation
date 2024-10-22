import subprocess
import os
import shutil

# List of np values to replace {%d}
np_values = [3, 4, 5, 6, 7]

# Get the current working directory
current_dir = os.getcwd()

for np_value in np_values:
    # Build the command with the current np_value
    cmd = f"mpirun -np {np_value} ./nParticleSim --mode=3 --cutoff_radius=45000 --input=../dataset/particles.csv --num_threads=10"
    print(f"Running command: {cmd}")

    # Record the list of files before execution
    before_files = set(os.listdir(current_dir))

    # Execute the command and wait for it to finish
    process = subprocess.run(cmd, shell=True)

    # Check if the command executed successfully
    if process.returncode != 0:
        print(f"Command failed with return code {process.returncode}. Exiting.")
        break  # Stop the script if the command fails

    # Record the list of files after execution
    after_files = set(os.listdir(current_dir))

    # Identify newly created files
    new_files = after_files - before_files

    # Create a directory for the current np_value if it doesn't exist
    output_dir = os.path.join(current_dir, str(np_value))
    os.makedirs(output_dir, exist_ok=True)

    # Move the new output files to the corresponding directory
    for file_name in new_files:
        src_path = os.path.join(current_dir, file_name)
        dst_path = os.path.join(output_dir, file_name)
        shutil.move(src_path, dst_path)
        print(f"Moved {file_name} to {output_dir}")

    print(f"Completed run with np={np_value}. Outputs moved to {output_dir}\n")

print("All commands executed successfully.")
