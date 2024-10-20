//
// Created by Damian Li on 2024-10-16.
//

#include <iostream>
#include <string>
#include <vector>
#include <filesystem>
#include <stdexcept>
#include "Sequential.h"
#include "Parallel.h"
#include "ParticleReader.h"
#include "Leader.h"


using namespace std;
namespace fs = std::filesystem;

// Function to parse command line arguments
void parseArguments(int argc, char* argv[], int& mode, double& cutoff_radius, std::string& input_file, int& num_threads, int& num_leaders) {
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];

    // Check if the argument starts with "--mode="
    if (arg.find("--mode=") == 0) {
      std::string mode_value = arg.substr(7);  // Extract the value after "--mode="
      try {
        mode = std::stoi(mode_value);  // Convert to int
        if (mode < 1 || mode > 3) {
          throw std::invalid_argument("Invalid mode value.");
        }
      } catch (std::exception& e) {
        std::cerr << "Error: Invalid mode value. Expected 1, 2, or 3." << std::endl;
        throw;
      }
    }
    // Check if the argument starts with "--cutoff_radius="
    else if (arg.find("--cutoff_radius=") == 0) {
      std::string radius_value = arg.substr(16);  // Extract the value after "--cutoff_radius="
      try {
        cutoff_radius = std::stod(radius_value);  // Convert to double
        if (cutoff_radius <= 0) {
          throw std::invalid_argument("Cutoff radius must be positive.");
        }
        cutoff_radius *= 1e-10;
      } catch (std::exception& e) {
        std::cerr << "Error: Invalid cutoff radius value." << std::endl;
        throw;
      }
    }
    else if (arg.find("--num_threads=") == 0) {
      std::string num_threads_value = arg.substr(14);
      try {
        num_threads = std::stoi(num_threads_value);  // Convert to int
        if (num_threads < 1) {
          throw std::invalid_argument("Invalid num_threads value.");
        }
      } catch (std::exception& e) {
        std::cerr << "Error: Invalid num_threads value." << std::endl;
        throw;
      }
    }
    else if (arg.find("--num_leaders=") == 0) {
      std::string num_leaders_value = arg.substr(14);
      try {
        num_leaders = std::stoi(num_leaders_value);  // Convert to int
        if (num_leaders < 1) {
          throw std::invalid_argument("Invalid num_threads value.");
        }
      } catch (std::exception& e) {
        std::cerr << "Error: Invalid num_threads value." << std::endl;
        throw;
      }
    }
    // Check if the argument starts with "--input="
    else if (arg.find("--input=") == 0) {
      input_file = arg.substr(8);  // Extract the value after "--input="

      // Check if the file exists using std::filesystem
      if (!fs::exists(input_file)) {
        std::cerr << "Error: Input file does not exist: " << input_file << std::endl;
        throw std::invalid_argument("Invalid input file.");
      }

      // Check if it's a regular file, not a directory
      if (!fs::is_regular_file(input_file)) {
        std::cerr << "Error: Input path is not a valid file: " << input_file << std::endl;
        throw std::invalid_argument("Invalid input file.");
      }

    }
    else {
      std::cerr << "Error: Unrecognized argument: " << arg << std::endl;
      throw std::invalid_argument("Invalid argument.");
    }
  }
}

int main(int argc, char* argv[]) {
  int mode = 1;
  double cutoff_radius = 0.0;
  int num_threads = 1;
  int num_leaders = 1;
  std::string input_file;

  try {
    parseArguments(argc, argv, mode, cutoff_radius, input_file, num_threads, num_leaders);

    // Display the parsed values
    std::cout << "Mode: " << mode << std::endl;
    std::cout << "Cutoff Radius: " << cutoff_radius << " meters" << std::endl;
    std::cout << "Input File: " << input_file << std::endl;

    ParticleReader reader;

    if (mode == 1) {
      // Sequential Computation
      // Read particles from the file
      auto particles = reader.readParticles(input_file);

      Sequential sequentialProcessor(cutoff_radius);

      // Compute the forces
      auto forces = sequentialProcessor.computeForces(particles);

      // Output the forces to a CSV file
      std::string outputFilePath = "output.csv";
      sequentialProcessor.writeForcesToFile(forces, outputFilePath);

    } else if (mode == 2) {
      // Evenly-Distributed Parallel Computation
      // Read particles from the file
      auto particles = reader.readParticles(input_file);

      Parallel parallelProcessor(cutoff_radius, num_threads);

      // Compute the forces in parallel
      auto forces = parallelProcessor.computeForces(particles);

      // Output the forces to a CSV file
      std::string outputFilePath = "output_parallel.csv";
      parallelProcessor.writeForcesToFile(forces, outputFilePath);


    } else if (mode == 3) {
      // Load-Balanced, Leader-Based Parallel Computation
      // Initialize MPI
      MPI_Init(&argc, &argv);
      // Read particles from the file
      auto particles = reader.readParticles(input_file);

      Leader leaderProcessor(cutoff_radius, num_threads, num_leaders);

      // Compute
      auto forces = leaderProcessor.computeForces(particles);

      // Output
      std::string outputFilePath = "output_leader.csv";
      leaderProcessor.writeForcesToFile(forces, outputFilePath);

      // Finalize MPI
      MPI_Finalize();
    }
  }
  catch (std::exception& e) {
    std::cerr << "Failed to parse arguments: " << e.what() << std::endl;
    return 1;
  }

  return 0;
}