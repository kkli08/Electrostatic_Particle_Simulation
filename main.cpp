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
#include <chrono>

using namespace std;
namespace fs = std::filesystem;

// Function to parse command line arguments
void parseArguments(int argc, char* argv[], int& mode, double& cutoff_radius, std::string& input_file, int& num_threads) {
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
    std::string input_file;

    try {
        parseArguments(argc, argv, mode, cutoff_radius, input_file, num_threads);

        // Display the parsed values
        std::cout << "Mode: " << mode << std::endl;
        std::cout << "Cutoff Radius: " << cutoff_radius << " meters" << std::endl;
        std::cout << "Input File: " << input_file << std::endl;

        if (mode == 1) {
            // Sequential Computation
            ParticleReader reader;
            auto particles = reader.readParticles(input_file);
            Sequential sequentialProcessor(cutoff_radius);
            auto forces = sequentialProcessor.computeForces(particles);
            std::string outputFilePath = "output.csv";
            sequentialProcessor.writeForcesToFile(forces, outputFilePath);
        } else if (mode == 2) {
            // Evenly-Distributed Parallel Computation
            ParticleReader reader;
            auto particles = reader.readParticles(input_file);
            Parallel parallelProcessor(cutoff_radius, num_threads);
            auto forces = parallelProcessor.computeForces(particles);
            std::string outputFilePath = "output_parallel.csv";
            parallelProcessor.writeForcesToFile(forces, outputFilePath);
        } else if (mode == 3) {
            // Initialize MPI
            MPI_Init(&argc, &argv);

            int world_rank = 0;
            int world_size = 1;
            MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
            MPI_Comm_size(MPI_COMM_WORLD, &world_size);

            if (world_size == 1) {
                if (world_rank == 0) {
                    std::cerr << "Error: Mode 3 requires running with multiple processes using mpirun or mpiexec." << std::endl;
                    std::cerr << "Example: mpirun -np <num_leaders> ./nParticleSim --mode=3 ..." << std::endl;
                }
                MPI_Finalize();
                return 1;
            }

            ParticleReader reader;
            std::vector<std::unique_ptr<Particle>> particles;
            size_t num_particles = 0;

            // Timing variables
            std::chrono::high_resolution_clock::time_point read_start, read_end;
            std::chrono::duration<double> read_duration(0);
            std::chrono::high_resolution_clock::time_point bcast_start, bcast_end;
            std::chrono::duration<double> bcast_duration(0);

            if (world_rank == 0) {
                // Start timing for reading particles
                read_start = std::chrono::high_resolution_clock::now();

                // Only the root process reads the particles from the file
                particles = reader.readParticles(input_file);
                num_particles = particles.size();

                // End timing for reading particles
                read_end = std::chrono::high_resolution_clock::now();
                read_duration = read_end - read_start;
                std::cout << "Time to read particles: " << read_duration.count() << " seconds." << std::endl;
            }

            // Start timing for broadcasting particles
            bcast_start = std::chrono::high_resolution_clock::now();

            // Broadcast the number of particles to all processes
            MPI_Bcast(&num_particles, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);

            if (world_rank != 0) {
                // Other processes allocate the vector of particles
                particles.resize(num_particles);
                for (size_t i = 0; i < num_particles; ++i) {
                    particles[i] = std::make_unique<Particle>();
                }
            }

            // Define MPI_Datatype for ParticleData
            MPI_Datatype MPI_ParticleData;
            int blocklengths[3] = {1, 1, 1};
            MPI_Aint offsets[3];
            offsets[0] = offsetof(ParticleData, x);
            offsets[1] = offsetof(ParticleData, y);
            offsets[2] = offsetof(ParticleData, charge);
            MPI_Datatype types[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
            MPI_Type_create_struct(3, blocklengths, offsets, types, &MPI_ParticleData);
            MPI_Type_commit(&MPI_ParticleData);

            // Prepare an array of ParticleData
            std::vector<ParticleData> particleDataArray(num_particles);

            if (world_rank == 0) {
                // Root process fills particleDataArray
                for (size_t i = 0; i < num_particles; ++i) {
                    particleDataArray[i] = particles[i]->getData();
                }
            }

            // Broadcast the array of ParticleData to all processes
            MPI_Bcast(particleDataArray.data(), num_particles, MPI_ParticleData, 0, MPI_COMM_WORLD);

            if (world_rank != 0) {
                // Other processes set their particle data
                for (size_t i = 0; i < num_particles; ++i) {
                    particles[i]->setData(particleDataArray[i]);
                }
            }

            // Free the MPI datatype
            MPI_Type_free(&MPI_ParticleData);

            // End timing for broadcasting particles
            bcast_end = std::chrono::high_resolution_clock::now();
            bcast_duration = bcast_end - bcast_start;

            // Output timing information
            if (world_rank == 0) {
                std::cout << "Time to broadcast particles: " << bcast_duration.count() << " seconds." << std::endl;
            }

            // Create leader processor
            Leader leaderProcessor(cutoff_radius, num_threads);

            // Compute
            auto forces = leaderProcessor.computeForces(particles);

            // Output
            std::string outputFilePath = "output_leader.csv";
            leaderProcessor.writeForcesToFile(forces, outputFilePath);

            // Finalize MPI
            MPI_Finalize();
        } else {
            std::cerr << "Error: Invalid mode selected. Mode should be 1, 2, or 3." << std::endl;
            return 1;
        }
    } catch (std::exception& e) {
        std::cerr << "Failed to parse arguments: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}