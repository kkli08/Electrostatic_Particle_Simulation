//
// Created by Damian Li on 2024-10-19.
//

#include "Leader.h"
#include <iostream>
#include <mpi.h>
#include <cmath>
#include <fstream>
#include <chrono>  // Include chrono for timing
#include <thread>

// Constructor to initialize with cutoff radius, number of threads, and number of leaders
Leader::Leader(double cutoff_radius, int num_threads, int num_leaders)
    : cutoff_radius(cutoff_radius), num_threads(num_threads), num_leaders(num_leaders) {}

// Method to compute the forces on each particle using MPI and threads, with time tracking
std::vector<double> Leader::computeForces(const std::vector<std::unique_ptr<Particle>>& particles) {
    // Initialize MPI
    int world_rank, world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);  // Get the number of processes (leaders)
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);  // Get the rank of the process (leader ID)

    size_t num_particles = particles.size();
    size_t particles_per_leader = num_particles / world_size;
    size_t start_index = world_rank * particles_per_leader;
    size_t end_index = (world_rank == world_size - 1) ? num_particles : start_index + particles_per_leader;

    std::vector<double> forces(num_particles, 0.0);  // Forces for all particles

    // Dynamically generate the output filename based on cutoff radius, num_threads, and leader rank
    std::string outputFileName = "timing_cutoff_" + std::to_string(cutoff_radius) +
                                 "_threads_" + std::to_string(num_threads) +
                                 "_leader_" + std::to_string(world_rank) + ".txt";
    std::ofstream timingFile(outputFileName);

    if (!timingFile.is_open()) {
        std::cerr << "Error: Could not open timing file " << outputFileName << std::endl;
        return forces;
    }

    // Record start time for total computation
    auto total_start = std::chrono::high_resolution_clock::now();

    // Distribute the work among the threads in this leader
    std::vector<std::thread> threads;
    size_t chunk_size = (particles_per_leader + num_threads - 1) / num_threads;

    // Create and launch threads
    for (int t = 0; t < num_threads; ++t) {
        size_t thread_start = start_index + t * chunk_size;
        size_t thread_end = std::min(thread_start + chunk_size, end_index);

        if (thread_start < end_index) {
            threads.emplace_back(&Leader::computeForceRange, this, std::cref(particles), std::ref(forces), thread_start, thread_end, std::ref(timingFile));
        }
    }

    // Join all threads
    for (auto& thread : threads) {
        thread.join();
    }

    // Record end time for total computation
    auto total_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> total_duration = total_end - total_start;

    // Print and write the total time taken by the leader
    std::cout << "Leader " << world_rank << " total computation time: " << total_duration.count() << " seconds." << std::endl;
    timingFile << "Leader " << world_rank << " total computation time: " << total_duration.count() << " seconds." << std::endl;

    // Combine forces from all leaders using MPI_Reduce (root process gathers results)
    std::vector<double> global_forces(num_particles, 0.0);
    MPI_Reduce(forces.data(), global_forces.data(), num_particles, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);


    // Close the timing file
    timingFile.close();

    // Only return the global forces from the root process (rank 0)
    if (world_rank == 0) {
        return global_forces;
    } else {
        return {};
    }
}

// Worker function to compute forces for a range of particles using threads, with time tracking
void Leader::computeForceRange(const std::vector<std::unique_ptr<Particle>>& particles,
                               std::vector<double>& forces,
                               size_t start_index, size_t end_index,
                               std::ofstream& timingFile) {
    // Record start time for each thread
    auto thread_start = std::chrono::high_resolution_clock::now();

    for (size_t i = start_index; i < end_index; ++i) {
        const Particle& p1 = *particles[i];

        // Compute the forces from all other particles
        for (size_t j = 0; j < particles.size(); ++j) {
            if (i == j) continue;  // Skip the same particle

            const Particle& p2 = *particles[j];
            double distance = forceCalculator.calculateDistance(p1, p2);

            // Only consider particles within the cutoff radius
            if (distance <= cutoff_radius) {
                double force = forceCalculator.computeCoulombForce(p1, p2);
                forces[i] += force;  // Sum up the scalar force (signed)
            }
        }
    }

    // Record end time for each thread
    auto thread_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> thread_duration = thread_end - thread_start;

    // Print and write the time taken by this thread
    std::cout << "Thread handling particles from " << start_index << " to " << end_index
              << " took " << thread_duration.count() << " seconds." << std::endl;
    timingFile << "Thread handling particles from " << start_index << " to " << end_index
               << " took " << thread_duration.count() << " seconds." << std::endl;
}

// Method to write the computed forces to a CSV file
void Leader::writeForcesToFile(const std::vector<double>& forces, const std::string& outputFilePath) {
    std::ofstream outputFile(outputFilePath);

    if (!outputFile.is_open()) {
        std::cerr << "Error: Could not open output file " << outputFilePath << std::endl;
        return;
    }

    // Write the force data, one value per line
    for (const auto& force : forces) {
        outputFile << force << "\n";  // Write force value
    }

    outputFile.close();
    std::cout << "Forces written to " << outputFilePath << std::endl;
}


