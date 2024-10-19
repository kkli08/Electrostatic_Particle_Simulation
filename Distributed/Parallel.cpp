//
// Created by Damian Li on 2024-10-16.
//

#include "Parallel.h"
#include <fstream>
#include <iostream>
#include <cmath>
#include <chrono>  // Include chrono for timing

// Constructor to initialize with cutoff radius and number of threads
Parallel::Parallel(double cutoff_radius, int num_threads)
    : cutoff_radius(cutoff_radius), num_threads(num_threads) {}

// Method to compute the forces on each particle in parallel
std::vector<double> Parallel::computeForces(const std::vector<std::unique_ptr<Particle>>& particles) {
    std::vector<double> forces(particles.size(), 0.0);  // Store the force sum for each particle
    std::vector<std::thread> threads;
    size_t num_particles = particles.size();
    size_t chunk_size = (num_particles + num_threads - 1) / num_threads;  // Calculate chunk size per thread

    // Dynamically generate the output filename based on cutoff radius and num_threads
    std::string outputFileName = "timing_cutoff_" + std::to_string(cutoff_radius) + "_threads_" + std::to_string(num_threads) + ".txt";
    std::ofstream timingFile(outputFileName);

    if (!timingFile.is_open()) {
        std::cerr << "Error: Could not open timing file " << outputFileName << std::endl;
        return forces;
    }

    // Record start time for total computation
    auto total_start = std::chrono::high_resolution_clock::now();

    // Create and launch threads
    for (int t = 0; t < num_threads; ++t) {
        size_t start_index = t * chunk_size;
        size_t end_index = std::min(start_index + chunk_size, num_particles);

        if (start_index < num_particles) {
            threads.emplace_back(&Parallel::computeForceRange, this, std::cref(particles), std::ref(forces), start_index, end_index, std::ref(timingFile));
        }
    }

    // Join all threads
    for (auto& thread : threads) {
        thread.join();
    }

    // Record end time for total computation
    auto total_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> total_duration = total_end - total_start;

    // Print and write the total time taken
    std::cout << "Total computation time: " << total_duration.count() << " seconds." << std::endl;
    timingFile << "Total computation time: " << total_duration.count() << " seconds." << std::endl;

    // Close the timing file
    timingFile.close();

    return forces;  // Return the summed forces for all particles
}

// Worker function to compute forces for a range of particles
void Parallel::computeForceRange(const std::vector<std::unique_ptr<Particle>>& particles,
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
void Parallel::writeForcesToFile(const std::vector<double>& forces, const std::string& outputFilePath) {
    std::ofstream outputFile(outputFilePath);

    if (!outputFile.is_open()) {
        std::cerr << "Error: Could not open output file " << outputFilePath << std::endl;
        return;
    }

    // Write the force data, one value per line
    for (const auto& force : forces) {
        outputFile << force << "\n";  // Write force value, no comma or index
    }

    outputFile.close();
    std::cout << "Forces written to " << outputFilePath << std::endl;
}
