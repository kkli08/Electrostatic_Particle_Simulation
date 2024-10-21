//
// Created by Damian Li on 2024-10-19.
//

#include "Leader.h"
#include <iostream>
#include <mpi.h>
#include <cmath>
#include <fstream>
#include <chrono>
#include <thread>

// Constructor to initialize with cutoff radius and number of threads
Leader::Leader(double cutoff_radius, int num_threads)
    : cutoff_radius(cutoff_radius), num_threads(num_threads) {}

// Method to compute the forces on each particle using MPI and threads, with time tracking
std::vector<double> Leader::computeForces(const std::vector<std::unique_ptr<Particle>>& particles) {
    int world_rank, world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);  // Get the number of processes (leaders)
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);  // Get the rank of the process (leader ID)

    size_t num_particles = particles.size();
    size_t particles_per_leader = num_particles / world_size;
    size_t start_index = world_rank * particles_per_leader;
    size_t end_index = (world_rank == world_size - 1) ? num_particles : start_index + particles_per_leader;

    std::vector<double> forces(num_particles, 0.0);  // Forces for all particles

    // Initialize the work queue with particle indices for this leader's partition
    {
        std::lock_guard<std::mutex> lock(queue_mutex);
        for (size_t i = start_index; i < end_index; ++i) {
            task_queue.push_back(i);
        }
    }

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

    // Create and launch worker threads
    std::vector<std::thread> threads;
    for (int t = 0; t < num_threads; ++t) {
        threads.emplace_back(&Leader::workerThreadFunction, this, std::cref(particles), std::ref(forces), std::ref(timingFile));
    }

    // Join all threads
    for (auto& thread : threads) {
        thread.join();
    }

    // Record end time for total computation
    auto total_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> total_duration = total_end - total_start;

    // Compute average time per thread
    double total_thread_time = 0.0;
    for (const auto& entry : thread_times) {
        total_thread_time += entry.second;
    }
    double avg_thread_time = total_thread_time / num_threads;

    // Compute average time per particle
    size_t num_particles_processed = end_index - start_index;
    double avg_particle_time = total_duration.count() / num_particles_processed;

    // Print and write the total time taken by the leader
    std::cout << "Leader " << world_rank << " total computation time: " << total_duration.count() << " seconds." << std::endl;
    timingFile << "Leader " << world_rank << " total computation time: " << total_duration.count() << " seconds." << std::endl;

    // Print and write average time per thread
    std::cout << "Leader " << world_rank << " average thread computation time: " << avg_thread_time << " seconds." << std::endl;
    timingFile << "Leader " << world_rank << " average thread computation time: " << avg_thread_time << " seconds." << std::endl;

    // Print and write average time per particle
    std::cout << "Leader " << world_rank << " average time per particle: " << avg_particle_time << " seconds." << std::endl;
    timingFile << "Leader " << world_rank << " average time per particle: " << avg_particle_time << " seconds." << std::endl;

    // Combine forces from all leaders using MPI_Reduce (root process gathers results)
    std::vector<double> global_forces(num_particles, 0.0);
    MPI_Reduce(forces.data(), global_forces.data(), num_particles, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    // Close the timing file
    timingFile.close();

    // Collect total computation times from all leaders
    double leader_total_time = total_duration.count();
    double total_leader_time = 0.0;
    MPI_Reduce(&leader_total_time, &total_leader_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    // Root process calculates average leader time
    if (world_rank == 0) {
        double avg_leader_time = total_leader_time / world_size;
        std::cout << "Average leader computation time: " << avg_leader_time << " seconds." << std::endl;

        // Write to a new summary file
        std::ofstream summaryFile("timing_summary.txt", std::ios::app);
        if (summaryFile.is_open()) {
            summaryFile << "Average leader computation time: " << avg_leader_time << " seconds." << std::endl;
            summaryFile.close();
        } else {
            std::cerr << "Error: Could not open timing summary file." << std::endl;
        }
    }

    // Only return the global forces from the root process (rank 0)
    if (world_rank == 0) {
        return global_forces;
    } else {
        return {};
    }
}

// Worker thread function to compute forces for particles from the queue
void Leader::workerThreadFunction(const std::vector<std::unique_ptr<Particle>>& particles,
                                  std::vector<double>& forces,
                                  std::ofstream& timingFile) {
    auto thread_start = std::chrono::high_resolution_clock::now();
    size_t num_particles_processed = 0;

    while (true) {
        size_t particle_index;

        // Fetch a task from the queue
        {
            std::lock_guard<std::mutex> lock(queue_mutex);
            if (task_queue.empty()) {
                break;  // Exit if no more tasks
            }
            particle_index = task_queue.front();
            task_queue.pop_front();
        }

        num_particles_processed++;

        // Record start time for this task
        auto task_start = std::chrono::high_resolution_clock::now();

        const Particle& p1 = *particles[particle_index];

        // Compute the forces from all other particles
        for (size_t j = 0; j < particles.size(); ++j) {
            if (particle_index == j) continue;  // Skip the same particle

            const Particle& p2 = *particles[j];
            double distance = forceCalculator.calculateDistance(p1, p2);

            // Only consider particles within the cutoff radius
            if (distance <= cutoff_radius) {
                double force = forceCalculator.computeCoulombForce(p1, p2);
                // Since only this thread updates forces[particle_index], no need for synchronization
                forces[particle_index] += force;  // Sum up the scalar force (signed)
            }
        }

        // Record end time for this task
        auto task_end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> task_duration = task_end - task_start;

        // Optional: Log the time taken for this task
        timingFile << "Thread " << std::this_thread::get_id()
                   << " processed particle " << particle_index
                   << " in " << task_duration.count() << " seconds." << std::endl;
    }

    auto thread_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> thread_duration = thread_end - thread_start;

    // Store the total computation time for this thread
    {
        std::lock_guard<std::mutex> lock(time_mutex);
        thread_times[std::this_thread::get_id()] = thread_duration.count();
    }

    // Optional: Log the total time taken by this thread
    timingFile << "Thread " << std::this_thread::get_id()
               << " total computation time: " << thread_duration.count() << " seconds."
               << " Processed " << num_particles_processed << " particles." << std::endl;
}

// Method to write the computed forces to a CSV file
void Leader::writeForcesToFile(const std::vector<double>& forces, const std::string& outputFilePath) {
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Only the root process writes the output
    if (world_rank == 0) {
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
}
