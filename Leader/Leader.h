//
// Created by Damian Li on 2024-10-19.
//

#ifndef LEADER_H
#define LEADER_H

#include "Particle.h"
#include "ParticleReader.h"
#include <vector>
#include <string>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <deque>
#include <map>
#include <mpi.h>
#include <fstream>

// Leader class definition
class Leader {
private:
    double cutoff_radius;  // Cutoff radius for force calculations
    int num_threads;       // Number of threads per leader
    CoulombForceCalculator forceCalculator;

    // For the dynamic work queue
    std::mutex queue_mutex;
    std::condition_variable queue_cv;
    std::deque<size_t> task_queue;  // Task queue containing indices of particles to process

    // For timing
    std::mutex time_mutex;  // Mutex to protect access to thread_times
    std::map<std::thread::id, double> thread_times;  // Map of thread IDs to computation times

public:
    // Constructor to initialize with cutoff radius and number of threads
    Leader(double cutoff_radius, int num_threads);

    // Method to compute the forces on each particle across leaders (MPI) and threads
    std::vector<double> computeForces(const std::vector<std::unique_ptr<Particle>>& particles);

    // Method to write the computed forces to a CSV file
    void writeForcesToFile(const std::vector<double>& forces, const std::string& outputFilePath);

private:
    // Worker thread function to compute forces for particles from the queue
    void workerThreadFunction(const std::vector<std::unique_ptr<Particle>>& particles,
                              std::vector<double>& forces,
                              std::ofstream& timingFile);
};

#endif // LEADER_H
