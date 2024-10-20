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
#include <mpi.h>
#include <fstream>

// Leader class definition
class Leader {
private:
    double cutoff_radius;  // Cutoff radius for force calculations
    int num_threads;       // Number of threads per leader
    int num_leaders;       // Number of leader processes (using MPI)
    CoulombForceCalculator forceCalculator;

public:
    // Constructor to initialize with cutoff radius, number of threads, and number of leaders
    Leader(double cutoff_radius, int num_threads, int num_leaders);

    // Method to compute the forces on each particle across leaders (MPI) and threads
    std::vector<double> computeForces(const std::vector<std::unique_ptr<Particle>>& particles);

    // Method to write the computed forces to a CSV file
    void writeForcesToFile(const std::vector<double>& forces, const std::string& outputFilePath);

private:
    // Thread worker function to compute forces for a range of particles
    void computeForceRange(const std::vector<std::unique_ptr<Particle>>& particles,
                           std::vector<double>& forces,
                           size_t start_index, size_t end_index,
                           std::ofstream& timingFile);

    // Method to handle MPI parallelization among leader processes
    void distributeParticlesMPI(const std::vector<std::unique_ptr<Particle>>& particles,
                                std::vector<double>& forces);
};

#endif // LEADER_H
