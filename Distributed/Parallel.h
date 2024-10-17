//
// Created by Damian Li on 2024-10-16.
//

#ifndef PARALLEL_H
#define PARALLEL_H

#include "Particle.h"
#include "ParticleReader.h"
#include <vector>
#include <string>
#include <thread>

// Parallel class definition
class Parallel {
private:
    double cutoff_radius;  // Cutoff radius for force calculations
    int num_threads;       // Number of threads to use
    CoulombForceCalculator forceCalculator;

public:
    // Constructor to initialize with cutoff radius and number of threads
    Parallel(double cutoff_radius, int num_threads);

    // Method to compute the forces on each particle in parallel
    std::vector<double> computeForces(const std::vector<std::unique_ptr<Particle>>& particles);

    // Method to write the computed forces to a CSV file
    void writeForcesToFile(const std::vector<double>& forces, const std::string& outputFilePath);

private:
    // Thread worker function to compute forces for a range of particles
    void computeForceRange(const std::vector<std::unique_ptr<Particle>>& particles,
                           std::vector<double>& forces,
                           size_t start_index, size_t end_index);
};

#endif // PARALLEL_H