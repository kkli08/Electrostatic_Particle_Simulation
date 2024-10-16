//
// Created by Damian Li on 2024-10-16.
//

#ifndef SEQUENTIAL_H
#define SEQUENTIAL_H

#include "Particle.h"
#include "ParticleReader.h"
#include <vector>

// Sequential class definition
class Sequential {
private:
    double cutoff_radius;  // Cutoff radius for force calculations
    CoulombForceCalculator forceCalculator;

public:
    // Constructor to initialize with cutoff radius
    Sequential(double cutoff_radius);

    // Method to compute the forces on each particle
    std::vector<double> computeForces(const std::vector<std::unique_ptr<Particle>>& particles);

    // Method to write the computed forces to a CSV file
    void writeForcesToFile(const std::vector<double>& forces, const std::string& outputFilePath);
};



#endif // SEQUENTIAL_H
