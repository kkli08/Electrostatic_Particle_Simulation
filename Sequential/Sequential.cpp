//
// Created by Damian Li on 2024-10-16.
//

#include "Sequential.h"
#include <iostream>
#include <fstream>

// Constructor to initialize with cutoff radius
Sequential::Sequential(double cutoff_radius) : cutoff_radius(cutoff_radius) {}

// Method to compute the forces on each particle
std::vector<double> Sequential::computeForces(const std::vector<std::unique_ptr<Particle>>& particles) {
    std::vector<double> forces(particles.size(), 0.0);  // Store the force sum for each particle

    // Iterate over each particle
    for (size_t i = 0; i < particles.size(); ++i) {
        const Particle& p1 = *particles[i];

        // Sum the forces from other particles within the cutoff radius
        for (size_t j = 0; j < particles.size(); ++j) {
            if (i == j) continue;  // Skip the same particle

            const Particle& p2 = *particles[j];
            double distance = forceCalculator.calculateDistance(p1, p2);

            // Only consider the particles within the cutoff radius
            if (distance <= cutoff_radius) {
                double force = forceCalculator.computeCoulombForce(p1, p2);
                forces[i] += force;  // Sum up the scalar force (signed)
            }
        }
    }

    return forces;  // Return the summed forces for all particles
}

// Method to write the computed forces to a CSV file
void Sequential::writeForcesToFile(const std::vector<double>& forces, const std::string& outputFilePath) {
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