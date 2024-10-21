//
// Created by Damian Li on 2024-10-16.
//

#ifndef PARTICLE_H
#define PARTICLE_H

#include <cmath>
#include <string>

// Constants
const double k = 8.99e9;                  // Coulomb constant in Nm^2/C^2
const double charge_proton = 1.60e-19;    // Charge of a proton in C
const double charge_electron = -1.60e-19; // Charge of an electron in C

// ParticleData struct for MPI communication
struct ParticleData {
    double x;
    double y;
    double charge;
};

// Particle class definition
class Particle {
public:
    double x;       // x-coordinate (in 1e-10 meters)
    double y;       // y-coordinate (in 1e-10 meters)
    double charge;  // Charge in coulombs (1.60e-19 for proton, -1.60e-19 for electron)

    // Constructors
    Particle();  // Default constructor
    Particle(double xCoord, double yCoord, char type);

    // Methods for MPI communication
    ParticleData getData() const;
    void setData(const ParticleData& data);
};

// CoulombForceCalculator class definition
class CoulombForceCalculator {
public:
    // Function to calculate the distance between two particles
    double calculateDistance(const Particle& p1, const Particle& p2) const;

    // Function to compute Coulomb's force between two particles
    double computeCoulombForce(const Particle& p1, const Particle& p2) const;
};

#endif // PARTICLE_H
