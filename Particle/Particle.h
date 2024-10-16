//
// Created by Damian Li on 2024-10-16.
//

#ifndef PARTICLE_H
#define PARTICLE_H



#include <cmath>
#include <string>

// Constants
const double k = 8.99e9;                // Coulomb constant in Nm^2/C^2
const double charge_proton = 1.60e-19;  // Charge of a proton in C
const double charge_electron = -1.60e-19; // Charge of an electron in C

// Particle class definition
class Particle {
public:
    double x;       // x-coordinate (in 1e-10 meters)
    double y;       // y-coordinate (in 1e-10 meters)
    double charge;  // Charge in coulombs (1.60e-19 for proton, -1.60e-19 for electron)

    // Constructor to initialize from dataset format (xyq)
    Particle(double xCoord, double yCoord, char type);
};

// CoulombForceCalculator class definition
class CoulombForceCalculator {
public:
    // Function to calculate the distance between two particles
    double calculateDistance(const Particle& p1, const Particle& p2) const;

    // Function to compute Coulomb's force between two particles
    double computeCoulombForce(const Particle& p1, const Particle& p2) const;
};



#endif //PARTICLE_H
