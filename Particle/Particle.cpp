//
// Created by Damian Li on 2024-10-16.
//

#include "Particle.h"

// Default constructor for Particle class
Particle::Particle() : x(0.0), y(0.0), charge(0.0) {}

// Constructor for Particle class
Particle::Particle(double xCoord, double yCoord, char type) {
    x = xCoord;
    y = yCoord;
    if (type == 'p') {
        charge = charge_proton;  // Proton
    } else if (type == 'e') {
        charge = charge_electron; // Electron
    }
}

// Get ParticleData for MPI communication
ParticleData Particle::getData() const {
    ParticleData data;
    data.x = x;
    data.y = y;
    data.charge = charge;
    return data;
}

// Set ParticleData from MPI communication
void Particle::setData(const ParticleData& data) {
    x = data.x;
    y = data.y;
    charge = data.charge;
}

// Calculate the distance between two particles
double CoulombForceCalculator::calculateDistance(const Particle& p1, const Particle& p2) const {
    return std::sqrt(std::pow(p2.x - p1.x, 2) + std::pow(p2.y - p1.y, 2)) * 1e-10;  // Distance in meters
}

// Compute the Coulomb force between two particles
double CoulombForceCalculator::computeCoulombForce(const Particle& p1, const Particle& p2) const {
    double r = calculateDistance(p1, p2);
    if (r == 0) return 0;  // Avoid division by zero

    double F = k * (std::abs(p1.charge * p2.charge)) / (r * r);

    // If the charges are opposite, the force is attractive (negative)
    if ((p1.charge > 0 && p2.charge < 0) || (p1.charge < 0 && p2.charge > 0)) {
        F = -F;
    }

    return F;
}
