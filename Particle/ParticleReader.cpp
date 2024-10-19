//
// Created by Damian Li on 2024-10-16.
//

#include "ParticleReader.h"
#include <fstream>
#include <sstream>
#include <iostream>

// Helper function to parse a CSV line and create a Particle
std::unique_ptr<Particle> parseParticle(const std::string& line) {
    std::stringstream ss(line);
    std::string xStr, yStr, qStr;

    // Parse the line
    if (std::getline(ss, xStr, ',') && std::getline(ss, yStr, ',') && std::getline(ss, qStr)) {
        double x = std::stod(xStr);
        double y = std::stod(yStr);
        char chargeType = qStr[0];

        // Create a Particle object
        return std::make_unique<Particle>(x, y, chargeType);
    }

    return nullptr; // Return null if the line is invalid
}

// Method to read particles from a CSV file
std::vector<std::unique_ptr<Particle>> ParticleReader::readParticles(const std::string& filePath) {
    std::vector<std::unique_ptr<Particle>> particles;
    std::ifstream file(filePath);

    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filePath << std::endl;
        return particles;
    }

    std::string line;
    while (std::getline(file, line)) {
        auto particle = parseParticle(line);
        if (particle) {
            particles.push_back(std::move(particle));
        }
    }

    file.close();
    return particles;
}

