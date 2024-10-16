//
// Created by Damian Li on 2024-10-16.
//

#ifndef PARTICLEREADER_H
#define PARTICLEREADER_H

#include "Particle.h"
#include <vector>
#include <string>

// ParticleReader class definition
class ParticleReader {
public:
    // Method to read particles from a CSV file
    std::vector<std::unique_ptr<Particle>> readParticles(const std::string& filePath);
};


#endif //PARTICLEREADER_H
