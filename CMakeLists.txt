# Specify the minimum CMake version required
cmake_minimum_required(VERSION 3.10)

# Set the project name
project(nParticleSim)

# Set the C++ standard to at least C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add the executable target
add_executable(nParticleSim main.cpp
        Sequential/Sequential.cpp
        Particle/Particle.cpp
        Particle/ParticleReader.cpp)

target_include_directories(nParticleSim PUBLIC
        ${PROJECT_SOURCE_DIR}/Particle
        ${PROJECT_SOURCE_DIR}/Sequential
)