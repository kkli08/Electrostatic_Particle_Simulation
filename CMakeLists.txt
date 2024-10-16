# Specify the minimum CMake version required
cmake_minimum_required(VERSION 3.10)

# Set the project name
project(nParticleSim)

# Set the C++ standard to at least C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add the executable target
add_executable(nParticleSim main.cpp)
