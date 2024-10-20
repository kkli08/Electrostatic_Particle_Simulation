# Specify the minimum CMake version required
cmake_minimum_required(VERSION 3.10)

# Set the project name
project(nParticleSim)

# Set the C++ standard to at least C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Find and link MPI
find_package(MPI REQUIRED)

# Add the executable target
add_executable(nParticleSim main.cpp
        Sequential/Sequential.cpp
        Particle/Particle.cpp
        Particle/ParticleReader.cpp
        Distributed/Parallel.cpp
        Leader/Leader.cpp)

# Include directories for the project
target_include_directories(nParticleSim PUBLIC
        ${PROJECT_SOURCE_DIR}/Particle
        ${PROJECT_SOURCE_DIR}/Sequential
        ${PROJECT_SOURCE_DIR}/Distributed
        ${PROJECT_SOURCE_DIR}/Leader
        ${MPI_INCLUDE_PATH}  # Include MPI headers
)

# Link MPI libraries to the target
target_link_libraries(nParticleSim PUBLIC MPI::MPI_CXX)

# Enable threading support (optional, if you use pthreads or std::thread)
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
target_link_libraries(nParticleSim PUBLIC Threads::Threads)
