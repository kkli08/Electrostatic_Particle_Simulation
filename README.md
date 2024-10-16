# Electrostatic_Particle_Simulation
A Range-Limited Electrostatic N-Particle Simulation.

## Mode 1: Sequential Computation
> This implementation should be entirely serial. No multithreading, just the approximation method you 
applied to compute the signed scalar force sums on every particle. The cutoff radius must be an input 
parameter for this mode.


## Mode 2: Evenly-Distributed Parallel Computation
> In this implementation, you will use the Pthread/thread execution model to create multiple threads and 
divide the computation among the threads. You must divide the dataset among the threads at the point in 
time when the thread is created, so that once it finishes its given portion of work, it returns. Divide 
the work as evenly as possible among the threads. The number of threads and the cutoff radius must be input 
parameters for this mode.


## Mode 3: Load-Balanced, Leader-Based Parallel Computation
> In this implementation, you will begin by creating leader processes using MPI. Each leader must be given an 
equal partition of the dataset. Each leader creates a pool of worker threads in the form of Pthreads/threads. 
Each leader’s partition of the data must be further partitioned into smaller chunks and placed into a queue that 
can be accessed by all of its worker threads. Worker threads must take one small chunk of data at a time, execute 
the necessary computation, and then return to the queue to take more work. Threads only return once the queue is 
empty and all work is done. The number of leader processes and worker threads and the cutoff radius must be input 
parameters for this mode.