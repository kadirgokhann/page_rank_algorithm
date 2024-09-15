# PageRank Algorithm Implementations

## 1. Sequential PageRank Implementation (C++)
This version of the algorithm processes the PageRank calculations in a single-threaded manner. The algorithm iteratively computes the ranking of web pages based on their incoming links. The computation continues until the difference between PageRank values of consecutive iterations falls below a set threshold. While this method is straightforward, its performance is limited by the sequential nature of the algorithm, making it inefficient for processing large-scale graphs due to the high computation time.

## 2. Parallel PageRank with MPI (Message Passing Interface)
In this version, the PageRank algorithm is parallelized using MPI, which is suitable for distributed-memory systems. The input graph is divided across multiple processors, each handling a portion of the computation. The PageRank values are computed locally, and at the end of each iteration, the values are exchanged among processors to ensure global consistency. This parallel approach reduces the computation time significantly, but its efficiency depends on the communication overhead between processors.

## 3. Parallel PageRank with OpenMPI
OpenMPI is a specific implementation of MPI that provides an efficient platform for parallel processing in both shared-memory and distributed-memory environments. The OpenMPI-based PageRank algorithm benefits from optimizations in communication and process management provided by the OpenMPI library. Similar to MPI, this version divides the graph and performs PageRank computations across multiple nodes or cores, improving scalability for large datasets while leveraging OpenMPI’s cross-platform compatibility and performance features.

## 4. Parallel PageRank with Thrust
In this implementation, the Thrust library (part of CUDA) is used to accelerate PageRank computations on GPUs. Thrust provides a high-level interface for parallel algorithms and is optimized for large-scale data processing. The PageRank algorithm is adapted to run on the GPU, taking advantage of the massive parallelism offered by modern GPUs to perform computations on the entire graph simultaneously. This approach significantly boosts performance for large-scale graphs, leveraging the GPU's ability to perform thousands of calculations in parallel.
