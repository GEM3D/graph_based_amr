
Parallel Adaotive Mesh Refinement using Unstructured Grid Connectivity Approach for Cartesian Domains

  - Cartesian Processor Topology for Generation of the Initial Mesh
  - Usesd Zoltan (Sandia National Lab) for Grid Partitioning 
  - The treeis implicit as all the modifications are done on the connectivity

# Required Libraries

  - HDF5
  - Cmake
  - Zoltan
  - Paraview

# Input Parameters specified in text file run.txt
 - Cartesian Process Topology, i.e. number of processes in x-y-z directions
 - Maximum level of refinement

# Run 
  - cd build
  - cmake ..
  - cd .. 
  - mpirun -np xxx  amrGem input/geometry.stl <run.txt
  - The result will be written in  /soln folder
