echo
echo ---------Serial SSSP is starting----------
./sssp
echo ---------Serial SSSP is done--------------
echo
echo ---------OpenMP SSSP is starting----------
./sssp_omp
echo ---------OpenMP SSSP is done--------------
echo
echo ---------MPI SSSP is starting-------------
mpiexec -n 3 ./sssp_mpi
echo ---------MPI SSSP is done-----------------
echo