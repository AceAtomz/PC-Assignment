echo
echo ---------Serial scan is starting----------
./sssp
echo ---------Serial scan is done--------------
echo
echo ---------OpenMP scan is starting----------
./sssp_omp
echo ---------OpenMP scan is done--------------
echo
echo ---------MPI scan is starting-------------
mpiexec -n 4 ./sssp_mpi
echo ---------MPI scan is done-----------------
echo