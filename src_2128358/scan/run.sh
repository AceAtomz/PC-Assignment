echo
echo ---------Serial scan is starting----------
./scan
echo ---------Serial scan is done--------------
echo
echo ---------OpenMP scan is starting----------
./scan_omp
echo ---------OpenMP scan is done--------------
echo
echo ---------MPI scan is starting-------------
mpiexec -n 4 ./scan_mpi
echo ---------MPI scan is done-----------------
echo