/* File:     sssp_mpi.c
 * Purpose:  Serial implementation of Dijkstra's SSSP operation
 *
 * Input:    graph.txt file holding the arrangement of the graph
 * Output:   If the sssp was validated or not.
 *           sssp.out holds average times for plotting
 *
 * Compile:  gcc -g -Wall sssp_mpi.c -o sssp_mpi -lm
 * Run:      mpiexec -n 4 ./sssp_mpi
 *
 * Algorithm:
 *    1.  For a weighted graph G = (V,E,w), where: V is the set of vertices,
 *                                                 E is the set of edges,
 *                                                 w contains the weights
 *    2.  Dijkstra’s SSSP algorithm finds the shorted paths from a source vertex 
 *        s ∈ V to all other vertices in V
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

int main(int argc, char *argv[]) {
    printf("Hello World!");
    return 0;
}