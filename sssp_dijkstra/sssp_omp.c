/* File:     sssp_omp.c
 * Purpose:  Serial implementation of Dijkstra's SSSP operation
 *
 * Input:    graph.txt file holding the arrangement of the graph
 * Output:   If sssp_omp was validated or not.
 *           sssp.out holds average times for plotting
 *
 * Compile:  gcc -g -Wall -fopen sssp_omp.c -o sssp_omp
 * Run:      ./sssp_omp
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

int main(int argc, char *argv[]) {
    printf("Hello World!");
    return 0;
}