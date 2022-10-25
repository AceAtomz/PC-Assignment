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
#include <stdbool.h>
#include <string.h>
#include <limits.h>
#include <omp.h>

#define MAX_LINE_LENGTH 25
#define n_times 8
#define runs 15
double times[n_times];
FILE *fptr, *fptr2;
int V=0, n;

int* dijkstra(int** graph, int src);

// A utility function to print the constructed distance array
void printSolution(int dist[]){   
    char* currN = malloc(2);
    snprintf(currN, 2, "%d", n);
    char* fileName = malloc(50);
    strcpy(fileName, "./input_graphs/graph_");
    strcat(fileName, currN);
    strcat(fileName, ".out");
    free(currN);

    fptr = fopen(fileName,"w");
    free(fileName);
    if(fptr == NULL){
        printf("Error writing dist");
        exit(1);
    }

    for(int i=1;i<V;i++){
        fprintf(fptr,"%d %d", i, dist[i]); //print to file
        if(i!=V-1) fprintf(fptr,"\n");
    }
    fclose(fptr);
}

void validate_Dijkstra(){
    char line[MAX_LINE_LENGTH], line2[MAX_LINE_LENGTH];

    char* currN = malloc(2);
    snprintf(currN, 2, "%d", n);

    char* fileName = malloc(50);
    char* fileName2 = malloc(50);
    strcpy(fileName, "./input_graphs/graph_");
    strcpy(fileName2, "./input_graphs/graph_");
    strcat(fileName, currN);
    strcat(fileName2, currN);
    strcat(fileName, ".out");
    strcat(fileName2, "_serial.out");
    free(currN);

    fptr = fopen(fileName,"r");
    fptr2 = fopen(fileName2,"r");
    free(fileName);
    free(fileName2);
    if(fptr == NULL || fptr2 == NULL){
        printf("Error validating\n");
        exit(1);
    }

    char *throw;
    while(fgets(line, MAX_LINE_LENGTH, fptr)){
        fgets(line2, MAX_LINE_LENGTH, fptr2);
        
        char * token = strtok(line, " ");
        int vertex = strtol(token, &throw, 10);
        token = strtok(NULL, " ");
        int d = strtol(token, &throw, 10);

        char * token2 = strtok(line2, " ");
        token2 = strtok(NULL, " ");
        int vertex2 = strtol(token2, &throw, 10);
        int d2 = strtol(token2, &throw, 10);

        if(vertex==vertex2 && d!=d2){
            printf("Validata failed. \nReturned: %d\nShould be: %d\nIn graph: %d\n", d, d2, n);
            fclose(fptr);
            fclose(fptr2);
            exit(1);
        }
    }
    printf("Validate passed.\n");
    fclose(fptr);
    fclose(fptr2);
}

void printTimes(){
    fptr = fopen("sssp.out","a");

    if(fptr == NULL){
        printf("Error writing times");
        exit(1);
    }
    
    fprintf(fptr,"\n"); //newline
    for(int i=0;i<n_times;i++){
        times[i] = times[i]*1000;
        fprintf(fptr,"%f",times[i]); //print to file
        if(i!=n_times-1) fprintf(fptr,",");
    }
    fclose(fptr);
}

int** readGraph(char* fileName){
    char line[MAX_LINE_LENGTH];
    fptr = fopen(fileName,"r");
    
    if(fptr == NULL){
        printf("Error reading %s", fileName);   
        exit(1);
    }

    char *throw;
    fgets(line, MAX_LINE_LENGTH, fptr); //get the array size from V
    char * token = strtok(line, " ");
    V = strtol(token, &throw, 10);

    int** graph = (int**)malloc(V * sizeof(int*));
    for (int i = 0; i < V; i++)
        graph[i] = (int*)malloc(V * sizeof(int));
    
    for(int i=0;i<V;i++){
        for(int j=0;j<V;j++){
            graph[i][j] =0;
        }
    }

    while(fgets(line, MAX_LINE_LENGTH, fptr)){
        token = strtok(line, " ");
        int row = strtol(token, &throw, 10);//row
        token = strtok(NULL, " ");
        int col  = strtol(token, &throw, 10); //col
        token = strtok(NULL, " ");
        int weight  = strtol(token, &throw, 10); //weight

        graph[row][col] = weight;
        graph[col][row] = weight;
    }
    fclose(fptr);
    return graph;
}

int main(int argc, char *argv[]) {
    if(argc>2){
		printf("Too many input parameters\n graph.txt"); //must not have any parameters
		return 1;
	}

    if(argc==2){ //if graph provided as input
        int **graph = readGraph(argv[1]);
        dijkstra(graph, 0); //performs the dijkstra sssp operation
    }
    else{ //if no input provided, run through example graphs
        for(n=0;n<n_times;n++){
            double total_time=0, average_time=0;
            
            char* currN = malloc(2);
            snprintf(currN, 2, "%d", n);
            char* fileName = malloc(50);
            strcpy(fileName, "./input_graphs/graph_");
            strcat(fileName, currN);
            strcat(fileName, ".txt");
            free(currN);

            int **graph = readGraph(fileName);
            int *dist = (int*)malloc(V * sizeof(int));
            free(fileName);
    
            //time the scan operation
            for(int i=0;i<runs;i++){
                double start_time=omp_get_wtime(); //get start time
                
                dist = dijkstra(graph, 0); //performs the dijkstra sssp operation

                double finish_time=omp_get_wtime(); //get finish time
                total_time+= finish_time-start_time;    //add to total running time
            }
            average_time = total_time / runs;   //calc average time of algorithm
            times[n] = average_time; //store average times
            
            printSolution(dist);
            validate_Dijkstra(); //validate that the prefix sum works correctly
            free(graph);
        }
        printTimes();
    }
    return 0;
}

// Function that implements Dijkstra's single source
// shortest path algorithm for a graph represented using
// adjacency matrix representation
int* dijkstra(int **graph, int src){
    int* dist = (int*)malloc(V * sizeof(int));; // The output array.  dist[i] will hold the
                                                // shortest
                                                // distance from src to i
    int i, md, mv;
    int my_first; // The first vertex that stores in one thread locally. 
    int my_id; // ID for threads
    int my_last; //The last vertex that stores in one thread locally. 
    int my_md; // local minimum distance
    int my_mv; // local minimum vertex
    int my_step; /* local vertex that is at the minimum distance from the source */
    int nth; /* number of threads */
 
    bool sptSet[V]; // sptSet[i] will be true if vertex i is
                    // included in shortest
    // path tree or shortest distance from src to i is
    // finalized
 
    // Initialize all distances as INFINITE and stpSet[] as
    // false
    for (i = 0; i < V; i++)
        dist[i] = INT_MAX, sptSet[i] = false;
 
    // Distance of source vertex from itself is always 0
    dist[src] = 0;
    //sptSet[src] = true;

    /* OpenMP parallelization starts here */
    #pragma omp parallel private ( my_first, my_id, my_last, my_md, my_mv, my_step )\
    shared ( src, md, dist, mv, nth, graph )
    {
        my_id = omp_get_thread_num ( );
        nth = omp_get_num_threads ( );
        my_first = (my_id * V) / nth;
        my_last = ((my_id + 1) * V) / nth - 1;
        for (my_step = 0; my_step < V-1; my_step++) {
            #pragma omp single
            {
                md = INT_MAX;
                mv = -1;
            }
            int k;
            my_md = INT_MAX;
            my_mv = -1;

            /* Each thread finds the minimum distance unconnected vertex inner of
            the graph */
            for (k = my_first; k <= my_last; k++) {
                if (!sptSet[k] && dist[k] <= my_md) {
                    my_md = dist[k];
                    my_mv = k;
                }
            }

            /* 'critical' specifies that code is only be executed on 
            * one thread at a time, because we need to determine the 
            * minimum of all the my_md here. */
            #pragma omp critical
            {
                if (my_md < md) {
                    md = my_md;
                    mv = my_mv;
                }
            }

            /* 'barrier' identifies a synchronization point at which threads in a 
            * parallel region will wait until all other threads in this section reach 
            * the same point. So that md and mv have the correct value. */
            #pragma omp barrier

            #pragma omp single
            {
                /* It means we find the vertex and set its status to true. */
                if (mv != - 1){
                    sptSet[mv] = true;
                }
            }

            # pragma omp barrier

            if ( mv != -1 ){
                int j;
                for (j = my_first; j <= my_last; j++) {
                    if (!sptSet[j] && graph[mv][j] && dist[mv] != INT_MAX
                        && dist[mv] + graph[mv][j] < dist[j])
                        dist[j] = dist[mv] + graph[mv][j];
                }
            }
            #pragma omp barrier
        }
    }
    // print the constructed distance array
    return dist;
}