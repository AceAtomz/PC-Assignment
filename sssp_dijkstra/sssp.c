/* File:     sssp.c
 * Purpose:  Serial implementation of Dijkstra's SSSP operation
 *
 * Input:    graph.txt file holding the arrangement of the graph
 * Output:   If sssp was validated or not.
 *           sssp.out holds average times for plotting
 *
 * Compile:  gcc -g -Wall -fopen sssp.c -o sssp
 * Run:      ./sssp
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

int minDistance(int *dist, bool sptSet[]);
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
            printf("Validata failed. \n%d \n%d", d, d2);
        }
    }
    printf("Validate passed.\n");
    fclose(fptr);
    fclose(fptr2);
}

void printTimes(){
    fptr = fopen("sssp.out","w");

    if(fptr == NULL){
        printf("Error writing times");
        exit(1);
    }
    
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

// A utility function to find the vertex with minimum
// distance value, from the set of vertices not yet included
// in shortest path tree
int minDistance(int *dist, bool sptSet[]){
    // Initialize min value
    int min = INT_MAX, min_index;
 
    for (int v = 0; v < V; v++)
        if (sptSet[v] == false && dist[v] <= min)
            min = dist[v], min_index = v;
 
    return min_index;
}

// Function that implements Dijkstra's single source
// shortest path algorithm for a graph represented using
// adjacency matrix representation
int* dijkstra(int **graph, int src){
    int* dist = (int*)malloc(V * sizeof(int));; // The output array.  dist[i] will hold the
                 // shortest
    // distance from src to i
 
    bool sptSet[V]; // sptSet[i] will be true if vertex i is
                    // included in shortest
    // path tree or shortest distance from src to i is
    // finalized
 
    // Initialize all distances as INFINITE and stpSet[] as
    // false
    for (int i = 0; i < V; i++)
        dist[i] = INT_MAX, sptSet[i] = false;
 
    // Distance of source vertex from itself is always 0
    dist[src] = 0;
 
    // Find shortest path for all vertices
    for (int count = 0; count < V - 1; count++) {
        // Pick the minimum distance vertex from the set of
        // vertices not yet processed. u is always equal to
        // src in the first iteration.
        int u = minDistance(dist, sptSet);
 
        // Mark the picked vertex as processed
        sptSet[u] = true;
 
        // Update dist value of the adjacent vertices of the
        // picked vertex.
        for (int v = 0; v < V; v++)
 
            // Update dist[v] only if is not in sptSet,
            // there is an edge from u to v, and total
            // weight of path from src to  v through u is
            // smaller than current value of dist[v]
            if (!sptSet[v] && graph[u][v]
                && dist[u] != INT_MAX
                && dist[u] + graph[u][v] < dist[v])
                dist[v] = dist[u] + graph[u][v];
    }
 
    // print the constructed distance array
    return dist;
}