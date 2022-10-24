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

#define n_times 8
#define runs 15
double times[n_times];
FILE *fptr;
int V=0;

int minDistance(int dist[], bool sptSet[]);
void dijkstra(int graph[V][V], int src);

// A utility function to print the constructed distance array
void printSolution(int dist[])
{   
    printf("The distance from 0 to each vertex is:\nv \t\t dist 0->v\n");
    for (int i = 0; i < V; i++)
        printf("%d \t\t\t %d\n", i, dist[i]);
}

void printTimes(){
    fptr = fopen("sssp.out","w");

    if(fptr == NULL){
        printf("Error!");   
        exit(1);
    }
    
    for(int i=0;i<n_times;i++){
        times[i] = times[i]*1000;
        fprintf(fptr,"%f",times[i]); //print to file
        if(i!=n_times-1) fprintf(fptr,",");
    }
    fclose(fptr);
}

void readGraph(char* fileName){
    char *text;
    long numbytes;
    fptr = fopen(fileName,"r");

    if(fptr == NULL){
        printf("Error!");   
        exit(1);
    }

    fseek(fptr, 0L, SEEK_END);
    numbytes = ftell(fptr);
    fseek(fptr, 0L, SEEK_SET); 

    text = (char*)calloc(numbytes, sizeof(char));   
    if(text == NULL){
        printf("Error reading from file");
        exit(1);
    }

    fread(text, sizeof(char), numbytes, fptr);
    fclose(fptr);
    printf("%s",text);
}

int main(int argc, char *argv[]) {
    int times_counter=0;

    if(argc>2){
		printf("Too many input parameters\n graph.txt"); //must not have any parameters
		return 1;
	}

    if(argc==2){ //if graph provided as input
        printf("Hello World!");
        return 0;
    }
    else{ //if no input provided, run through example graphs
        for(int n=0;n<n_times;n++){
            double total_time=0, average_time=0;
            
            char* currN = malloc(2);
            snprintf(currN, 2, "%d", n);
            char* fileName = malloc(50);
            strcat(fileName, "./input_graphs/graph_");
            strcat(fileName, currN);
            strcat(fileName, ".txt");
            free(currN);

            readGraph(fileName); //int *graph = 
            free(fileName);
            
            // int *graph[V][V] = malloc((V * V) * sizeof(int));
            // int *data=(int*)malloc(sizeof(int)*n);
            // for(int i=0; i<n; i++){     //generates sudo-random numbers
            //     data[i]=1.1*rand()*5000/RAND_MAX;
            // }
            
            // //time the scan operation
            // for(int i=0;i<runs;i++){
            //     double start_time=omp_get_wtime(); //get start time
                
            //     scan(outData, data, n); //performs the prefix sum operation

            //     double finish_time=omp_get_wtime(); //get finish time
            //     total_time+= finish_time-start_time;    //add to total running time
            // }
            // average_time = total_time / runs;   //calc average time of algorithm
            // times[times_counter] = average_time; //store average times
            
            // if(n==maxsize) validate_scan(n, outData); //validate that the prefix sum works correctly
            // free(data);
            // free(outData);
            // times_counter++;
        }

        printTimes();
    }
    return 0;
}

// A utility function to find the vertex with minimum
// distance value, from the set of vertices not yet included
// in shortest path tree
int minDistance(int dist[], bool sptSet[]){
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
void dijkstra(int graph[V][V], int src){
    int dist[V]; // The output array.  dist[i] will hold the
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
    printSolution(dist);
}