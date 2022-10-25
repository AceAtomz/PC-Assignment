/* File:     sssp_mpi.c
 * Purpose:  Serial implementation of Dijkstra's SSSP operation
 *
 * Input:    graph.txt file holding the arrangement of the graph
 * Output:   If the sssp was validated or not.
 *           sssp.out holds average times for plotting
 *
 * Compile:  gcc -g -Wall sssp_mpi.c -o sssp_mpi -lm
 * Run:      mpiexec -n 3 ./sssp_mpi
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
#include <mpi.h>

#define MAX_LINE_LENGTH 25
#define n_times 8
#define runs 15
double times[n_times];
FILE *fptr, *fptr2;
int V=0, n;

void SingleSource(int source, int *wgt, int *lengths, MPI_Comm comm);

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
    
    for(int i=0;i<V;i++){ //initialize all weights to INT_MAX
        for(int j=0;j<V;j++){
            graph[i][j] = INT_MAX;
            if(i==j) graph[i][j] = 0;
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
    if(argc!=1){
		printf("Too many input parameters\n"); //must not have any parameters
		return 1;
	}
     // run through example graphs
    int npes, myrank, nlocal;
    int i, j, k;
    int *localWeight; /*local weight array*/
    int *localDistance; /*local distance vector*/

    MPI_Init(&argc, &argv);
    for(n=0;n<n_times;n++){
        double total_time=0, average_time=0;
        int *dist;
        char* currN = malloc(2);
        snprintf(currN, 2, "%d", n);
        char* fileName = malloc(50);
        strcpy(fileName, "./input_graphs/graph_");
        strcat(fileName, currN);
        strcat(fileName, ".txt");
        if (myrank == 0) free(currN);

        int **graph = readGraph(fileName); /*adjacency matrix*/
        if (myrank == 0) free(fileName);
        
        int sendbuf[V*V]; /*local weight to distribute*/

        MPI_Comm_size(MPI_COMM_WORLD, &npes);
        MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

        nlocal = V/npes; /* Compute the number of elements to be stored locally per thread. */
        //time the scan operation
        for(int r=0;r<runs;r++){
            dist = (int*)malloc(V * sizeof(int)); /*distance vector*/
            localWeight = (int *)malloc(nlocal*V*sizeof(int));
            localDistance = (int *)malloc(nlocal*sizeof(int));

            

            /*prepare send data */
            for(k=0; k<npes; ++k) {
                for(i=0; i<V;++i) {
                    for(j=0; j<nlocal;++j) {
                        sendbuf[k*V*nlocal+i*nlocal+j]=graph[i][k*nlocal+j];
                    }
                }
            }

            /*distribute data*/
            MPI_Scatter(sendbuf, nlocal*V, MPI_INT, localWeight, nlocal*V, MPI_INT,
            0, MPI_COMM_WORLD); 
            
            double start_time=MPI_Wtime(); //get start time

            /*Implement the single source dijkstra's algorithm*/
            SingleSource(0, localWeight, localDistance, MPI_COMM_WORLD);

            double finish_time=MPI_Wtime(); //get finish time   
            
            /*collect local distance vector at the source process*/
            MPI_Gather(localDistance, nlocal, MPI_INT, dist, nlocal, MPI_INT, 
            0, MPI_COMM_WORLD);

             
            if (myrank == 0) {
                total_time+= finish_time-start_time;    //add to total running time
                free(localWeight);
                free(localDistance);
                if(r!=runs-1) free(dist);
            }
        }
        if (myrank == 0){
            average_time = total_time / runs;   //calc average time of algorithm
            times[n] = average_time; //store average times

            printSolution(dist);
            validate_Dijkstra(); //validate that the prefix sum works correctly
            free(graph);
        } 
    }
    MPI_Finalize();
    if (myrank == 0) printTimes();
    return 0;
}


/*single source Dijkstra's Algorithm*/
/* @param source: rank of the root
 @param wgt: points to locally stored portion of the weight adjacency matrix of the graph;
 @param lengths: points to a vector that will store the distance of the shortest paths from the
source to the locally stored vertices;
*/
void SingleSource(int source, int *wgt, int *lengths, MPI_Comm comm) {
    int i, j;
    int nlocal; /* The number of vertices stored locally */
    int *marker; /* Used to mark the vertices belonging to Vo */
    int firstvtx; /* The index number of the first vertex that is stored locally */
    int lastvtx; /* The index number of the last vertex that is stored locally */
    int u, udist;
    int lminpair[2], gminpair[2];
    int npes, myrank;

    MPI_Comm_size(comm, &npes);
    MPI_Comm_rank(comm, &myrank);

    nlocal = V / npes;
    firstvtx = myrank*nlocal;
    lastvtx = firstvtx + nlocal - 1;

    /* Set the initial distances from source to all the other vertices */
    for (j = 0; j<nlocal; j++) {
        lengths[j] = wgt[source*nlocal + j];
    }
    /* This array is used to indicate if the shortest part to a vertex has been found or not. */
    /* if marker [v] is one, then the shortest path to v has been found, otherwise it has not. */
    marker = (int *)malloc(nlocal*sizeof(int));
    for (j = 0; j<nlocal; j++) {
        marker[j] = 1;
    }

    /* The process that stores the source vertex, marks it as being seen */
    if (source >= firstvtx && source <= lastvtx) {
        marker[source - firstvtx] = 0;
        //lengths[source] = 0;
    }

    /* The main loop of Dijkstra's algorithm */
    for (i = 1; i<V; i++) {
        /* Step 1: Find the local vertex that is at the smallest distance from source */
        lminpair[0] = INT_MAX; /* set it to an architecture dependent large number */
        lminpair[1] = -1;
        for (j = 0; j<nlocal; j++) {
            if (marker[j] && lengths[j] < lminpair[0]) {
                lminpair[0] = lengths[j];
                lminpair[1] = firstvtx + j;
            }
        }

        /* Step 2: Compute the global minimum vertex, and insert it into Vc */
        MPI_Allreduce(lminpair, gminpair, 1, MPI_2INT, MPI_MINLOC, comm);
        udist = gminpair[0];
        u = gminpair[1];
        /* The process that stores the minimum vertex, marks it as being seen */
        if (u == lminpair[1]) {
            marker[u - firstvtx] = 0;
        }

        /* Step 3: Update the distances given that u got inserted */
        for (j = 0; j<nlocal; j++) {
            if (marker[j]  && wgt[u*nlocal + j] != INT_MAX && udist != INT_MAX
                && ((udist + wgt[u*nlocal + j]) < lengths[j])) {
                lengths[j] = udist + wgt[u*nlocal + j];
            }
        }
    }
    if (myrank == 0) free(marker);
}