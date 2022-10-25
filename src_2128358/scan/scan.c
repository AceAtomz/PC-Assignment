/* File:     scan.c
 * Purpose:  Serial implementation of the Prefix Sum operation (Blelloch 1990)
 *
 * Input:    
 * Output:   If the scan was validated or not.
 *           scan.out holds average times for plotting
 *
 * Compile:  gcc -g -Wall -fopen scan.c -o scan
 * Run:      ./scan
 *
 * Algorithm:
 *    1.  Given a set of elements [a0,a1,··· ,an−1].
 *    2.  The scan operation returns the set [a0,(a0+a1),···,(a0+a1 +···+an−1)].
 *    eg. The input set [2,1,4,0,3,7,6,3],
 *        yields [2,3,7,7,10,17,23,26].
 */

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#define size 1000
#define maxsize 50000
#define n_times (maxsize/size)
#define runs 5
double times[n_times];
FILE *fptr;

void scan(int out[], int in[], int N){
    out[0] = in[0];

    for(int i=1; i<N; i++) {
        out[i] = in[i] + out[i-1];
    }
}

void printTimes(){
    fptr = fopen("scan.out","w");

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

void validate_scan(int n, int *data){
	int i;
	for(i=0;i<n-1;i++){
		if(data[i] > data[i+1]){
			printf("Validata failed. \n");
		}
	}
	printf("Validate passed.\n");
}

int main(int argc, char *argv[]) {
    int times_counter=0;
    
    if(argc != 1){
		printf("Too many input parameters\n"); //must not have any parameters
		return 1;
	}

    for(int n=size;n<=maxsize;n+=size){
        double total_time=0, average_time=0;
        srand(n); //sets seed for random generation

        int *data=(int*)malloc(sizeof(int)*n);
        int *outData=(int*)malloc(sizeof(int)*n);
        for(int i=0; i<n; i++){     //generates sudo-random numbers
            data[i]=1.1*rand()*5000/RAND_MAX;
        }
        
        //time the scan operation
        for(int i=0;i<runs;i++){
            double start_time=omp_get_wtime(); //get start time
            
            scan(outData, data, n); //performs the prefix sum operation

            double finish_time=omp_get_wtime(); //get finish time
            total_time+= finish_time-start_time;    //add to total running time
        }
        average_time = total_time / runs;   //calc average time of algorithm
        times[times_counter] = average_time; //store average times
        
        if(n==maxsize) validate_scan(n, outData); //validate that the prefix sum works correctly
        free(data);
        free(outData);
        times_counter++;
    }

    printTimes();
}