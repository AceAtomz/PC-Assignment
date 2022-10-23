/* File:     scan.c
 * Purpose:  Serial implementation of the Scan operation/All-prefix-sums operation (Blelloch 1990)
 *
 * Input:    
 * Output:   Estimate of the integral from a to b of f(x)
 *           using the trapezoidal rule and n trapezoids.
 *
 * Compile:  gcc -g -Wall -fopen scan.c -o scan
 * Run:      ./scan
 *
 * Algorithm: (See Lec 8 slides)
 *    1.  Each process calculates "its" interval of integration.
 *    2.  Each process estimates the integral of f(x)
 *        over its interval using the trapezoidal rule.
 *    3a. Each process != 0 sends its integral to 0.
 *    3b. Process 0 sums the calculations received from
 *        the individual processes and prints the result.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

int size = 1000;
int maxsize = 20000;
double times[20];
int runs = 5;

void scan(int out[], int in[], int N){
    out[0] = in[0];

    for(int i=1; i<N; i++) {
        out[i] = in[i] + out[i-1];
    }
}

int main(int argc, char *argv[]) {
    int times_counter=0;
    
    if(argc != 1){
		printf("Wrong input parameters\nscan.out\n"); //must specify .out file
		return 1;
	}

    for(int n=size;n<=maxsize;n+=size){
        double total_time=0, average_time=0;
        srand(n); //sets seed for random generation

        int data[n];
        for(int i=0; i<n; i++){     //generates sudo-random numbers
            data[i]=1.1*rand()*5000/RAND_MAX;
        }
        
        //time the scan operation
        for(int i=0;i<runs;i++){
            double start_time=omp_get_wtime(); //get start time
            
            double finish_time=omp_get_wtime(); //get finish time
            total_time+= finish_time-start_time;    //add to total running time
        }
        average_time = total_time / runs;   //calc average time of algorithm
        times[times_counter] = average_time;

        //scan();
        times_counter++;
    }
}