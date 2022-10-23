/* File:     scan.c
 * Purpose:  Serial implementation of the Scan operation/All-prefix-sums operation (Blelloch 1990)
 *
 * Input:    
 * Output:   Estimate of the integral from a to b of f(x)
 *           using the trapezoidal rule and n trapezoids.
 *
 * Compile:  gcc scan.c -o scan
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

int size = 1000;
int maxsize = 20000;
int times[20];
int runs = 5;

void scan(int out[], int in[], int N){
    out[0] = in[0];

    for(int i=1; i<N; i++) {
        out[i] = in[i] + out[i-1];
    }
}

int main(int argc, char *argv[]) {
    double start, s_time=0, p_time=0;
    
    if(argc != 1){
		printf("Wrong input parameters\nscan.out\n"); //must specify .out file
		return 1;
	}

    for(int n=size;n<=maxsize;n+=size){
        srand(n); //sets seed for random generation

        int data[n];
        for(int i=0; i<n; i++){     //generates sudo-random numbers
            data[i]=1.1*rand()*5000/RAND_MAX;
        }
        

        //scan();
    }
}