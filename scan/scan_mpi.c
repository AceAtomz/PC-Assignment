/* File:     scan_omp.c
 * Purpose:  Parallel implementation of the Prefix Sum operation (Blelloch 1990)
 *           implemented with OpenMP
 *
 * Input:    
 * Output:   If the scan was validated or not.
 *           scan.out holds average times for plotting
 *
 * Compile:  gcc -g -Wall -fopen scan_omp.c -o scan_omp
 * Run:      ./scan_omp
 *
 * Algorithm:
 *    1.  Given a set of elements [a0,a1,··· ,an−1].
 *    2.  The scan operation returns the set [a0,(a0+a1),···,(a0+a1 +···+an−1)].
 *    eg. The input set [2,1,4,0,3,7,6,3],
 *        yields [2,3,7,7,10,17,23,26].
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define size 1000
#define maxsize 50000
#define n_times (maxsize/size)
#define runs 10
double times[n_times];
FILE *fptr;

int find_sum(int* numbers, size_t mysize);
void start_find_sum(int rank, int mysize, int* in,
                    size_t num_per_proc, int* overall_sum);

void start_find_psum(int rank, int mysize, int* in, 
                     size_t num_per_proc, int sum);

void printTimes(){
    fptr = fopen("scan.out","a");

    if(fptr == NULL){
        printf("Error!");   
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
    int my_rank, comm_size;
    
    if(argc != 1){
		printf("Too many input parameters\n"); //must not have any parameters
		return 1;
	}

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size); 
    MPI_Barrier(MPI_COMM_WORLD);
    

    for(int n=size;n<=maxsize;n+=size){
        double total_time=0.0, average_time=0.0;
        srand(n); //sets seed for random generation

        int *data=(int*)malloc(sizeof(int)*n);
        int *outData=(int*)malloc(sizeof(int)*n);
        for(int i=0; i<n; i++){     //generates sudo-random numbers
            data[i]=1.1*rand()*5000/RAND_MAX;
        }
        
        //time the scan operation
        for(int i=0;i<runs;i++){
            MPI_Barrier(MPI_COMM_WORLD);

            size_t num_per_proc = n / comm_size;
            int sum;

            double start_time=MPI_Wtime(); //get start time
            start_find_sum(my_rank, comm_size, data, num_per_proc,
                           &sum);
            start_find_psum(my_rank, comm_size, data, num_per_proc,
                            sum);
            
            MPI_Barrier(MPI_COMM_WORLD); //ensures all processes have finished computing prefix sums

            double finish_time=MPI_Wtime(); //get finish time
            
            if (my_rank == 0) {
                total_time+= finish_time-start_time;    //add to total running time
            }
            
        }
        if (my_rank == 0) {
            average_time = total_time / runs;   //calc average time of algorithm
            times[times_counter] = average_time; //store average times
            
            if(n==maxsize) validate_scan(n, outData); //validate that the prefix sum works correctly
            free(data);
            free(outData);
            times_counter++;
        }
    }
    MPI_Finalize();
    if (my_rank == 0) printTimes();
    return 0;
}

// returns the sum of all the elements in `numbers`, an array of size `size`
int find_sum(int* numbers, size_t mysize){
  int sum = 0;
  for (size_t i = 0; i < mysize; i++) {
    sum += numbers[i];
  }
  return sum;
}

/*
 * rank - the current process
 * size - the total number of processes
 * in - the set of numbers of the current process to find the max of.
 * num_per_proc - the length of `in` array
 * overall_sum - the sum of this rank will be written
 *
 * Find the max among all processes, each having a different
 * `in`
 */
void start_find_sum(int rank, int mysize, int* in, 
                    size_t num_per_proc, int* overall_sum){
  MPI_Status status;

  int sum = find_sum(in, num_per_proc);

  int still_alive = 1;
  int level;

  /* `level` is the current level of the complete binary
   * tree. At level 0, there are n nodes (processes); each
   * of these nodes has its sum already, computed by
   * find_sum.  A node with label `rank` that is in an
   * even-numbered position on the current level becomes a
   * parent on the next level; its new sum is the sum of
   * its current sum and the sum of node `rank + 2^level`.
   * A node with label `rank` in an odd-numbered position
   * on the current level is responsible for sending its
   * sum to the parent, `rank - 2^level`; after that, the
   * node becomes inactive (still_alive is set to 0).
   * After log2(size) levels have been created, node with
   * rank 0 contains the sum.
   */

  for (level = 0; level < (int)log2(mysize); level++) {
    if (still_alive) {
      int position = rank / (int)pow(2, level);

      if (position % 2 == 0) {
        // I am a receiver
        int sender_sum;
        int sending_rank = rank + (int)pow(2, level);

        MPI_Recv(&sender_sum, 1, MPI_INT, sending_rank,
                 0, MPI_COMM_WORLD, &status);

        sum += sender_sum;
      }

      else {
        // I am a sender
        int receiving_rank = rank - (int)pow(2, level);

        MPI_Send(&sum, 1, MPI_INT, receiving_rank, 0,
                 MPI_COMM_WORLD);
        still_alive = 0;
      }
    }

    MPI_Barrier(MPI_COMM_WORLD);
  }

  *overall_sum = sum;
}

/*
 *             Description of start_find_psum
 * A process with rank R on level i that is an even-numbered
 * position as its sum set to the sum of this process R on
 * level i + 1, rather than level i. This sum was created by
 * adding the old sum from level i with the sum of process
 * R + 2^i (the sibling of process R). The correct sum for
 * process R on level i can be restored by subtracting away
 * the sum of the sibling process of R.
 *
 * (The above about restoring the correct sum isn't necessary and is
 * commented out.)
 *
 * Each node (process) in the tree has an associated prefix
 * sum that represents the prefix sum of all the numbers up
 * to the last number of the rightmost leaf process in the
 * tree rooted at that node. Rank 0 is the root node and so
 * its psum is set to sum.
 *
 * Iteration is from the level under the root to the bottom
 * level. A node in an odd-numbered position on the current
 * level receives its prefix sum from its parent.  Because
 * even-numbered nodes become parents on the next level, the
 * parent of any odd numbered node of rank R on level i is R
 * - 2^i. A node in an even-numbered position on level i has
 * its prefix sum set as the prefix sum of its parent minus
 * the regular sum of its sibling.  Because a node in an
 * even-numbered position is also its own parent, it already
 * has the prefix sum of its parent. It also has access to
 * the regular sum of its sibling from the time it fixed its
 * own sum, so that value can be reused.
 *
 * Now every node of rank R has its prefix sum; this
 * prefix sum is the sum of all of the `random_numbers` of
 * nodes from Rank 0 to Rank R. We will overwrite
 * `random_numbers` to be the prefix sums for
 * num_per_proc*rank  TO   num_per_proc*(rank + 1) - 1,
 * from a global perspective, if we consider each
 * random_numbers as part of a distributed array.
 */
void start_find_psum(int rank, int mysize, int* in,
                    size_t num_per_proc, int sum){
  int psum;
  int level;
  MPI_Status status;

  if (rank == 0) {
    psum = sum;
  }
  for (level = (int)log2(mysize) - 1; level >= 0; level--) {

    // only trigger the processes on the current level
    if (level == 0 || rank % (int)pow(2, level) == 0) {
      int position = rank / (int)pow(2, level);

      if (position % 2 == 0) {
        int sender_sum;
        int sending_rank = rank + (int)pow(2, level);

        // in the previous level, this process was the parent
        // of `sending_rank`, and had it as its right child,
        // so this psum value is the psum value of the parent
        // of our sibling
        MPI_Send(&psum, 1, MPI_INT,
                 sending_rank, // RIGHT CHILD
                 0, MPI_COMM_WORLD);

        MPI_Recv(&sender_sum, 1, MPI_INT,
                 sending_rank, 0, MPI_COMM_WORLD, &status);

        // fix the sum to be the correct value:
	// this isn't required, so it's commented out
        //sum -= sender_sum;

        // psum <- (prefix sum of parent) - (sum of sibling)
        psum -= sender_sum;
      }

      else{
        int receiving_rank = rank - (int)pow(2, level);

        MPI_Recv(&psum, 1, MPI_INT,
                 receiving_rank, // PARENT
                 0, MPI_COMM_WORLD, &status);

        // send sum to receiving_rank so it can fix its sum
        MPI_Send(&sum, 1, MPI_INT,
                 receiving_rank, 0, MPI_COMM_WORLD);
      }

    }

    MPI_Barrier(MPI_COMM_WORLD);
  }

  // put the prefix sums associated with this node in random_numbers
  int next_sum = in[num_per_proc-1];
  in[num_per_proc-1] = psum;

  for (int j = num_per_proc - 2; j >= 0; j--) {
    int next_sum_tmp = in[j];

    in[j] = in[j+1] - next_sum;

    next_sum = next_sum_tmp;
  }
}