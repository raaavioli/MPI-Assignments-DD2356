
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <mpi.h>

#define SEED     921
#define NUM_ITER 1000000000

int main(int argc, char* argv[])
{
    int size, rank, provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int count = 0;
    double x, y, z, pi, start;

    srand(SEED * rank); // Important: Multiply SEED by "rank" when you introduce MPI!
    
    int max_iter = NUM_ITER / size;
    // Calculate PI following a Monte Carlo method
    for (int iter = 0; iter < max_iter; iter++)
    {
        // Generate random (X,Y) points
        x = (double)random() / (double)RAND_MAX;
        y = (double)random() / (double)RAND_MAX;
        z = sqrt((x*x) + (y*y));
        
        // Check if point is in unit circle
        if (z <= 1.0)
        {
            count++;
        }
    }

    // Warm up
    if (rank != 0) {
        MPI_Send(&count, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    } else {
	for (int i = 1; i < size; i++) {
            int o_count;
	    MPI_Recv(&o_count, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
    }
    // End warmup

    start = MPI_Wtime();
    if (rank != 0) {
        MPI_Send(&count, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    } else {
	for (int i = 1; i < size; i++) {
            int o_count;
	    MPI_Recv(&o_count, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    count += o_count;
	}
        // Estimate Pi and display the result
        pi = ((double)count / (double)(max_iter * size)) * 4.0;
    
	double time = MPI_Wtime() - start;    	
        printf("The result is %f\n", pi);
	printf("Execution time: %f\n", time);
    }    

    MPI_Finalize();
    return 0;
}

