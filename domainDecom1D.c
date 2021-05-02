
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char *argv[]){

    int rank, size, i, provided;
    
    // number of cells (global)
    int nxc = 128; // make sure nxc is divisible by size
    double L = 2 * 3.1415; // Length of the domain
    
    MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // number of nodes (local to the process): 0 and nxn_loc-1 are ghost cells 
    int nxn_loc = nxc/size + 3; // number of nodes is number cells + 1; we add also 2 ghost cells
    double L_loc = L/((double) size);
    double dx = L / ((double) nxc);
    
    // define out function
    double *f = calloc(nxn_loc, sizeof(double)); // allocate and fill with z
    double *dfdx = calloc(nxn_loc, sizeof(double)); // allocate and fill with z
    for (i=1; i<(nxn_loc-1); i++)
      f[i] = sin(L_loc*rank + (i-1) * dx);

    // need to communicate and fill ghost cells f[0] and f[nxn_loc-1]
    // communicate ghost cells

		int next_process = (size + rank + 1) % size;
		int prev_process = (size + rank - 1) % size;
		
		// Send first ghost cell back to previous process
		MPI_Send(&f[2], 1, MPI_DOUBLE, prev_process, 0, MPI_COMM_WORLD);		
		// Recevive last ghost cell from next process
		MPI_Recv(&f[nxn_loc-1], 1, MPI_DOUBLE, next_process, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 

		// Send last ghost cell forward to next process
		MPI_Send(&f[nxn_loc-3], 1, MPI_DOUBLE, next_process, 0, MPI_COMM_WORLD);
		// Receive last ghost cell from previous process
		MPI_Recv(&f[0], 1, MPI_DOUBLE, prev_process, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		
    // calculate first order derivative using central difference
    // here we need to correct value of the ghost cells!
    for (i=1; i<(nxn_loc-1); i++)
      dfdx[i] = (f[i+1] - f[i-1]) / (2.0*dx);
    
    // Print f values
    if (rank==0){ // print only rank 0 for convenience
        printf("My rank %d of %d\n", rank, size );
        printf("Here are my values for f including ghost cells\n");
        for (i=0; i<nxn_loc; i++) {
					if (i == 0 || i == nxn_loc - 1) {
	      		printf("f (sin): %f, dfdx: NAN, cos: %f\n", f[i],  cos(L_loc*rank + (i-1) * dx));
					} else {
						double real_cos = cos(L_loc*rank + (i-1) * dx);
	      		printf("f (sin): %f, dfdx: %f, cos: %f, diff: %f\n", f[i], dfdx[i], real_cos, (dfdx[i] - real_cos));
					}
				}
				printf("\n");
    }   

    MPI_Finalize();
}






