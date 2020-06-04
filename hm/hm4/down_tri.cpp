#include<stdio.h>
#include<mpi.h>
#include<stdlib.h>

int main(void){
	int n=10;
	int my_rank,comm_sz, i, j;
	double A[n][n];
	double T[n][n];
	int displacements[n];
	int block_lengths[n];
	MPI_Datatype index_mpi_t;
	MPI_Status status;
	
	MPI_Init(NULL,NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
	
	for(i = 0;i<n;i++){
		block_lengths[i]= i+1;
		displacements[i]= n*i;
		
	}
	MPI_Type_indexed(n, block_lengths, displacements, MPI_DOUBLE, &index_mpi_t);
	MPI_Type_commit(&index_mpi_t);
	
	if(my_rank == 0){
		for(i=0;i<n;i++){
			for(j=0;j<n;j++){
				A[i][j]=rand()%10;
			}
		}
		MPI_Send(A, 1, index_mpi_t, 1,0, MPI_COMM_WORLD);
	}else{
		MPI_Recv(T, 1, index_mpi_t, 0,0, MPI_COMM_WORLD, &status);
		for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++)
                printf("%4.1f ", T[i][j]);
            printf("\n");}
	
	}
		
		
	MPI_Finalize();
	return 0;
}