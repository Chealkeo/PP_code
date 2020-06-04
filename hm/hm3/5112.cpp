#include <stdio.h>
#include <mpi.h>

int main(void){
	int x, y, z;
	int my_rank, comm_sz; 
	MPI_Status status;

	MPI_Init(NULL,NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
	
	switch(my_rank){
		case 0: x=0; y=1; z=2;
			MPI_Bcast(&x, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Send(&y, 1, MPI_INT, 2, 43, MPI_COMM_WORLD);
			MPI_Bcast(&z, 1, MPI_INT, 1, MPI_COMM_WORLD);
			break;
		case 1: x=3; y=4; z=5;
			MPI_Bcast(&x, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast(&y, 1, MPI_INT, 1, MPI_COMM_WORLD);
			break;
		case 2: x=6; y=7; z=8;
			MPI_Bcast(&z, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Recv(&x, 1, MPI_INT, 0, 43, MPI_COMM_WORLD, &status);
			MPI_Bcast(&y, 1, MPI_INT, 1, MPI_COMM_WORLD);
			break;
				
	}
	
	for(; my_rank <3; my_rank ++){
		printf("In process %d, x is %d, y is %d, z is %d;\n", my_rank, x, y, z);
	MPI_Finalize();
	return 0;
}
}

