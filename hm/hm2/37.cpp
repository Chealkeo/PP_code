#include<stdio.h>
#include<string.h>
#include<mpi.h>

const int MAX_STRING = 100;

int main(void){
	char greeting [MAX_STRING];
	int comm_sz;
	int my_rank;
	
	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	
	for(;my_rank <= comm_sz-1; my_rank++) {
		sprintf(greeting, "Greetings from process %d of %d!", my_rank, comm_sz);
		MPI_Send(greeting, strlen(greeting)+1, MPI_CHAR, (my_rank+1)%comm_sz, 0, MPI_COMM_WORLD);
		MPI_Recv(greeting, MAX_STRING, MPI_CHAR, (my_rank+comm_sz-1)%comm_sz, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		printf("%s\n", greeting);
	}
	
	
	MPI_Finalize();
	return 0;
}