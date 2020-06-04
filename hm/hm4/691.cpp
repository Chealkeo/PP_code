#include<stdio.h>
#include<mpi.h>


double f(double x){
		double return_val;
		return_val = x * x;
		return return_val;
	}
	
double Trap(double left_end, double right_end, int trap_count, double base_len){
		double estimate, x;
		int i;
		
		estimate = (f(left_end) + f(right_end))/2.0;
		for(i = 1; i <= trap_count-1; i++){
			x = left_end + i * base_len;
			estimate += f(x);
		}
		estimate *= base_len;
		
		return estimate;
	}
	
void Build_derived_type(double* a_ptr, double* b_ptr, int* n_ptr, MPI_Datatype* mpi_new_ptr){
	int block_lengths[3];
	MPI_Aint displacements[3];
	MPI_Datatype typelist[3];
	
	MPI_Aint start_address;
	MPI_Aint address;
	
	block_lengths[0] = block_lengths[1] = block_lengths[2] = 1;
	
	typelist[0]= MPI_DOUBLE;
	typelist[1]= MPI_DOUBLE;
	typelist[2]= MPI_INT;
	
	displacements[0]= 0;
	MPI_Address(a_ptr, &start_address);
	MPI_Address(b_ptr, &address);
	displacements[1]= address-start_address;
	
	MPI_Address(n_ptr, &address);
	displacements[2]= address-start_address;
	
	MPI_Type_struct(3, block_lengths, displacements, typelist, mpi_new_ptr); 
	MPI_Type_commit(mpi_new_ptr);
} 
void Get_data3(double* a_ptr, double* b_ptr, int* n_ptr, int my_rank){
	MPI_Datatype mpi_new_ptr;
	
	if(my_rank==0){
		printf("Enter a, b and n\n");
		scanf("%lf %lf %d", a_ptr, b_ptr, n_ptr);
	}
	
	Build_derived_type(a_ptr, b_ptr, n_ptr, &mpi_new_ptr);
	MPI_Bcast(a_ptr, 1, mpi_new_ptr, 0, MPI_COMM_WORLD);
	
}
/*
void Get_data(double* a_ptr, double* b_ptr, int* n_ptr, int my_rank, int p){
	int source = 0;
	int dest;
	int tag;
//	MPI_Status status;
	
	if(my_rank == 0){
		printf("enter a, b, and n\n");
		scanf("%lf %lf %d", a_ptr, b_ptr, n_ptr);
		for(dest=1; dest < p; dest++){
			tag = 0;
			MPI_Send(a_ptr, 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
			tag = 1;
			MPI_Send(b_ptr, 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
			tag = 2;
			MPI_Send(n_ptr, 1, MPI_INT, dest, tag, MPI_COMM_WORLD);
		}
	}else{
		tag = 0;
		MPI_Recv(a_ptr, 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		tag = 1;
		MPI_Recv(b_ptr, 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		tag = 2;
		MPI_Recv(n_ptr, 1, MPI_INT, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
}

*/	
int main(void){
	int my_rank, n, comm_sz, local_n;
	double h, a, b, local_a, local_b;
	double local_int, total_int;
	int source;

	
	MPI_Init(NULL,NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
	
	Get_data3(&a, &b, &n, my_rank);

	h = (b-a)/n;
	local_n = n/comm_sz;
	
	local_a = a + my_rank*h*local_n;
	local_b = local_a + h*local_n;
	local_int = Trap(local_a, local_b, local_n, h);
	
	if (my_rank != 0){
		MPI_Send(&local_int, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}else{
		total_int = local_int;
		for(source = 1;source < comm_sz; source++){
			MPI_Recv(&local_int, 1, MPI_DOUBLE, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			total_int += local_int;
		}
	}
	if(my_rank == 0){
		printf("with n = %d trapezoids, our estimate\n", n);
		printf("of the integral from %f to %f = %f\n", a, b, total_int);
		
	}
	MPI_Finalize();
	return 0;
} 
