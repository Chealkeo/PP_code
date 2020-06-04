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
	

	
int main(void){
	int my_rank, comm_sz, n=1024, local_n;
	double a = 0.0, b = 1.0, h, local_a, local_b;
	double local_int, total_int;
	int source;
	
	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
	
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
		printf("of the integral from %f to %f = %.15e\n", a, b, total_int);
		
	}
	MPI_Finalize();
	return 0;
} 