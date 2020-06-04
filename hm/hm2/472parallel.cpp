#include<stdio.h>
#include<mpi.h>


double f(double x){
		double return_val;
		return_val = x * x;
		return return_val;
	}
	
double Simpson(double left_end, double right_end, int n){
		double estimate, x, h;
		int i;
		
		estimate = f(left_end) + f(right_end);
		h = (right_end-left_end)/n;
		x = left_end;
		for(i=1;i < n; i++){
			x += h;
			if(i%2==1){
				estimate += 4*f(x);
			}else{
				estimate += 2*f(x);
			}
	}
		estimate *= h/3;
		
		return estimate;
	}
	
void Get_data(double* a_ptr, double* b_ptr, int* n_ptr, int my_rank, int p){
	int source = 0;
	int dest;
	int tag;

	
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

	
int main(void){
	int my_rank, n, comm_sz, local_n;
	double h, a, b, local_a, local_b;
	double local_int, total_int;
	int source;

	
	MPI_Init(NULL,NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
	
	Get_data(&a, &b, &n, my_rank, comm_sz);

	h = (b-a)/n;
	local_n = n/comm_sz;
	
	local_a = a + my_rank*h*local_n;
	local_b = local_a + h*local_n;
	local_int = Simpson(local_a, local_b, local_n);
	
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
