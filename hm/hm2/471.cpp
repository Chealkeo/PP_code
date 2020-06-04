#include<stdio.h>
#include<math.h>
#include<mpi.h>


double f_1(double x){
		double return_val;
		return_val = x * x;
		return return_val;
	}

double f_2(double x){
		double return_val;
		return_val = exp(x);
		return return_val;
	}
	
double f_3(double x){
		double return_val;
		return_val = atan(x);
		return return_val;
	}
	
	
double Trap(double left_end, double right_end, int trap_count, double base_len, int c){
		double estimate, x;
		int i;
		
		if(c==1){
			estimate = (f_1(left_end) + f_1(right_end))/2.0;
			for(i = 1; i <= trap_count-1; i++){
			x = left_end + i * base_len;
			estimate += f_1(x);
		}	
		}else if(c==2){
			estimate = (f_2(left_end) + f_2(right_end))/2.0;
			for(i = 1; i <= trap_count-1; i++){
				x = left_end + i * base_len;
				estimate += f_2(x);
		}
		}
		else if(c==3){
			estimate = (f_3(left_end) + f_3(right_end))/2.0;
			for(i = 1; i <= trap_count-1; i++){
				x = left_end + i * base_len;
				estimate += f_3(x);
		}
		}
		estimate *= base_len;
		
		return estimate;
	}
	
void Get_data(double* a_ptr, double* b_ptr, int* n_ptr, int* c_ptr, int my_rank, int p){
	int source = 0;
	int dest;
	int tag;
	MPI_Status status;
	
	
	if(my_rank == 0){
		printf("please chooose the wanted function:\n");
		printf("f1:f(x)=x^2\n");
		printf("f2:f(x)=e^x\n");
		printf("f3:f(x)=actan x\n");
		printf("enter the number:c\n");
		printf("enter a, b, and n\n");
		scanf("%d %lf %lf %d",c_ptr, a_ptr, b_ptr, n_ptr);
		for(dest=1; dest < p; dest++){
			tag = 0;
			MPI_Send(a_ptr, 1, MPI_FLOAT, dest, tag, MPI_COMM_WORLD);
			tag = 1;
			MPI_Send(b_ptr, 1, MPI_FLOAT, dest, tag, MPI_COMM_WORLD);
			tag = 2;
			MPI_Send(n_ptr, 1, MPI_INT, dest, tag, MPI_COMM_WORLD);
			tag = 3;
			MPI_Send(c_ptr, 1, MPI_INT, dest, tag, MPI_COMM_WORLD);
		}
	}else{
		tag = 0;
		MPI_Recv(a_ptr, 1, MPI_FLOAT, source, tag, MPI_COMM_WORLD, &status);
		tag = 1;
		MPI_Recv(b_ptr, 1, MPI_FLOAT, source, tag, MPI_COMM_WORLD, &status);
		tag = 2;
		MPI_Recv(n_ptr, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
		tag = 3;
		MPI_Recv(c_ptr, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
	}
}

	
int main(void){
	int my_rank, n, c, comm_sz, local_n;
	double h, a, b, local_a, local_b;
	double local_int, total_int;
	int source;

	int* c_ptr;

	MPI_Status status;
	
	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
	
	Get_data(&a, &b, &n, &c, my_rank, comm_sz);
	
	h = (b-a)/n;
	local_n = n/comm_sz;
	
	local_a = a + my_rank*h*local_n;
	local_b = local_a + h*local_n;
	local_int = Trap(local_a, local_b, local_n, h, c);
	
	if (my_rank == 0){
		total_int = local_int;
		for(source = 1;source < comm_sz; source++){
			MPI_Recv(&local_int, 1, MPI_DOUBLE, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			total_int += local_int;}
		
	}else{
		MPI_Send(&local_int, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		
	}
	if(my_rank == 0){
		printf("with n = %d trapezoids, our estimate\n", n);
		printf("of the integral from %f to %f = %f\n", a, b, total_int);
		
	}
	MPI_Finalize();
	return 0;
} 