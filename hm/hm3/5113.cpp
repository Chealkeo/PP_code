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
	
int Ceiling_log2(int x){
	unsigned temp = (unsigned)x - 1;
	int result = 0;
	
	while(temp != 0){
		temp = temp >> 1;
		result += 1;
	}
	return result;
}

int I_receive(int stage, int my_rank, int* source_ptr){
	int power_2_stage;
	
	power_2_stage = 1 << stage;
	if((my_rank >= power_2_stage) && (my_rank < 2*power_2_stage)){
		*source_ptr = my_rank - power_2_stage;
		return 1;
	}else return 0;
}	

int I_send(int stage, int my_rank, int p, int* dest_ptr){
	int power_2_stage;
	
	power_2_stage = 1 << stage;
	if(my_rank < power_2_stage){
		*dest_ptr = my_rank + power_2_stage;
		if(*dest_ptr >= p) return 0;
		else return 1;
	}else return 0;
}

void Send(double a, double b, int n, int dest){
	MPI_Send(&a, 1, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD);
	MPI_Send(&b, 1, MPI_DOUBLE, dest, 1, MPI_COMM_WORLD);
	MPI_Send(&n, 1, MPI_INT, dest, 2, MPI_COMM_WORLD);
}

void Receive(double* a_ptr, double* b_ptr, int* n_ptr, int source){
	MPI_Recv(a_ptr, 1, MPI_DOUBLE, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	MPI_Recv(b_ptr, 1, MPI_DOUBLE, source, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	MPI_Recv(n_ptr, 1, MPI_INT, source, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}	


void Get_data1(double* a_ptr, double* b_ptr, int* n_ptr, int my_rank, int p){
	int source = 0;
	int dest;
	int stage;

	
	if(my_rank == 0){
		printf("enter a, b, and n\n");
		scanf("%lf %lf %d", a_ptr, b_ptr, n_ptr);
	}
	for(stage = 0;stage < Ceiling_log2(p);stage++)
		if(I_receive(stage, my_rank, &source))
			Receive(a_ptr, b_ptr, n_ptr, source);
		else if(I_send(stage, my_rank, p, &dest))
			Send(*a_ptr, *b_ptr, *n_ptr, dest);
	
}


int main(void){
	int my_rank, n, comm_sz, local_n;
	double h, a, b, local_a, local_b;
	double local_int, total_int;
	int source;

	
	MPI_Init(NULL,NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
	
	Get_data1(&a, &b, &n, my_rank, comm_sz);

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
