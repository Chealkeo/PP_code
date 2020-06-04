#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>

double f(double x)
{
   double return_val;
   return_val = x*x;
   return return_val;
}

double Local_trap(double a, double b, int n) {
   double  h, x, my_result;
   double  local_a, local_b;
   int  i, local_n;
   int my_rank = omp_get_thread_num();
   int thread_count = omp_get_num_threads();

   h = (b-a)/n;
   local_n = n/thread_count;
   local_a = a + my_rank*local_n*h;
   local_b = local_a + local_n*h;
   my_result = (f(local_a) + f(local_b))/2.0;
   for (i = 1; i <= local_n-1; i++) {
     x = local_a + i*h;
     my_result += f(x);
   }
   my_result = my_result*h;

   return my_result;
}

void Trap(double a, double b, int n, double* global_result_p);

int main(int argc, char* argv[])
{
	double  global_result = 0.0;  
	double  startTime,finishTime;
	
	double a = strtol(argv[1],NULL, 10);
	double b = strtol(argv[2],NULL, 10);
	int n = strtol(argv[3],NULL, 10);
	int version = strtol(argv[4],NULL, 10);
	int threadCount = strtol(argv[5],NULL, 10);
//	double  a,b;
//	int n;
//	printf("enter a, b and n\n");
//	scanf("%lf %lf %d", &a, &b, &n);
	
	if(version==1)
	{
		startTime = omp_get_wtime();
		#  pragma omp parallel num_threads(threadCount)
		{
		   Trap(a, b, n, &global_result);
		}
		finishTime = omp_get_wtime();
		printf("ThreadCount: %d,versionchoose: %d\n", threadCount,version);
		printf("With n = %d trapezoids, our estimate\n", n);
		printf("of the integral from %f to %f = %f\n",a, b, global_result);
		printf("Time used in parallel block = %e seconds\n", finishTime-startTime);
	}
	else if(version==2)
	{
		startTime = omp_get_wtime();
		#  pragma omp parallel num_threads(threadCount)
		   {
		      double my_result = 0.0;
		      my_result += Local_trap(a,b, n);
		#     pragma omp critical
		      global_result += my_result;
		   }
  		finishTime = omp_get_wtime();
		printf("ThreadCount: %d,versionchoose: %d\n", threadCount,version);
		printf("With n = %d trapezoids, our estimate\n", n);
		printf("of the integral from %f to %f = %f\n",a, b, global_result);
		printf("Time used in parallel block = %e seconds\n", finishTime-startTime);
	}
	else
	{
		startTime = omp_get_wtime();
  		#  pragma omp parallel num_threads(threadCount) \
		      reduction(+: global_result)
		   {
		      global_result += Local_trap(a,b,n);
		   }
		finishTime = omp_get_wtime();
		printf("ThreadCount: %d,versionchoose: %d\n", threadCount,version);
		printf("With n = %d trapezoids, our estimate\n", n);
		printf("of the integral from %f to %f = %f\n",a, b, global_result);
		printf("Time used in parallel block = %e seconds\n", finishTime-startTime);
		
	}
	return 0;
}

void Trap(double a, double b, int n, double* global_result_p) {
   double  h, x, my_result;
   double  local_a, local_b;
   int  i, local_n;
   int my_rank = omp_get_thread_num();
   int thread_count = omp_get_num_threads();

   h = (b-a)/n;
   local_n = n/thread_count;
   local_a = a + my_rank*local_n*h;
   local_b = local_a + local_n*h;
   my_result = (f(local_a) + f(local_b))/2.0;
   for (i = 1; i <= local_n-1; i++) {
     x = local_a + i*h;
     my_result += f(x);
   }
   my_result = my_result*h;
   #  pragma omp critical
   *global_result_p += my_result;
}


