#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "my_rand.h"


void Get_args(char* argv[], int* thread_count_p, 
      long long int* number_of_tosses_p);
void Usage(char* prog_name);

long long int Count_hits(long long int number_of_tosses, int thread_count);


int main(int argc, char* argv[]) {
   double pi_estimate;
   int thread_count;
   long long int number_in_circle;
   long long int number_of_tosses;
   
   if (argc != 3) Usage(argv[0]);
   Get_args(argv, &thread_count, &number_of_tosses);
   
   number_in_circle = Count_hits(number_of_tosses, thread_count);

   pi_estimate = 4*number_in_circle/((double) number_of_tosses);
   printf("Estimated pi: %e\n", pi_estimate);

   return 0;
}

double my_drand(unsigned* seed_p);

long long int Count_hits(long long int number_of_tosses, int thread_count) {

   long long int number_in_circle = 0;
   
#  pragma omp parallel num_threads(thread_count) \
      default(none) reduction(+: number_in_circle) \
      shared(number_of_tosses, thread_count)
   {
      int my_rank = omp_get_thread_num();
      unsigned seed = my_rank + 1;
      long long int toss;
      double x, y, distance_squared;

#     pragma omp for
      for(toss = 0; toss < number_of_tosses; toss++) {
         x = 2*my_drand(&seed) - 1;
         y = 2*my_drand(&seed) - 1;
         distance_squared = x*x + y*y;
         if (distance_squared <= 1) number_in_circle++;
#        ifdef DEBUG
         printf("Th %d > toss = %lld, x = %.3f, y = %.3f, dist = %.3f\n",
            my_rank, toss, x, y, distance_squared);
#        endif
      }
   }  /* pragma omp parallel */

#  ifdef DEBUG
   printf("Total number in circle = %lld\n", number_in_circle);
#  endif
   
   return number_in_circle;
}  

void Usage(char prog_name[] /* in */) {
   fprintf(stderr, "usage: %s ", prog_name); 
   fprintf(stderr, "<number of threads> <total number of tosses>\n");
   exit(0);
}  

void Get_args(
           char*           argv[]              /* in  */,
           int*            thread_count_p      /* out */,
           long long int*  number_of_tosses_p  /* out */) {
   
   *thread_count_p = strtol(argv[1], NULL, 10);  
   *number_of_tosses_p = strtoll(argv[2], NULL, 10);
} 

