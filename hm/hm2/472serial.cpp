#include<stdio.h>


double f(double x){
		double return_val;
		return_val = x * x;
		return return_val;
	}
		
int main(void){
	int  n, i;
	double a, b, h, x;
	double  total_int;

	
	printf("enter a, b, n:");
	
	scanf("%lf %lf %d", &a, &b, &n);
	

	
	h = (b-a)/(double)n;
	total_int = f(a)+f(b);
	x = a;
	for(i=1;i < n; i++){
		x += h;
		if(i%2==1){
			total_int += (double)4*f(x);
		}else{
			total_int += (double)2*f(x);
		}
	}
	
    total_int *= h/3;
	printf("with n = %d trap, our estimate\n", n);
	printf("of the integral from %f to %f = %f\n", a, b, total_int);
	return 0;
} 