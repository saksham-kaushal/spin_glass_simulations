#include<stdio.h>
#include<stdlib.h>
#include<omp.h>

int main(void)
{
	int i;
	int iseed[]={100,200,300,400,500};
	float rans[10];
	#pragma omp parallel num_threads(5) \
	 firstprivate(iseed,i,rans)

	{
	for(i=0;i<10;i++) printf("%f  ",((float)rand_r(&iseed[omp_get_thread_num()]))/(RAND_MAX));
	}
}
