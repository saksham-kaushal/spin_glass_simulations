#include<stdio.h>
#include<stdlib.h>
#include<omp.h>

#define N 10
#define NT 3

double check(int itemp, double occ[N][N])
{
	int a=0,b=0;
	#pragma omp critical (sectiona)
	{
	printf("thread %d entered check.\n",itemp);
	}
	for(a=0;a<N;a++)
	{
		for(b=0;b<N;b++)
		{
			occ[a][b]*=itemp;
		}
	}
	#pragma omp critical (sectionb)
	{
	printf("thread %d exited check.\n",itemp);
	}
	return occ[0][0];
}

int main(void)
{
	double betaa[]={1.0,2.0,3.0};
	int nmeas[]={300,200,100};
	int nskip[]={150,100,50};
	double rand_arr[N][N];
	int i,j,iseed=1234,iseeds_arr[NT];
	for(i=0;i<N;i++)
	{
		for(j=0;j<N;j++)
		{
			rand_arr[i][j]=(double)rand_r(&iseed)/RAND_MAX;
		}
	}
	for(i=0;i<NT;i++)
	{
		iseeds_arr[i]=rand_r(&iseed)%1000;
	}
	double condition_arr[NT]={0};
	int itemps_arr[NT]={0};
	#pragma omp parallel num_threads(NT) \
	 default(none) \
	 shared(betaa,nmeas,nskip,iseeds_arr,condition_arr,itemps_arr) \
	 firstprivate(rand_arr)
	{
		int k,id=omp_get_thread_num();
		itemps_arr[id]=id;
		int itemp=itemps_arr[id];
		double condition;
		condition=check(itemp,rand_arr);
		condition_arr[itemp]=condition;
		#pragma omp barrier
		#pragma omp master
		{
			printf("all threads completed.\nOld itemps: %d, %d, %d.\nNow exchanging ...\n",itemps_arr[0],itemps_arr[1],itemps_arr[2]);
			for(k=0;k<NT;k++)
			{
				printf("condition_arr[%d]=%lf\n",k,condition_arr[k]);
				if(condition_arr[k]<condition_arr[k+1])
				{
					double temp_condition;
					temp_condition=itemps_arr[k];
					itemps_arr[k]=itemps_arr[k+1];
					itemps_arr[k+1]=temp_condition;
				}
			}
			printf("Exchange complete. New itemps: %d, %d, %d.\n",itemps_arr[0],itemps_arr[1],itemps_arr[2]);
		}
		#pragma omp barrier
		itemp=itemps_arr[id];
		condition=check(itemp,rand_arr);
	}
}
