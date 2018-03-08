#include <stdio.h>
#include<stdlib.h>
#include<omp.h>

void edit (int i[], int *j, int threadid)
{
	i[threadid]+=threadid;
	*j+=threadid;
}

int main(void) {
	int i[]={1,2,3,4,5},j=7,k;
    #pragma omp parallel num_threads(5) shared(i) firstprivate(j)
	{edit(i,&j,omp_get_thread_num());
	printf("in thread %d j=%d\n",omp_get_thread_num(),j);
	}
	for (k=0;k<5;k++)
	{
		printf("%d",i[k]);
	}
	printf("outside threads j=%d",j);
	return 0;
}

