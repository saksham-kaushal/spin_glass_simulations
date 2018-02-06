#include<stdio.h>

double ran2_ts(long idum);

int main(void)
{
	int i;
	long iseed=-123456789;
	for(i=0;i<200;i++)
	{
		printf("%f\t",(float)ran2_ts(iseed));
	}
}
