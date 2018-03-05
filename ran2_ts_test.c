#include<stdio.h>

double ran2(long seed);

int main(void)
{
	long iseed=-123456789;
	int i;
	for(i=0;i<10;i++) printf("%f",(float) ran2(iseed));
}
