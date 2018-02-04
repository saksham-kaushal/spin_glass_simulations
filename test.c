#include <stdio.h>
#include<stdlib.h>

int main(void) {
    int i;
    for(i=0;i<5;i++)
    {
    	int seed = 5;
    	srand(seed);
    	printf("%f",(float)rand());
    	printf("%f",(float)rand());
    	printf("%f",(float)rand());
    	printf("%f",(float)rand());
    	printf("\n");
    }
	return 0;
}

