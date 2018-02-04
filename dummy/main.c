#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define TEMP 3
#define N 16

int main(void)
{
	int i, j, k;
	FILE *fp1, *fp2;

	fp1 = fopen("./fort.10","r");
	if (fp10 == NULL) {fprintf(stderr,"error opening fort.10\n");exit(EXIT_FAILURE);}

	fp2 = fopen("./fort.10","r");
	if (fp17 == NULL) {fprintf(stderr,"error opening fort.10\n");exit(EXIT_FAILURE);}
	
	float beta[TEMP], totens[TEMP], occs[TEMP][N][N], stens[TEMP][N][N];
	for(i=0;i<TEMP;i++)
	{
		fscanf(fp1,"%f", &betas[i]);
		fscanf(fp1,"%f", &totens[i]);
		for(j=0;j<N;j++)
		{
			for(k=0;k<N;k++)
				{
					fscanf(fp1,"%f", &occs[i][j][k]);
				}
		}
		for(j=0;j<N;j++)
		{
			for(k=0;k<N;k++)
				{
					fscanf(fp1,"%f", &stens[i][j][k]);
				}
		}
		
	}
	
	for(i=0;i<TEMP-1;i++)
	{
		float b, a;
		float occtemp[N][N];
		if((toten[i]-toten[i+1])!=0)
		{
			b = exp(-(betas[i+1]-betas[i])*(toten[i]-toten[i+1]));
			a = ((float)rand()/(float)RAND_MAX);
		}
		if((toten[i]<toten[i+1]) || (b<a))
		{
			memcpy(occtemp,occs[i+1],sizeof(occtemp));
			for(j=0;j<N;j++)
			{
				for(k=0;k<N;k++)
				{
					occs[i+1][j][k]=occs[i][j][k];
					occs[i][j][k]=occtemp[i][j];
				}
			}
		}
		else
		{
			continue;
		}
	}
	
	for(i=0;i<TEMP;i++)
	{
		fprintf(fp1,"%f\n", &betas[i]);
		fprintf(fp1,"%f\n", &totens[i]);
		for(j=0;j<N;j++)
		{
			for(k=0;k<N;k++)
				{
					fprintf(fp1,"%f", &occs[i][j][k]);
				}
		}
		fprintf(fp1,"\n");
		for(j=0;j<N;j++)
		{
			for(k=0;k<N;k++)
				{
					fscanf(fp1,"%f", &stens[i][j][k]);
				}
		}
		fprintf(fp1,"\n");
	}

}
