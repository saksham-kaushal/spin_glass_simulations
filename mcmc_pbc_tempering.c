#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<omp.h>

#define N 16
#define NDIM 2
#define NS 256
#define N2 NS/2
#define NCONFIG 1
#define NSQ2 NS/2
#define TEMP 10

//========================================
//========================================

typedef struct
{
	int ie;
	int je;
	int ih;
	int jh;
} sites_struct;

//========================================
//========================================

double ran2(long *idum);
void distpbc(float distinv[N][N],float dpbc[N][N]);
void initialocc(float occ2[NS], long iseed);
void occupancy(float occ2[NS],float occ1[N][N]);
void rocc(float occ[3*N][3*N],float occ1[N][N]);
void siteenergy(float occ[3*N][3*N], float phi[N][N], float sten1[N][N]);
void hamiltonian(float occ[3*N][3*N],float phi[N][N], float* toten, float *totenst);
sites_struct choosesite(float occ1[N][N], long iseed);
float delsten(float sten1[N][N], int ie, int je, int ih, int jh, float distinv[N][N]);
void updatehe(int ie, int je, int ih, int jh, float sten1[N][N], float distinv[N][N]);

//=========================================
//=========================================

int main(void)
{
	FILE *fp1;
	fp1 = fopen("./fort.1","r");
	if (fp1==NULL) {fprintf(stderr,"error opening fort.1\n"); exit(EXIT_FAILURE);}
	
	FILE *fp7;
	fp7 = fopen("./fort.7","r");
	if (fp7 == NULL) {fprintf(stderr,"error opening fort.7\n");exit(EXIT_FAILURE);}
	
	FILE *fp3;
	fp3 = fopen("./fort.3","w");
	if (fp3 == NULL) {fprintf(stderr,"error opening fort.3\n");exit(EXIT_FAILURE);}
	
//------------------------------------------
	
	int it, jt, k;
	long iseed;
	
	int nmeas1[TEMP], nskip1[TEMP];
	
	float betaa[TEMP],distinv[N][N]={0},dpbc[N][N]={0}, phi[N][N]={0}, occ2[NS]={0}, occ1[N][N]={0}, occ[3*N][3*N]={0}, sten1[N][N]={0}, toten=0.0, totenst=0.0;

//------------------------------------------

	for(it=0;it<TEMP;it++)
	{
		fscanf(fp1,"%f",&betaa[it]);
		fscanf(fp7,"%d\t%d",&nmeas1[it], &nskip1[it]);
	}
	
/*	for(it=0;it<TEMP;it++)
	{
		printf("%f\n",betaa[it]);
		printf("%d\t%d\n",nmeas1[it],nskip1[it]);
	}
*/
	
	distpbc(distinv,dpbc);
	
/*	for(it=0;it<N;it++)
	{
		for(jt=0;jt<N;jt++)
		{
			printf("%f  ",distinv[it][jt]);
		}
		printf("\n");
	}
	printf("\n");
	for(it=0;it<N;it++)
	{
		for(jt=0;jt<N;jt++)
		{
			printf("%f  ",dpbc[it][jt]);
		}
		printf("\n");
	}
*/
	iseed = -1234567;
	for (k=0;k<NCONFIG;k++)
	{
		//evaluate and write rsiteen to fp2
		initialocc(occ2, iseed);
/*								//for getting occ2
		for(it=0;it<NS;it++)
		{
			printf("%f\t",occ2[it]);
		}
		printf("\n\n");
*/
		
		occupancy(occ2,occ1);
/*		for(it=0;it<N;it++)
		{
			for(jt=0;jt<N;jt++)
			{
				printf("%f  ",occ1[it][jt]);
			}
			printf("\n");
		}
*/		
		
		rocc(occ,occ1);
/*		for(it=0;it<3*N;it++)
		{
			for(jt=0;jt<3*N;jt++)
			{
				printf("%f  ",occ[it][jt]);
			}
			printf("\n");
		}
*/
		for (it=0;it<N;it++)
		{
			for(jt=0;jt<N;jt++)
			{
				fprintf(fp3,"%8.5f\n",occ[it][jt]);
			}
		}
		
		siteenergy(occ,phi,sten1);
		
/*		for(it=0;it<N;it++)
		{
			for(jt=0;jt<N;jt++)
			{
				printf("%f  ",sten1[it][jt]);
			}
			printf("\n");
		}
*/
		printf("%f\t%f\n", toten, totenst);

		hamiltonian(occ,phi,&toten,&totenst);
		
		printf("%f\t%f\n", toten, totenst);

		
		float totmin=0.0, beta=0.0, occmin[N][N]={0}, delst=0, a=0, c=0;
		int nmeas=0, nskip=0, itemp=0, imeas=0, iskip=0, nmcs=0, i=0, j=0;
		sites_struct sites;
		
		
		#pragma omp parallel for \
		 firstprivate(iseed,toten,totmin,nmeas,nskip,occ1,sten1,occmin,itemp,imeas,iskip,nmcs,sites,i,j,delst,a,c,beta) \
		 shared(betaa,nmeas1,nskip1,distinv)

		{
			for(itemp=0;itemp<TEMP;itemp++)
			{
				
			}
		}	//omp parallel for block
	}	//config loop end

//-----------------------------------------
	
	fclose(fp1);
	fclose(fp7);
	fclose(fp3);

}			//main end


//=========================================
//=========================================

void distpbc(float distinv[N][N],float dpbc[N][N])
{
	int ix, iy, j, k;
	for(ix=0;ix<N;ix++)
	{
		for(iy=0;iy<N;iy++)
		{
			k=ix;
			j=iy;
			if(ix>N/2) {k=abs(ix-N);}
			if(iy>N/2) {j=abs(iy-N);}
			
			dpbc[ix][iy]=sqrt((float)((k*k)+(j*j)));
			distinv[ix][iy]=1/dpbc[ix][iy];
		}
	}
	distinv[0][0]=0.0;
}

//++++++++++++++++++++++++++++++++++++++++

void initialocc(float occ2[NS], long iseed)
{
	int i,j;
	double k;
	long seed;
	seed = -abs(iseed);
	
	for(i=0;i<NS;i++)
	{
		occ2[i]=0.5;
	}
	i=0;
	while(i<N2)
	{
		#pragma omp critical (random_no)
		k=ran2(&seed);
		j=(int)(k*NS);
		if(occ2[j]!=-0.5)
		{
			occ2[j] = -0.5;
			i++;
		}
	}
}

//+++++++++++++++++++++++++++++++++++++++++

void occupancy(float occ2[NS],float occ1[N][N])
{
	int i, x1, ix, iy;
	for(i=0;i<NS;i++)
	{
		x1=(i+1)%N;
		ix=x1-1;
		if(x1==0)
			{
				ix = N-1;
				iy = (i+1)/N-1;
			}
		else
			{
				iy=(i+1)/N;
				
			}
		occ1[ix][iy]=occ2[i];
	}
}

//+++++++++++++++++++++++++++++++++++++++

void rocc(float occ[3*N][3*N],float occ1[N][N])
{
	int i, j, nx, ny, ix, jx;
	for(i=0;i<N;i++)
	{
		for(j=0;j<N;j++)
		{
			for(nx=0;nx<3;nx++)
			{
				for(ny=0;ny<3;ny++)
				{
					ix = i+nx*N;
					jx = j+ny*N;
					
					occ[ix][jx]=occ1[i][j];
				}
			}
		}
	}
}

//+++++++++++++++++++++++++++++++++++++++++

void siteenergy(float occ[3*N][3*N], float phi[N][N], float sten1[N][N])
{
	int ix, iy, jx, jy, nrc;
	float dist1;
	for(ix=0;ix<N;ix++)
	{
		for(iy=0;iy<N;iy++)
		{
			sten1[ix][iy]=phi[ix][iy];
			nrc = N/2;
			for(jx=0;jx<ix+ix+2;jx++)
			{
				for(jy=0;jy<iy+iy+2;jy++)
				{
					if(!((jx==nrc-1)&&(jy==nrc-1)))
					{
						dist1=sqrt((float)((-jx+nrc-1)*(-jx+nrc-1))+(float)((-jy+nrc-1)*(-jy+nrc-1)));
						sten1[ix][iy]+=((occ[jx+ix-nrc+1][jy+iy-nrc+1])/dist1);
					}
				}
			}
		}
	}
}

//+++++++++++++++++++++++++++++++++++++++++

void hamiltonian(float occ[3*N][3*N],float phi[N][N], float* toten, float *totenst)
{
	int nrc, i, j, ix, iy, jx, jy;
	float toten1=0.0, dist1;
	nrc = N/2;
	
	for(i=0;i<N;i++)
	{
		for(j=0;j<N;j++)
		{
			*toten+=phi[i][j]*occ[i][j];
		}
	}
	for(ix=0;ix<N;ix++)
	{
		for(iy=0;iy<N;iy++)
		{
			for(jx=0;jx<ix+ix+2;jx++)
			{
				for(jy=0;jy<iy+iy+2;jy++)
				{
					if(!((jx==nrc-1)&&(jy==nrc-1)))
					{
						dist1=sqrt((float)((-jx+nrc-1)*(-jx+nrc-1))+(float)((-jy+nrc-1)*(-jy+nrc-1)));
						toten1+=((occ[ix][iy]*occ[jx+ix-nrc+1][jy+iy-nrc-1])/dist1);
					}
				}
			}
		}
	}
	*toten+=(toten1/2.0);
	*totenst=*toten/NS;
}

//+++++++++++++++++++++++++++++++++++++++++

