#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

//______________global constants definition_________________

#define N 16
#define NDIM 2
#define NS 256						//N**NDIM
#define N2 NS/2						//NS/2
#define NCONFIG 1
#define NSQ2 NS/2					//NS/2
#define TEMP 8


//___________________ prototypes ___________________________

void distpbc(int n,float distinv[n][n], float dpbc[n][n]);
void rsiteen(int n, float phi[n][n]);
void initialocc(int ns, int n2, float occ2[ns]);
void occupancy(int n, int ns, float occ2[ns], float occ1[n][n]);
void rocc(int n, float occ[3*n][3*n], float occ1[n][n], int nsq2);
void siteenergy(int n, int ns, float occ[3*n][3*n], float phi[n][n], float sten1[n][n]);
void hamiltonian(int n, int ns, float occ[3*n][3*n], float phi[n][n], float *toten, float *totenst);
void choosesite(int n, int ns, float occ1[n][n], int *ie, int *je, int *ih, int *jh);
float delsten(int n, float sten1[n][n], int ie, int je, int ih, int jh, float distinv[n][n]);
void updatehe(int n, int ie, int je, int ih, int jh, float sten1[n][n], float distinv[n][n]);

//___________________ main program _________________________

int main(void)
{

//========== variables
	int iseed, it, nmeas1[TEMP], nskip1[TEMP], i, j, k, ie, je, ih, jh, itemp, occint1[N][N]={0};
	float betaa[TEMP], distinv[N][N]={0}, dpbc[N][N]={0}, phi[N][N]={0}, occ2[NS]={0}, occ1[N][N]={0}, occ[3*N][3*N]={0}, sten1[N][N]={0}, occmin[N][N]={0}, toten=0.0, totenst=0.0;
	
	int ntime, nmeas, nskip, nmcs, imeas, iskip;
	float beta, totmin, delst, a, c;


//========= open files

							//fort.1
	FILE *fp1;
	fp1 = fopen("./fort.1","r");
	if (fp1 == NULL) {fprintf(stderr,"error opening fort.1\n");exit(EXIT_FAILURE);}

							//fort.7
	FILE *fp7;
	fp7 = fopen("./fort.7","r");
	if (fp7 == NULL) {fprintf(stderr,"error opening fort.7\n");exit(EXIT_FAILURE);}

							//fort.2.2
	FILE *fp2;
	fp2 = fopen("./fort.2.2","w");
	if (fp2 == NULL) {fprintf(stderr,"error opening fort.2.2\n");exit(EXIT_FAILURE);}

							//fort.3.2
	FILE *fp3;
	fp3 = fopen("./fort.3.2","w");
	if (fp3 == NULL) {fprintf(stderr,"error opening fort.3.2\n");exit(EXIT_FAILURE);}
	
							//fort.4.2
	FILE *fp4;
	fp4 = fopen("./fort.4.2","w");
	if (fp4 == NULL) {fprintf(stderr,"error opening fort.4.2\n");exit(EXIT_FAILURE);}

							//fort.8.2
	FILE *fp8;
	fp8 = fopen("./fort.8.2","w");
	if (fp8 == NULL) {fprintf(stderr,"error opening fort.8.2\n");exit(EXIT_FAILURE);}

							//fort.10
	FILE *fp10;
	fp10 = fopen("./fort.10","w");
	if (fp10 == NULL) {fprintf(stderr,"error opening fort.10\n");exit(EXIT_FAILURE);}

	
//==============

	iseed =-1234567;
	srand(iseed);
	
	for(it=0;it<TEMP;it++)
	{
		fscanf(fp1,"%f", &betaa[it]);
		fscanf(fp7,"%d\t%d", &nmeas1[it], &nskip1[it]);
	}
/*								//for checking distinv and dpbc values, uncomment
	int jt;
	for(it=0;it<N;it++)
	{	for(jt=0;jt<N;jt++)
	{
		printf("%f %f\n", distinv[it][jt], dpbc[it][jt]);
	}}
*/
	distpbc(N,distinv,dpbc);
/*
	printf("\nviolaaa\n\n");	
	for(it=0;it<N;it++)
	{	for(jt=0;jt<N;jt++)
	{
		printf("%f %f\n", distinv[it][jt], dpbc[it][jt]);
	}}
*/

	for(k=0;k<NCONFIG;k++)
	{
		rsiteen(N,phi);

		for(i=0;i<N;i++)
		{
			for(j=0;j<N;j++)
			{
				fprintf(fp2,"%f\n",phi[i][j]);
			}
		}

		initialocc(NS,N2,occ2);

/*								//for getting occ2
		for(i=0;i<NS;i++)
		{
			printf("%f",occ2[i]);
		}
		printf("\n");
*/
		occupancy(N,NS,occ2,occ1);
		
								//for getting occ1
		for(i=0;i<N;i++)
		{
			for(j=0;j<N;j++)
			{
				printf("%f ",occ1[i][j]);
			}
		}
		printf("\n");

		
		rocc(N,occ,occ1,NSQ2);
		
/*								//for getting occ
		for(i=0;i<3*N;i++)
		{
			for(j=0;j<3*N;j++)
			{
				printf("%f ",occ[i][j]);
			}
		}
		printf("\n");
*/
		for (i=0;i<N;i++)
		{
			for(j=0;j<N;j++)
			{
				fprintf(fp3,"%8.5f\n",occ[i][j]);
			}
		}
		
		siteenergy(N,NS,occ,phi,sten1);
		
/*								//for getting sten1
		for(i=0;i<N;i++)
		{
			for(j=0;j<N;j++)
			{
				printf("%f ",sten1[i][j]);
			}
			printf("\n");
		}
		printf("\n");
*/

		hamiltonian(N,NS,occ,phi,&toten,&totenst);

//		printf("%f\n",toten);	//for getting toten
//		printf("%f\n",totenst);	//for getting totenst


		for(itemp=0;itemp<TEMP;itemp++)
		{
			beta=betaa[itemp];
			clock_t start = clock();
			nmeas = nmeas1[itemp];
			nskip = nskip1[itemp];
			fprintf(fp4,"0\t%8.2f\n",toten);
			ntime = 0;
			totmin = 0.0;
			
//========= Monte-Carlo Simulation
			
			for(imeas=0;imeas<nmeas;imeas++)
			{
				for(iskip=0;iskip<nskip;iskip++)
				{
					for(nmcs=0;nmcs<N2;nmcs++)
					{
						choosesite(N,NS,occ1,&ie,&je,&ih,&jh);
						delst = delsten(N, sten1,ie,je,ih,jh,distinv);
						
						if(delst<0.0)
						{
							occ1[ie][je]=-occ1[ie][je];
							occ1[ih][jh]=-occ1[ih][jh];
							toten+=delst;
							updatehe(N,ie,je,ih,jh,sten1,distinv);
							ntime++;
							if (toten<totmin)
							{
								totmin = toten;
								for(i=0;i<N;i++)
								{
									for(j=0;j<N;j++)
									{
										occmin[i][j] = occ1[i][j];
									}
								}
							}
						}
						else
						{
							a = exp(-delst*beta);
							c = (float)rand()/(float)RAND_MAX;
							if(c<a)
							{
								occ1[ie][je]=-occ1[ie][je];
								occ1[ih][jh]=-occ1[ih][jh];
								toten+=delst;
								updatehe(N,ie,je,ih,jh,sten1,distinv);
								ntime++;
								if(toten<totmin)
								{
									totmin = toten;
									for(i=0;i<N;i++)
									{
										for(j=0;j<N;j++)
										{
											occmin[i][j] = occ1[i][j];
										}
									}
								}
							}
						}
					} //nmcs
				} //iskip
			} //imeas
			fprintf(fp4,"%14.4f\t%14.4f\n", toten, totmin);
			for(i=0;i<N;i++)
			{
				for(j=0;j<N;j++)
				{
					occint1[i][j]=(int)(occmin[i][j]+0.5);
					fprintf(fp8,"%d\n", occint1[i][j]);
				}
			}
			clock_t finish = clock();
			printf("k=%d\ttemp=%f\ttime=%f\trun=%d\ttoten=%f\n",k,beta,((float)(finish-start)/(float)CLOCKS_PER_SEC),ntime,toten);
			
		fprintf(fp10,"%f\n", beta);
		fprintf(fp10,"%f\n", toten);
		for(j=0;j<N;j++)
		{
			for(k=0;k<N;k++)
				{
					fprintf(fp10,"%f", occ[j][k]);
				}
		}
		fprintf(fp10,"\n");
		for(j=0;j<N;j++)
		{
			for(k=0;k<N;k++)
				{
					fscanf(fp10,"%f", &sten1[j][k]);
				}
		}
		fprintf(fp10,"\n");

			
		} //itemp
	} //iconfig

//======= close file descriptors

fclose(fp1);
fclose(fp7);
fclose(fp2);
fclose(fp3);
fclose(fp4);
fclose(fp8);
fclose(fp10);

}

//________________________________________________________

//		distpbc function
//		calculates 1/r term using pbc

void distpbc(int n,float distinv[n][n], float dpbc[n][n])
{
	int ix, iy,k ,j;
	for(ix=0;ix<n;ix++)
	{
		for(iy=0;iy<n;iy++)
		{
			k=ix;
			j=iy;
			if(ix>n/2) {k=abs(ix-n);}
			if(iy>n/2) {j=abs(iy-n);}
			dpbc[ix][iy]=sqrt((float)(k*k)+(float)(j*j));
			if(ix==0.0 && iy==0.0) {distinv[ix][iy]=0.0;}
			else {distinv[ix][iy]=1/dpbc[ix][iy];}
		}
	}

}

//_________________________________________________________


//		rsiteen function			
//		initializes phi matrix to 0.0

void rsiteen(int n, float phi[n][n])
{
	int i, j;
	for (i=0;i<n;i++)
	{
		for (j=0;j<n;j++)
		{
			phi[i][j] = 0.0;
		}
	}

}

//_________________________________________________________

//		initialocc function
//		random initial occ [e:0.5 , h:-0.5]


void initialocc(int ns, int n2, float occ2[ns])
{
	int i, j;
	for(i=0;i<ns;i++)
	{
		occ2[i] = 0.5;
	}

	i=0;
	while(i<n2)
	{
		j = (int)(ns*((float)rand()/(float)RAND_MAX));
		if(occ2[j]!=-0.5)
		{
			occ2[j] = -0.5;
			i++;
		}
	}
}

//__________________________________________________________

//		occupancy function
//		coordinates form

void occupancy(int n, int ns, float occ2[ns], float occ1[n][n])
{
	int i, x1, ix, iy;
	for(i=0;i<ns;i++)
	{
		x1=i%n;
		ix=x1-1;
		if(x1==0)
			{
				ix = n-1;
				iy = i/n-1;
			}
		else
			{
				iy=i/n;
			}
		occ1[ix][iy]=occ2[i];
	}
}

//__________________________________________________________

//		rocc function
//		minimum images
//						---translated indices of occ from -n:2n to 0:3n

void rocc(int n, float occ[3*n][3*n], float occ1[n][n], int nsq2)
{
	int i, j, nx, ny, ix, jx;
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			for(nx=0;nx<3;nx++)
			{
				for(ny=0;ny<3;ny++)
				{
					ix = i+nx*n;
					jx = j+ny*n;
					
					occ[ix][jx]=occ1[i][j];
				}
			}
		}
	}
}

//__________________________________________________________

//		siteenergy function
//		calculation of hartree energy
//						---translated indices of occ from -n:2n to 0:3n

void siteenergy(int n, int ns, float occ[3*n][3*n], float phi[n][n], float sten1[n][n])
{
	int ix, iy, jx, jy, nrc;
	float dist1;
	for(ix=0;ix<N;ix++)
	{
		for(iy=0;iy<N;iy++)
		{
			sten1[ix][iy]=phi[ix][iy];
			nrc = n/2;
			for(jx=ix-nrc+1;jx<ix+nrc+1;jx++)
			{
				for(jy=iy-nrc+1;jy<iy+nrc+1;jy++)
				{
					if(!((jx==ix)&&(jy==iy)))
					{
						dist1=sqrt((float)((ix-jx)*(ix-jx))+(float)((iy-jy)*(iy-jy)));
						sten1[ix][iy]+=((occ[jx+n][jy+n])/dist1);
					}
				}
			}
		}
	}
}

//_________________________________________________________

//		hamiltonian function
//		calculation of Hamiltonian

void hamiltonian(int n, int ns, float occ[3*n][3*n], float phi[n][n], float *toten, float *totenst)
{
	int nrc, i, j, ix, iy, jx, jy;
	float toten1, dist1;
	nrc = n/2;
	*toten = 0.0;
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			*toten+=(phi[i][j]*occ[i+n][j+n]);
		}
	}
	toten1=0.0;
	
	for(ix=0;ix<n;ix++)
	{
		for(iy=0;iy<n;iy++)
		{
			for(jx=ix-nrc+1;jx<ix+nrc+1;jx++)
			{
				for(jy=iy-nrc+1;jy<iy+nrc+1;jy++)
				{
					if (!(ix==jx)&&(iy==jy))
					{
						dist1=sqrt((float)((ix-jx)*(ix-jx))+(float)((iy-jy)*(iy-jy)));
						toten1+=((occ[ix+n][iy+n]*occ[jx+n][jy+n])/dist1);
					}
				}
			}
		}
	}
	*toten+=(toten1/2.0);
	*totenst=*toten/ns;
}

//________________________________________________________

//		choosesite function
//		randomly choose 2 sites for exchange

void choosesite(int n, int ns, float occ1[n][n], int *ie, int *je, int *ih, int *jh)
{
	int kk, x1, x2, ic, kc, n4;
	n4=pow(n,4);
	for(kk=0;kk<n4;kk++)
	{
		ic = (int)(((float)rand()/(float)RAND_MAX)*ns+1);
		x1 = ic%n;
		*ie = x1-1;
		
		if(x1==0)
		{
			*ie = n-1;
			*je = ic/n-1;
		}
		else
		{
			*je = ic/n;
		}
		
		if((occ1[*ie][*je])>0.0)
		{
			break;
		}
	}
	
	for(kk=0;kk<n4;kk++)
	{
		kc = (int)(((float)rand()/(float)RAND_MAX)*ns+1);
		x2 = kc%n;
		*ih = x2-1;
		if(x2==0)
		{
			*ih = n-1;
			*jh = kc/n-1;
		}
		else
		{
			*jh = kc/n;
		}
		
		if(occ1[*ih][*jh]<0.0)
		{
			break;
		}
	}
}

//__________________________________________________________

//		delsten function
//		calculation of change in hartree energy

float delsten(int n, float sten1[n][n], int ie, int je, int ih, int jh, float distinv[n][n])
{
	int ix, iy;
	float delst;
	ix = abs(ie-ih);
	iy = abs(je-jh);
	
	delst = sten1[ih][jh]-sten1[ie][je];
	if(!((ih==ie)&&(jh==je)))
	{
		delst-=distinv[ix][iy];
	}
	return delst;
}

//_________________________________________________________

//		updatehe function
//		update hartree energy after exchange

void updatehe(int n, int ie, int je, int ih, int jh, float sten1[n][n], float distinv[n][n])
{
	int ix, iy, ik, ij, il, im;

	for(ix=0;ix<n;ix++)
	{
		for(iy=0;iy<n;iy++)
		{
			ik=abs(ix-ih);
			ij=abs(iy-jh);
			il=abs(ix-ie);
			im=abs(iy-je);
			
			if((ix==ie)&&(iy==je)) {sten1[ix][iy]+=distinv[ik][ij];}
			if((ix==ih)&&(iy==jh)) {sten1[ix][iy]-=distinv[il][im];}
			if(!(((ix==ie)&&(iy==je))||((ix==ih)&&(iy==jh))))
			{
				sten1[ix][iy]+=(distinv[ik][ij]-distinv[il][im]);
			}
		}
	}
}

///////////////////////////////////////////////////////////

