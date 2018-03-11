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
#define TEMP 21	//>3
#define NUM_EXCHANGES 1

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

//double ran2(long *idum);
void distpbc(double distinv[N][N],double dpbc[N][N]);
void initialocc(double occ2[NS], int iseed);
void occupancy(double occ2[NS],double occ1[N][N]);
void rocc(double occ[3*N][3*N],double occ1[N][N]);
void siteenergy(double occ[3*N][3*N], double phi[N][N], double sten1[N][N]);
void hamiltonian(double occ[3*N][3*N],double phi[N][N], double* toten, double *totenst);
sites_struct choosesite(double occ1[N][N], int iseed);
double delsten(double sten1[N][N], int ie, int je, int ih, int jh, double distinv[N][N]);
void updatehe(int ie, int je, int ih, int jh, double sten1[N][N], double distinv[N][N]);
double montecarlo(int iseeds_arr[TEMP],double occ1[N][N],double sten1[N][N],double occmin[N][N],sites_struct sites,double delst,double betaa[TEMP],int nmeas1[TEMP],int nskip1[TEMP],double distinv[N][N],int itemp, double toten_arr[TEMP], double totmin_arr[TEMP],double toten);

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
	int iseed;
	int iseeds_arr[TEMP];
	
	int nmeas1[TEMP], nskip1[TEMP];
	
	double betaa[TEMP],distinv[N][N]={0},dpbc[N][N]={0}, phi[N][N]={0}, occ2[NS]={0}, occ1[N][N]={0}, occ[3*N][3*N]={0}, sten1[N][N]={0}, toten=0.0, totenst=0.0, toten_arr[TEMP]={0.0},totmin_arr[TEMP]={0};

//------------------------------------------

	for(it=0;it<TEMP;it++)
	{
		fscanf(fp1,"%lf",&betaa[it]);
		fscanf(fp7,"%d\t%d",&nmeas1[it], &nskip1[it]);
	}
	
/*	for(it=0;it<TEMP;it++)
	{
		printf("%lf\n",betaa[it]);
		printf("%d\t%d\n",nmeas1[it],nskip1[it]);
	}
*/
	
	distpbc(distinv,dpbc);
	
/*	for(it=0;it<N;it++)
	{
		for(jt=0;jt<N;jt++)
		{
			printf("%lf  ",distinv[it][jt]);
		}
		printf("\n");
	}
	printf("\n");
	for(it=0;it<N;it++)
	{
		for(jt=0;jt<N;jt++)
		{
			printf("%lf  ",dpbc[it][jt]);
		}
		printf("\n");
	}
*/
	iseed = 12345;
	
	for (k=0;k<NCONFIG;k++)
	{
		//evaluate and write rsiteen to fp2
		initialocc(occ2, iseed);
/*								//for getting occ2
		for(it=0;it<NS;it++)
		{
			printf("%lf\t",occ2[it]);
		}
		printf("\n\n");
*/
		int x;
		for (x=0;x<TEMP;x++)
		{
			iseeds_arr[x]=rand_r(&iseed);
		}
		occupancy(occ2,occ1);
/*		for(it=0;it<N;it++)
		{
			for(jt=0;jt<N;jt++)
			{
				printf("%lf  ",occ1[it][jt]);
			}
			printf("\n");
		}
*/		
		
		rocc(occ,occ1);
/*		for(it=0;it<3*N;it++)
		{
			for(jt=0;jt<3*N;jt++)
			{
				printf("%lf  ",occ[it][jt]);
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
				printf("%lf  ",sten1[it][jt]);
			}
			printf("\n");
		}
*/
//		printf("%lf\t%lf\n", toten, totenst);

		hamiltonian(occ,phi,&toten,&totenst);
		
//		printf("%lf\t%lf\n", toten, totenst);

		
		double totmin=0.0, occmin[N][N]={0}, delst=0.0;
		int nmeas=0, nskip=0, nmcs=0, i=0;
		sites_struct sites;
		
		int exchange;
		
		#pragma omp parallel num_threads(TEMP) \
		 default(none) \
		 firstprivate(iseeds_arr,toten,totmin,nmeas,nskip,occ1,sten1,occmin,sites,delst, \
		 exchange) \
		 shared(betaa,nmeas1,nskip1,distinv,toten_arr,totmin_arr,iseed) \
 
		{
			int itemp=omp_get_thread_num();	
			toten=montecarlo(iseeds_arr, occ1, sten1, occmin, sites, delst, betaa, nmeas1, nskip1, distinv, itemp, toten_arr, totmin_arr, toten);
			
			for (exchange=0;exchange<NUM_EXCHANGES;exchange++)
			{
				int i;
				#pragma omp barrier
				#pragma omp master
				{
					printf("\n");
					for (i=0;i<TEMP;i++)
					{
						printf("%14.10lf\t",toten_arr[i]);
					}
					printf("\n\nNow exchanging..\n\n");
					int x=0;
					for(x=TEMP-1;x>0;x--)
					{
						double del_e=toten_arr[x-1]-toten_arr[x];
						if(del_e<0)
						{
							double temp_temp=toten_arr[x];
							toten_arr[x]=toten_arr[x-1];
							toten_arr[x-1]=temp_temp;
						}
						else
						{
							if(((double)rand_r(&iseed)/RAND_MAX)<exp(del_e*(betaa[x-1]-betaa[x])))
							{
								double temp_temp=toten_arr[x];
								toten_arr[x]=toten_arr[x-1];
								toten_arr[x-1]=temp_temp;
							}
						}
					}
					for (i=0;i<TEMP;i++)
					{
						printf("%14.10lf\t",toten_arr[i]);
					}
					printf("\n\nExchanged..\n\n");
				}				
				#pragma omp barrier
				
				i=0;
				while(i<=TEMP)
				{
					if (toten_arr[i]==toten)
					{
						itemp=i;
						break;
					}
					i++;
				}
				if(i==TEMP)
				{
					printf("\nNot found\n");
					exit(1);
				}
				printf("itemp=%4d\t\ttoten=%14.10lf\n",itemp,toten);
				
				#pragma omp barrier
				#pragma omp master
				{
					printf("\n");
				}
				
				toten=montecarlo(iseeds_arr, occ1, sten1, occmin, sites, delst, betaa, nmeas1, nskip1, distinv, itemp, toten_arr, totmin_arr, toten);
				#pragma omp barrier
				#pragma omp master
				{
					printf("\n");
					for (i=0;i<TEMP;i++)
					{
						printf("%14.10lf\t",toten_arr[i]);
					}
					printf("\n");
				}
			}	//tempering: exchange block
		}	//tempering: parallel block
	}	//config loop end

//-----------------------------------------
	
	fclose(fp1);
	fclose(fp7);
	fclose(fp3);

}			//main end


//=========================================
//=========================================

double montecarlo(int iseeds_arr[TEMP],double occ1[N][N],double sten1[N][N],double occmin[N][N], sites_struct sites,double delst,double betaa[TEMP],int nmeas1[TEMP],int nskip1[TEMP],double distinv[N][N],int itemp, double toten_arr[TEMP], double totmin_arr[TEMP], double toten)
{
//	printf("%d",omp_get_num_threads());
	double beta=betaa[itemp], totmin=0.0;
	int nmeas=nmeas1[itemp];
	int nskip=nskip1[itemp];
	
	double a=0.0, c=0.0;
	int i=0, j=0;
			
//	printf("%d\n",iseeds_arr[itemp]);
//	double delst, a, c;
//	int imeas, iskip, nmcs, i, j;
//	printf("%ld %lf %lf %d %lf %d %d\n",iseed,toten,totmin, itemp, beta, nmeas, nskip);
	int imeas=0, iskip=0, nmcs=0;
	for(imeas=0;imeas<nmeas;imeas++)
	{
		for(iskip=0;iskip<nskip;iskip++)
		{
			for(nmcs=0;nmcs<N2;nmcs++)
			{
				{
					sites=choosesite(occ1,iseeds_arr[itemp]);
					delst=delsten(sten1, sites.ie, sites.je, sites.ih, sites.jh,distinv);
					if(delst<=0.0)
					{
						occ1[sites.ie][sites.je] = -occ1[sites.ie][sites.je];
						occ1[sites.ih][sites.jh] = -occ1[sites.ih][sites.jh];
						toten+=delst;
						updatehe(sites.ie, sites.je, sites.ih, sites.jh,sten1,distinv);
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
						//printf("%lf  ", beta);
						a = exp(-delst*beta);
//						#pragma omp critical (random_no)
						c = (double)rand_r(&iseeds_arr[itemp])/RAND_MAX;
						if(c<a)
						{
							occ1[sites.ie][sites.je] = -occ1[sites.ie][sites.je];
							occ1[sites.ih][sites.jh] = -occ1[sites.ih][sites.jh];
							toten+=delst;
							updatehe(sites.ie, sites.je, sites.ih, sites.jh,sten1,distinv);
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
				}
			}	// for nmcs loop
		}	//for nskip loop
	}	//omp for nmeas block
	toten_arr[itemp]=toten;
	totmin_arr[itemp]=totmin;
	printf("temp=%10.5lf\t\ttoten=%14.10lf\n",beta,toten);
	return toten;
}	//omp parallel temp block


//+++++++++++++++++++++++++++++++++++++++++

void distpbc(double distinv[N][N],double dpbc[N][N])
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
			
			dpbc[ix][iy]=sqrt((double)((k*k)+(j*j)));
			distinv[ix][iy]=1/dpbc[ix][iy];
		}
	}
	distinv[0][0]=0.0;
}

//++++++++++++++++++++++++++++++++++++++++

void initialocc(double occ2[NS], int iseed)
{
	int i,j;
	double k;
	int seed;
	seed = abs(iseed);
	
	for(i=0;i<NS;i++)
	{
		occ2[i]=0.5;
	}
	i=0;
	while(i<N2)
	{
//		#pragma omp critical (random_no)
		k=(double)rand_r(&seed)/RAND_MAX;
		j=(int)(k*NS);
		if(occ2[j]!=-0.5)
		{
			occ2[j] = -0.5;
			i++;
		}
	}
}

//+++++++++++++++++++++++++++++++++++++++++

void occupancy(double occ2[NS],double occ1[N][N])
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

void rocc(double occ[3*N][3*N],double occ1[N][N])
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

void siteenergy(double occ[3*N][3*N], double phi[N][N], double sten1[N][N])
{
	int ix, iy, jx, jy, nrc;
	double dist1;
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
						dist1=sqrt((double)((-jx+nrc-1)*(-jx+nrc-1))+(double)((-jy+nrc-1)*(-jy+nrc-1)));
						sten1[ix][iy]+=((occ[jx+ix-nrc+1][jy+iy-nrc+1])/dist1);
					}
				}
			}
		}
	}
}

//+++++++++++++++++++++++++++++++++++++++++

void hamiltonian(double occ[3*N][3*N],double phi[N][N], double* toten, double *totenst)
{
	int nrc, i, j, ix, iy, jx, jy;
	double toten1=0.0, dist1;
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
						dist1=sqrt((double)((-jx+nrc-1)*(-jx+nrc-1))+(double)((-jy+nrc-1)*(-jy+nrc-1)));
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

sites_struct choosesite(double occ1[N][N], int iseed)
{
	int kk, x1, x2, ic, kc, n4, ie, je, ih, jh;
	n4=N*N*N*N;		//faster than pow
	int seed = abs(iseed);
	for(kk=0;kk<n4;kk++)
	{
//		#pragma omp critical (random_no)
		
		ic = (int)(((double)rand_r(&seed)/RAND_MAX)*NS+1);
		
		x1 = ic%N;
		ie = x1-1;
		
		if(x1==0)
		{
			ie = N-1;
			je = ic/N-1;
		}
		else
		{
			je = ic/N;
		}
		
		if((occ1[ie][je])>0.0)
		{
			break;
		}
		else {;}
	}
	
	for(kk=0;kk<n4;kk++)
	{
//		#pragma omp critical (random_no)
		
		kc = (int)(((double)rand_r(&seed)/RAND_MAX)*NS+1);
		
		x2 = kc%N;
		ih = x2-1;
		if(x2==0)
		{
			ih = N-1;
			jh = kc/N-1;
		}
		else
		{
			jh = kc/N;
		}
		
		if(occ1[ih][jh]<0.0)
		{
			break;
		}
		else {;}
	}
	sites_struct sites = {ie,je,ih,jh};
	return sites;
}

//++++++++++++++++++++++++++++++++++++++++++++

double delsten(double sten1[N][N], int ie, int je, int ih, int jh, double distinv[N][N])
{
	int ix, iy;
	double delst;
	ix = abs(ie-ih);
	iy = abs(je-jh);
	
	delst = sten1[ih][jh]-sten1[ie][je];
	if(!((ih==ie)&&(jh==je)))
	{
		delst-=distinv[ix][iy];
	}
	return delst;
}

//+++++++++++++++++++++++++++++++++++++++++++

void updatehe(int ie, int je, int ih, int jh, double sten1[N][N], double distinv[N][N])
{
	int ix, iy, ik, ij, il, im;

	for(ix=0;ix<N;ix++)
	{
		for(iy=0;iy<N;iy++)
		{
			ik=abs(ix-ih);
			ij=abs(iy-jh);
			il=abs(ix-ie);
			im=abs(iy-je);
			
			if((ix==ie)&&(iy==je)) 
			
				{sten1[ix][iy]+=distinv[ik][ij];}
				
			if((ix==ih)&&(iy==jh)) 
			
				{sten1[ix][iy]-=distinv[il][im];}
				
			if(!(((ix==ie)&&(iy==je))||((ix==ih)&&(iy==jh))))
			
				{sten1[ix][iy]+=(distinv[ik][ij]-distinv[il][im]);}
		}
	}
}







