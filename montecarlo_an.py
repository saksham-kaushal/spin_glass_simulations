import numpy as np
from time import process_time
import pdb

class const(object):
	n=16
	ndim=2
	ns=n**ndim
	n2=ns/2
	nconfig=1
	nsq2=ns/2
	temp=21

	def __setattr__(self, *_): 
		raise AttributeError("Assignment to a constant is not permitted.")

##########################################################

def distpbc(n,distinv):
	dpbc=np.zeros((n,n))
	for ix in range(n):
		for iy in range(n):
			k=ix
			j=iy
			if (ix > n/2):
				k=abs(ix-n)
			if (iy > n/2):
				j=abs(iy-n)
			dpbc[ix][iy]=np.sqrt(float(k*k)+float(j*j))
			if ((ix == 0) and (iy == 0)):
				distinv[ix][iy] = float(0)
			else :
				distinv[ix][iy] = 1/dpbc[ix][iy]
	return distinv
			
##########################################################

"""
def rsiteen(n,iseed,phi):
	for i in range(n):
		for j in range(n):
			phi[i][j]=0
	return phi
"""

#phi matrix already initialized to zero. Consider return in main, if value returned by function changes from zero to some other value.

#########################################################


def initialocc(ns,n2,iseed,occ2):
	occ2=np.full((ns),0.5)
	i=0
	while (i<n2):
		j=int(ns*np.random.rand())
		if (occ2[j] == -0.5):
			continue
		occ2[j] = -0.5
		i+=1
	return occ2
	
#########################################################

def occupancy(n,ns,occ2,occ1):
	for i in range(1,ns+1):
		x1=i%n
		ix = x1-1
		if (x1==0):
			ix = n-1
			iy = (i/n)-1
		else:
			iy=i/n
		occ1[int(ix)][int(iy)] = occ2[i-1]
	return occ1
	
########################################################

def rocc(n,iseed,occ,occ1,nsq2):
	for i in range(n):
		for j in range(n):
			for nx in range(-1,2):
				for ny in range(-1,2):
					ix=i+nx*n
					jx=j+ny*n
					occ[ix][jx]=occ1[i][j]
	return occ
	
########################################################

def siteenergy(n,ns,occ,phi,sten1):
	for ix in range(n):
		for iy in range(n):
			sten1[ix][iy]=phi[ix][iy]
			nrc=n/2
			for jx in range(int(ix-nrc+1),int(ix+nrc+1)):
				for jy in range(int(iy-nrc+1),int(iy+nrc+1)):
					if((ix==jx) and (iy==jy)):
						sssssss=1.0
					else:
						dist1=np.sqrt(float((ix-jx)**2)+float((iy-jy)**2))
						sten1[ix][iy]+=(occ[jx][jy]/dist1)
	return sten1
						
########################################################

def hamiltonian(n,ns,occ,phi):
	nrc=n/2
	toten=0.0
	for i in range(n):
		for j in range(n):
			toten+=phi[i][j]*occ[i][j]
			
	toten1=0.0
	for ix in range(n):
		for iy in range(n):
			for jx in range(int(ix-nrc+1),int(ix+nrc+1)):
				for jy in range(int(iy-nrc+1),int(iy+nrc+1)):
					if ((ix==jx) and (iy==jy)):
						sssssss=1.0
					else:
						dist1=np.sqrt(float((ix - jx)**2) + float((iy - jy)**2))
						toten1+=((occ[ix][iy]*occ[jx][jy])/dist1)
	toten+=toten1/2.0
	totenst=toten/ns
	return toten,totenst
	
########################################################

def choosesite(n,ns,iseed,occ1):
	
	for kk in range(n**4):
		ic=int(ns*np.random.rand()+1)
		x1=ic%n
		ie=int(x1-1)
		if(x1==0):
			ie=int(n-1)
			je=int(ic/n-1)
		else:
			je=int(ic/n)
		if(occ1[ie][je]>0.0):
			ss=1.0
			break
			
	for kk in range(n**4):
		kc=int(ns*np.random.rand()+1)
		x2=kc%n
		ih=int(x2-1)
		if(x2==0):
			ih=int(n-1)
			jh=int(kc/n-1)
		else:
			jh=int(kc/n)
		if(occ1[ih][jh]<0.0):
			sss=1.0
			break
			
	return ie,je,ih,jh
	
#########################################################

def delsten(sten1,ie,je,n,ih,jh,distinv):
	ix=abs(ie-ih)
	iy=abs(je-jh)
	delst=sten1[ih][jh]-sten1[ie][je]
	if(ih==ie and jh==je):
		ssss=1.0
	else:
		delst-=distinv[ix][iy]
		
	return delst
	
#########################################################

def updatehe(n,ie,je,ih,jh,sten1,distinv):
	nrc=n/2
	for ix in range(n):
		for iy in range(n):
			ik=abs(ix-ih)
			ij=abs(iy-jh)
			il=abs(ix-ie)
			im=abs(iy-je)
			
			if(ix==ie and iy==je):
				sten1[ix][iy]+=distinv[ik][ij]
			elif(ix==ih and iy==jh):
				sten1[ix][iy]-=distinv[il][im]
			else:
				sten1[ix][iy]=sten1[ix][iy]+(distinv[ik][ij]-distinv[il][im])
			sssssss=1.0
				
	return sten1


#######################################################

#def updateocc not addded.







##########             main               ###############

const=const() #declaration of constants is in const class


iseed = -1234567		#unused for random number generator
with open('fort.1','r') as f1:
	betaa=np.fromfile(f1, dtype=float, count=const.temp, sep='\n\r')

print (betaa)

with open('fort.7','r') as f7:
	nm_ns=np.fromfile(f7, dtype=float, count=2*const.temp, sep='\n\r')

print (nm_ns)

nmeas1=nm_ns[0::2]
nskip1=nm_ns[1::2]

print (nmeas1)
print (nskip1)

phi=np.zeros((const.n,const.n))

distinv=np.zeros((const.n,const.n))
distinv=distpbc(const.n,distinv)

print (distinv)

occ2=np.zeros((const.ns))
occ1=np.zeros((const.n,const.n))
sten1=np.zeros((const.n,const.n))
occmin=np.zeros((const.n,const.n))
occint1=np.zeros((const.n,const.n))

with open('fort.21','wb') as f2, open('fort.31','wb') as f3, open('fort.41','w') as f4,open('fort.81','w') as f8:
	for k in range(const.nconfig):
		#rsiteen(const.n,iseed,phi)
		float_formatter = lambda x: "%08.5f" % x
		np.set_printoptions(formatter={'float_kind':float_formatter})
		np.savetxt(f2,phi,delimiter='\n',fmt='%f')
		np.set_printoptions()
		
		occ2=initialocc(const.ns,const.n2,iseed,occ2)

		occ1=occupancy(const.n,const.ns,occ2,occ1)
		
		occ=np.zeros((3*const.n-1,3*const.n-1))
		occ=rocc(const.n,iseed,occ,occ1,const.nsq2)
		
		np.set_printoptions(formatter={'float_kind':float_formatter})
		np.savetxt(f3,occ,delimiter='\n',fmt='%f')
		np.set_printoptions()
		
		sten1=siteenergy(const.n,const.ns,occ,phi,sten1)
		
		toten,totenst=hamiltonian(const.n,const.ns,occ,phi)#,toten,totenst)
		
		for itemp in range(const.temp):
			beta=betaa[itemp]
			start=process_time()
			nmeas=nmeas1[itemp]
			nskip=nskip1[itemp]
			toten=float("{0:08.2f}".format(toten))
			f4.write("0\t%f\n" % (toten))		#not serialized. Read as string if ever, the file is read by any program. Can be serialized using pickle; retain datatype if read later.
			ntime=0
			totmin=0.0
			
			#monte Carlo Simulation
			
			for imeas in range(1,int(nmeas)+1):
				#test1=process_time()
				for iskip in range(1,int(nskip)+1):
					for nmcs in range(1,int(const.n2)+1):
						ie,je,ih,jh=choosesite(const.n,const.ns,iseed,occ1)#,ie,je,ih,jh)
						#print (ie,je,ih,jh,'\n')
						delst=delsten(sten1,ie,je,const.n,ih,jh,distinv)
						#print (delst,'\n')
						if(delst<0.0):
							occ1[ie][je]=-occ1[ie][je]
							occ1[ih][jh]=-occ1[ih][jh]
							toten+=delst
							sten1=updatehe(const.n,ie,je,ih,jh,sten1,distinv)
							ntime+=1
							if(toten<totmin):
								totmin=toten
								occmin=occ1
						else:
							a=np.exp(-delst*beta)
							c=np.random.rand()
							if(c<a):
								occ1[ie][je]=-occ1[ie][je]
								occ1[ih][jh]=-occ1[ih][jh]
								toten+=delst
								sten1=updatehe(const.n,ie,je,ih,jh,sten1,distinv)
								ntime+=1
								if(toten<totmin):
									totmin=toten
									occmin=occ1
				#test2=process_time()
				#print(test2-test1,' ',)

			#print ("ntime ", ntime, nmeas, nskip,'\n')
			
			toten=float("{0:014.4f}".format(toten))
			totmin=float("{0:014.4f}".format(totmin))
			f4.write("\t%f\t%f\n" % (toten,totmin))
			
			for p in range(const.n):
				for s in range(const.n):
					occint1[p][s]=int(occmin[p][s]+0.5)
					f8.write("%i\n" % (occint1[p][s]))
					
			finish=process_time()
			
			print ('k=',k,'\ttemp=',beta,'\ttime=',(finish-start),'\trun',ntime,'\ttoten=',toten)
	
			
			

