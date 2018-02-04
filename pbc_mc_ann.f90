! 28.9.2015
! MC annealing
! new pbc program to solve RAM problem
! changes are made in the distance program
! hartree energy subroutine included
! instead of delH we are using ei-ej-1/rij

PROGRAM mc_pbc

IMPLICIT NONE
INTEGER,PARAMETER::n=16,ndim=2,ns=n**ndim,n2=ns/2,nconfig=1,nsq2 = ns/2
INTEGER,PARAMETER::temp=21
INTEGER::iseed,i,j,x1,iy,ix,ie,je,ih,jh,imeas,nmeas,iskip,nskip,nmcs,k
INTEGER::it,itemp,occint1(0:n-1,0:n-1),ntime, nmeas1(temp),nskip1(temp)
REAL(KIND=8)::phi(0:n-1,0:n-1),occ(-n:2*n -1,-n:2*n -1),occ2(ns)
REAL(KIND=8)::occ1(0:n-1,0:n-1),toten,totenst,beta,a,c,betaa(temp)
REAL(KIND=8)::distinv(0:n-1,0:n-1),delen,ran2,totmin,occmin(0:n-1,0:n-1)
REAL(KIND=8)::sten1(0:n-1,0:n-1),delst,start,finish,totenh

	iseed = -1234567


	DO it = 1,temp
	   READ(1,*) betaa(it)
	   READ(7,*) nmeas1(it),nskip1(it)
	ENDDO
	
	CALL distpbc(n,distinv)

	DO k = 1,nconfig
	
	CALL rsiteen(n,iseed,phi)

	do i=0,n-1
	   do j=0,n-1
		WRITE(2,46) phi(i,j)
46      FORMAT (1X,F8.5)
	   enddo
	enddo


	CALL intialocc(ns,n2,iseed,occ2)

	CALL occupancy(n,ns,occ2,occ1)
	print *,occ1

	CALL rocc(n,iseed,occ,occ1,nsq2)

	do i=0,n-1
	   do j=0,n-1
		WRITE(3,47) occ(i,j)
47      FORMAT (1X,F8.5)
	   enddo
	enddo

	CALL siteenergy(n,ns,occ,phi,sten1)

	CALL Hamiltonian(n,ns,occ,phi,toten,totenst)


	DO itemp = 1,temp
	   beta = betaa(itemp)
	call cpu_time (start)
	   nmeas = nmeas1(itemp)
	   nskip = nskip1(itemp)

	write(4,51) 0, toten
51      FORMAT (1X,I1,1X,F8.2)
	ntime=0
	totmin=0.0
!====================================================
! Monte Carlo Simulation

	DO imeas = 1,nmeas
	   Do iskip = 1,nskip
	      do nmcs = 1,n2
	   
		CALL choosesite(n,ns,iseed,occ1,ie,je,ih,jh)

		CALL delsten(sten1,ie,je,n,ih,jh,delst,distinv)
		
				if (delst .le. 0.0) then
					occ1(ie,je) = -occ1(ie,je)
					occ1(ih,jh) = -occ1(ih,jh)
					toten=toten+delst

				call updatehe(n,ie,je,ih,jh,sten1,distinv)
				ntime=ntime+1
					if (toten .lt. totmin ) then
						totmin=toten
  						occmin=occ1
					endif

				else
					a=exp(-delst*beta)
					c=ran2(iseed)

					if (c .le. a) then
						occ1(ie,je) = -occ1(ie,je)	
						occ1(ih,jh) = -occ1(ih,jh)
						toten=toten+delst

                                        call updatehe(n,ie,je,ih,jh,sten1,distinv)	    
					ntime=ntime+1
					if (toten .lt. totmin ) then
						totmin=toten
  						occmin=occ1
					endif

					endif
				endif
	      enddo				! nmcs loop
	   ENDDo				! iskip loop
	ENDDO					! imeas loop

	write(4,22) toten,totmin
22      FORMAT (1X,F14.4,1X,F14.4)

	do i=0,n-1
	   do j=0,n-1
		occint1(i,j) = INT(occmin(i,j) + 0.5)
	        write(8,23) occint1(i,j)
23              FORMAT (1X,I1)
	   enddo
	enddo

		call cpu_time (finish)

!	print*, "k=",k,"temp",beta,"time",(start-finish)/60.0,"run",ntime,"toten=",toten

	ENDDO					! temp loop

	ENDDO					! configuration loop

END PROGRAM mc_pbc

! ============================================
! disorder [-0.09,0.09]
SUBROUTINE rsiteen(n,iseed,phi)
IMPLICIT NONE
INTEGER:: n,iseed,i,j			!iseed purpose
REAL(KIND=8):: phi(0:n-1,0:n-1),RAN2	!RAN2 purpose?

	do i=0,n-1
	   do j=0,n-1
	       phi(i,j)=0.0

	   enddo
	enddo

RETURN
END
!==========================================
! random initial occ [e:0.5 , h:-0.5]
SUBROUTINE intialocc(ns,n2,iseed,occ2)
IMPLICIT NONE
INTEGER::n,ns,n2			!randomly making half of lattice sites -0.5 and other half is already 0.5?
INTEGER::i,j,iseed
REAL(KIND=8)::occ2(ns),ran2
	
	do i=1,ns
	   occ2(i)=.5
	enddo

	DO i = 1, n2
11	   j = INT(1 + ns * RAN2(ISEED))
	   IF (occ2(j).EQ.-.5)GOTO 11
	   occ2(j) = -.5
	ENDDO

RETURN
END
! =============================================================
! coordinates form
SUBROUTINE occupancy(n,ns,occ2,occ1)
IMPLICIT NONE
REAL(KIND=8):: occ1(0:n-1,0:n-1),occ2(ns)
INTEGER::i,j,ix,iy,jx,jy,n,ns,x1,j1

   DO i = 1,ns
      x1 = mod(i,n)
      ix = x1-1
      if(x1 .eq. 0) then
       ix = n-1
       iy = i/n-1
      else
       iy = i/n
      endif
      occ1(ix,iy) = occ2(i)
   ENDDO

RETURN
END

!================================================================
! minimum images
	subroutine rocc(n,iseed,occ,occ1,nsq2)
	IMPLICIT NONE
	INTEGER::n,iseed,nsq2,i,j,nx,ny,ix,jx,kkk	
	REAL(KIND=8):: occ(-n:2*n -1,-n:2*n -1),RAN2
	REAL(KIND=8):: occ1(0:n-1,0:n-1)
		
	do i = 0,n-1
	do j = 0,n-1
	
	do nx = -1,1
	do ny = -1,1
	
	
	   ix = i + nx * n
	   jx = j + ny * n
	   
	   occ(ix,jx) = occ1(i,j)
	
	enddo
	enddo
	
	enddo
	enddo
		
RETURN 
END	

!=============================================================
! calculation of Hamiltonian

subroutine Hamiltonian(n,ns,occ,phi,toten,totenst)
IMPLICIT NONE
INTEGER::n,nrc,i,j,ix,iy,jx,jy,ns
REAL(KIND=8):: phi(0:n-1,0:n-1),occ(-n:2*n -1,-n:2*n -1)
REAL(KIND=8):: toten,totenst,toten1,dist1,sssssss 

nrc = n/2
toten =0.0

	do i=0,n-1
	  do j=0,n-1
	     toten=toten+phi(i,j)*occ(i,j)
	  enddo
	enddo

toten1 =0.0

	do ix = 0,n-1
	  do iy = 0,n-1
	  	do jx = ix - nrc+1 ,ix + nrc
	       do jy = iy - nrc+1 ,iy + nrc


		if ((ix .eq. jx) .and. (iy .eq. jy)) goto 45
		  dist1 = sqrt (real((ix - jx)**2) + real((iy - jy)**2))
		  toten1=toten1+( (occ(ix,iy)*occ(jx,jy)) / dist1)

45 sssssss=1.0

		enddo
	      enddo
            enddo
	enddo
toten=toten+toten1/2.0
totenst=toten/ns

return
end
!=======================================================================
! hartree energy calculation

subroutine siteenergy(n,ns,occ,phi,sten1)
IMPLICIT NONE
INTEGER::n,nrc,i,j,ix,iy,jx,jy,ns
REAL(KIND=8):: phi(0:n-1,0:n-1),occ(-n:2*n -1,-n:2*n -1)	
REAL(KIND=8):: sten1(0:n -1,0:n -1),sssssss,dist1

	

	do ix = 0,n-1
	do iy = 0,n-1

	sten1(ix,iy) = phi(ix,iy)

	nrc = n/2

	do jx = ix - nrc+1 ,ix + nrc
		do jy = iy - nrc+1 ,iy + nrc
	
	if ((ix .eq. jx) .and. (iy .eq. jy)) goto 4

	dist1 =  sqrt (real((ix - jx)**2) + real((iy - jy)**2))

   	sten1(ix,iy)=sten1(ix,iy)+( (occ(jx,jy)) / dist1)

4	sssssss=1.0

	enddo
	enddo
	
	enddo
	enddo
	
	return
	end	
!=======================================================================
! randomly choose 2 sites for exchange 
subroutine choosesite(n,ns,iseed,occ1,ie,je,ih,jh)
IMPLICIT NONE
INTEGER::n,iseed,ie,je,ih,jh,kk,ns,ic,x1,x2,kc
REAL(KIND=8)::occ1(0:n-1,0:n-1),ss,sss,ran2


	do kk = 1,n**4
        	ic=INT(ran2(iseed)*ns+1)
		      x1 = mod(ic,n)
		      ie = x1-1
		      if(x1 .eq. 0) then
		       ie = n-1
		       je = ic/n-1
		      else
		       je = ic/n
		      endif

	if (occ1(ie,je) .gt. 0.0)then
	goto 46
	endif
	enddo		!jumps the loop? what's happening in this if clause?same in next loop below.

46      ss=1.0

	do kk = 0,n**4
	        kc=INT(ran2(iseed)*ns+1)
		      x2 = mod(kc,n)
		      ih = x2-1
		      if(x2 .eq. 0) then
		       ih = n-1
		       jh = kc/n-1
		      else
		       jh = kc/n
		      endif

	if (occ1(ih,jh) .lt. 0.0) then
	goto 47
	endif
	enddo

47	sss=1.0

return
end

!=======================================================================
! distance calculation using pbc
!=============================================================
! calculates 1/r term using pbc
subroutine distpbc(n,distinv)
IMPLICIT NONE
INTEGER::n,k,j,itest,ix,iy
REAL(KIND=8)::distinv(0:n-1,0:n-1),dpbc(0:n-1,0:n-1)

	DO ix = 0,n-1
	DO iy = 0,n-1

	   k = ix
	   j = iy

	If (ix .gt. n/2) then
	   k=abs(ix-n)
	endif
	If (iy .gt. n/2) then
	    j=abs(iy-n)
	endif

	dpbc(ix,iy) = sqrt(real(k*k)+ real(j*j))
	if ((ix .eq. 0.0) .and. (iy .eq. 0.0)) then
	    distinv(ix,iy) = 0.0
	else
	    distinv(ix,iy) = 1/dpbc(ix,iy)
	endif

	ENDDO
	ENDDO

return
end
!=======================================================================

!=======================================================================
! change in hartree energy calculation

subroutine delsten(sten1,ie,je,n,ih,jh,delst,distinv)
IMPLICIT NONE
REAL(KIND=8):: sten1(0:n -1,0:n -1),delst,dpbc,ssss
INTEGER::n,ih,jh,ie,je,ix,iy
REAL(KIND=8)::distinv(0:n -1,0:n -1)

	ix=abs(ie-ih)
	iy=abs(je-jh)
	
	delst = sten1 (ih,jh) - sten1 (ie,je)
	
	if ((ih .eq. ie) .and. (jh .eq. je)) goto 20

	delst = delst - distinv(ix,iy) 

20	ssss=1.0				!are s* having anyuse apart from halping gotos?
	
	return
	end		

!===============================================================
! updates the images after each exchange
subroutine updateocc(n,ie,je,ih,jh,occ)  		!used where?
IMPLICIT NONE
INTEGER::n,ie,je,ih,jh,nx,ny,ixe,jxe,ixh,jxh
REAL(KIND=8):: occ(-n:2*n -1,-n:2*n -1)

	do nx = -1,1
	  do ny = -1,1

	     ixe = ie + nx*n
	     jxe = je + ny*n
	     ixh = ih + nx*n
	     jxh = jh + ny*n

	     occ(ixe , jxe ) = occ(ie,je)
	     occ(ixh , jxh ) = occ(ih,jh)

	   enddo
	  enddo

return
end
!=============================================================

! update of hartree energy after each exchange.
SUBROUTINE updatehe(n,ie,je,ih,jh,sten1,distinv)
IMPLICIT NONE
INTEGER::n,ie,je,ih,jh,ix,iy,jx,jy,nrc,i,j,nx,ny
INTEGER::ik,ij,il,im
REAL(KIND=8)::distinv(0:n -1,0:n -1),sssssss
REAL(KIND=8):: sten1(0:n-1,0:n-1)

nrc = n/2

	do ix = 0,n-1
	do iy = 0,n-1

	ik=abs(ix-ih)
	ij=abs(iy-jh)
	il=abs(ix-ie)
	im=abs(iy-je)

		if ((ix .eq. ie) .and. (iy .eq. je)) then
		  sten1(ix,iy)=sten1(ix,iy)+distinv(ik,ij)

		endif
		
		if  ((ix .eq. ih) .and. (iy .eq. jh)) then
		  sten1(ix,iy)=sten1(ix,iy)-distinv(il,im)

		endif
	

           if (((ix .eq. ie) .and. (iy .eq. je)) .or. ((ix .eq. ih) .and. (iy .eq. jh))) goto 4

	   	sten1(ix,iy)=sten1(ix,iy)+(distinv(ik,ij)-distinv(il,im))

4	sssssss=1.0

	enddo
	enddo
	
RETURN
END

!=============================================================================
! Random number generator.

	  DOUBLE PRECISION FUNCTION ran2(idum)
!	  FUNCTION ran2(idum)
	  INTEGER :: idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
	  REAL(KIND=8) :: AM,EPS,RNMX
	  PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1)
          PARAMETER(IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791)
          PARAMETER(NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
	  INTEGER :: idum2,j,k,iv(NTAB),iy
	  SAVE iv,iy,idum2
	  DATA idum2/123456789/, iv/NTAB*0/, iy/0/
	  IF (idum.LE.0) THEN
	  idum=MAX(-idum,1)
	  idum2=idum
	  DO 11 j=NTAB+8,1,-1
	    k=idum/IQ1
	    idum=IA1*(idum-k*IQ1)-k*IR1
	    IF (idum.LT.0) idum=idum+IM1
	    IF (j.LE.NTAB) iv(j)=idum
  11      CONTINUE
	  iy=iv(1)
	  ENDIF
	  k=idum/IQ1
	  idum=IA1*(idum-k*IQ1)-k*IR1
	  IF (idum.LT.0) idum=idum+IM1
	  k=idum2/IQ2
	  idum2=IA2*(idum2-k*IQ2)-k*IR2
	  IF (idum2.LT.0) idum2=idum2+IM2
	  j=1+iy/NDIV
	  iy=iv(j)-idum2
	  iv(j)=idum
	  IF(iy.LT.1)iy=iy+IMM1
	  ran2=MIN(AM*iy,RNMX)
	  ran2=DBLE(ran2)

	  RETURN
	  END
