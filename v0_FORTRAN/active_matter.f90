



!****************************************************************************************

!                      active particles in cirular box

!****************************************************************************************
program test
implicit none
character(len=70) :: fn
real*8, dimension (1:5000) ::x,y,vx,vy,dumx,dumy,direcx,direcy,vxdum,vydum,dumvx,dumvy
real*8:: eta,L,l0,v0,t,dt,tau,r,r1,pi,x2,y2,Xji,Yji,a,theta,bound,alpha,fric,cosrv,vc,v,radRV,r2d
real*8:: eff,a0,er1,er2,f1,fpot,f3x,f3y,direcabs,stef,hig,KBT,gasdev,ran1
integer:: i,j,N,nn,pp,idum

real*8:: rr_2,rad_theta,rad_plus,rad_minus, vx1, vy1, vx2, vy2

f1(r)=4.0d0*eff*(12.0d0*(a/r)**12.0d0 - 6.0d0*(a/r)**6.0d0)/r
fpot(r)=-hig*stef/((cosh(stef*(r-bound)))**2.0d0)


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
stef=0.10d0 ! stiffness
v0=0.40d0     ! activity
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

N=1000              !number of particles
tau=80             ! time
idum=967369025
alpha=1.0d0
fric=10.0d0           ! gamma (friction)
KBT=0.00010d0         ! Temperature 
a=1.0d0               ! Interparticle interaction potential parameter 
a0=a*(2.0d0**(1.0d0/6.0d0))   ! cut off for interparticle interaction (WCA)  
eff=1.0d0             ! Interparticle interaction potential parameter
hig=100.0d0           ! Trap height          
bound=80.0d0!25.819          ! Trap radius  
eta=1.0d0             ! Angular noise span 
L=bound-2.0d0         ! Initial condition parameter, so that particles are inside
l0=2.0d0               ! Vicsek radius
dt=0.001d0            ! Time increment
t=0.0d0               ! Initial Time
pi=4.0d0*atan(1.0d0)       ! pi value
vc=pi*(1.0d0/6.0d0)   ! vision cone angle


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


 
 		open(unit=1,file='output',form='unformatted')
 	
 	open(unit=8,file="program_parameters.dat",form="formatted")
 	
 	write(unit=8,fmt=*)vc,stef,v0

i=1
      call random_number(r1)
      theta=2.0d0*pi*r1

      call random_number(r1)
      x(i)=L*r1*cos(theta)
      y(i)=L*r1*sin(theta)
      call random_number(r1)
      theta=2.0d0*pi*r1
      vx(i)=cos(theta)
      vy(i)=sin(theta)
      i=i+1


10 if(i.le.N)then
   write(*,*)i
      call random_number(r1)
      theta=2.0d0*pi*r1

      call random_number(r1)
      x(i)=L*r1*cos(theta)
      y(i)=L*r1*sin(theta)

      do j=1,i-1
            x2=x(i)-x(j)
         
            y2=y(i)-y(j)
        
            r=sqrt(x2*x2  + y2*y2)
            if(r.lt.a0) go to 10
      end do

      call random_number(r1)
      theta=2.0d0*pi*r1
      vx(i)=cos(theta)
      vy(i)=sin(theta)
!      write(*,*)i,x(i),y(i)
      i=i+1
      go to 10
   end if





nn=1
pp=1
20 if(t.lt.tau)then
   
   ! write(*,*)pp
    pp=pp+1

     !$OMP PARALLEL
     !$OMP DO PRIVATE(i)
     do i=1,N
       direcx(i)=vx(i)
       direcy(i)=vy(i)
       v=sqrt(vx(i)*vx(i) + vy(i)*vy(i))
       f3x=0.0d0
       f3y=0.0d0
!       write(*,*)t,i,vy(i),vx(i),atan(vy(i),vx(i)),theta(i)
       
       do j=1,N
     
           if(j.ne.i)then
             
            x2=x(i)-x(j)
            y2=y(i)-y(j)

            r=sqrt(x2*x2  + y2*y2)
	   
            


! Direction choose
  if(r.lt.l0)then

      Xji=-1.0d0*x2           
      Yji=-1.0d0*y2
      
      cosrv=((Xji*vx(i))+(Yji*vy(i)))/(r*v)
      
      if(cosrv.gt.cos(vc))then

               direcx(i)=direcx(i)+vx(j)
               direcy(i)=direcy(i)+vy(j)
            
      endif             
      endif   

! WCA interaction
 
            if(r.lt.a0)then
               f3x=f3x+f1(r)*(x2/r)
               f3y=f3y+f1(r)*(y2/r)
            endif

          endif
	end do
        
        direcabs =sqrt(direcx(i)*direcx(i) + direcy(i)*direcy(i))
        direcx(i)=direcx(i)/direcabs
        direcy(i)=direcy(i)/direcabs

        call random_number(r1) 
        r1=r1-0.50d0
        r1=r1*eta

        vxdum(i)=(direcx(i)*cos(r1)-direcy(i)*sin(r1))
        vydum(i)=(direcx(i)*sin(r1)+direcy(i)*cos(r1))
       
        r=sqrt(x(i)*x(i) + y(i)*y(i))

        dumvx(i)=-fric*vx(i)*dt+v0*vxdum(i)*dt+alpha*(fpot(r)*x(i)/r+f3x)*dt+sqrt(2.0d0*fric*KBT*dt)*gasdev(idum)
        dumvy(i)=-fric*vy(i)*dt+v0*vydum(i)*dt+alpha*(fpot(r)*y(i)/r+f3y)*dt+sqrt(2.0d0*fric*KBT*dt)*gasdev(idum)


        dumx(i)=vx(i)*dt
        dumy(i)=vy(i)*dt
             
      end do

      !$OMP END DO
      !$OMP END PARALLEL

      do i=1,N

         vx(i)=vx(i)+dumvx(i)
         vy(i)=vy(i)+dumvy(i)

         x(i)=x(i)+dumx(i)
         y(i)=y(i)+dumy(i)         

      end do

      t=t+dt
      

       if(t.gt.(nn*2.0d0))then
      
      	!write(*,*)nn
        
         do i=1,N           
        	 write(unit=1)x(i),y(i),vx(i),vy(i)
	 end do
         
          !close(unit=1000+nn)
          nn=nn+1
      end if
 
      go to 20
   endif
   




end program test

	!********************************************
	!Gaussian Random Number generation
	!********************************************

	FUNCTION gasdev(idum)
	INTEGER idum
	real*8 gasdev
!      	USES ran1
	INTEGER iset

	real*8 fac,gset,rsq,v1,v2,ran1
	SAVE iset,gset
	DATA iset/0/
	if (iset.eq.0) then
1       v1=2.*ran1(idum)-1.
	v2=2.*ran1(idum)-1.
	rsq=v1**2+v2**2
	if(rsq.ge.1..or.rsq.eq.0.)goto 1
	fac=sqrt(-2.*log(rsq)/rsq)
	gset=v1*fac
	gasdev=v2*fac
	iset=1
	else
	gasdev=gset
	iset=0
	endif
	return
	END

	FUNCTION ran1(idum)
	INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
	real*8 ran1,AM,EPS,RNMX
	PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836)
	PARAMETER (NTAB=32)
	PARAMETER (NDIV=1+(IM-1)/NTAB)
	PARAMETER (EPS=1.2e-7)
	PARAMETER (RNMX=1.-EPS)

	INTEGER j,k,iv(NTAB),iy
	SAVE iv,iy
	DATA iv /NTAB*0/, iy /0/
	if (idum.le.0.or.iy.eq.0) then
	idum=max(-idum,1)
	do 11 j=NTAB+8,1,-1
	k=idum/IQ
	idum=IA*(idum-k*IQ)-IR*k
	if (idum.lt.0) idum=idum+IM
	if (j.le.NTAB) iv(j)=idum
11      continue
	iy=iv(1)
	endif
	k=idum/IQ
	idum=IA*(idum-k*IQ)-IR*k
	if (idum.lt.0) idum=idum+IM
	j=1+iy/NDIV
	iy=iv(j)
	iv(j)=idum
	ran1=min(AM*iy,RNMX)
	return
	END
