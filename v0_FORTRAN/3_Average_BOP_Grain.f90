!************************************************************************************
!///////////////////////////// BOND ORDER PARAMETER ////////////////////////////////
!//////////////////////////////////////////////////////////////////////////////////
!************************************************************************************


        program bondorderparameter
        implicit none
      
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        
	integer :: pnum,i,j,k,num,nfiles,initial_nfiles 
	real*8  :: r0,pi,xij,yij,theta
        real*8  :: bopp,sinsum,cossum,r,bop_thresh,v0,stef,vc
	real*8, allocatable,dimension(:)::x,y,vx,vy,bop,bopavg	
	character(len=1024) :: filein,format1,fileout,format2,format3
	character(len=1024) :: fileout1        
	r0= 1.220d0
	pnum=1000
	bop_thresh=0.80d0
	pi=4.0d0*atan(1.0d0)
        initial_nfiles = 1000 ! No of files that has to be ignored 
	!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	allocate(x(1:pnum),y(1:pnum),vx(1:pnum),vy(1:pnum),bop(1:pnum))
	open(unit=1, file = 'output',form="unformatted", action='read')
        open(unit=2,file="nfiles.dat",form="formatted", action='read')
!	open(unit=3,file="bop.dat")
        open(unit=5,file="program_parameters.dat")

	open(unit=4,file="avg_bop.dat")
	read(unit=2,fmt=*)nfiles
	allocate(bopavg(1:nfiles))
        read(unit=5,fmt=*)vc,stef,v0
       !~~~~~~~Read first 50 files ie till system reaches steady state ~~~~~~~~~~~~~
       !                  so read 50 files but ignore it
        
        do k= 1,initial_nfiles    
      
          do i=1,pnum
             read(unit=1)x(i),y(i),vx(i),vy(i)
          end do
       
        end do
    
      !    write(*,*)"First 50 files ignored"
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	do k=initial_nfiles+1,nfiles                  ! Loop over all the configuration files 
                    
!                write(*,*)k
	            do i=1,pnum
                      read(unit=1)x(i),y(i),vx(i),vy(i)
                    enddo  
	
	          
	         do i=1,pnum       ! 1 to particle number

	             num=0
                       sinsum=0.0d0
                       cossum=0.0d0

                      do j=1,pnum
                        if(j.ne.i)then
                           xij=x(i)-x(j)
                           yij=y(i)-y(j)
                           r=sqrt(xij*xij+yij*yij)
                           if(r.lt.r0)then
                             num=num+1
                             theta=datan(yij/xij)
                        
                             if(xij.lt. 0.0d0)theta=theta+pi
                        
                             theta=theta*6.0d0
                             cossum=cossum+cos(theta)
                             sinsum=sinsum+sin(theta)
                          endif
                        endif
                      enddo   ! j loop END
 
                      bopp=sqrt(cossum**2.0d0+sinsum**2.0d0)

                      if(num.eq.0)then
                         bop(i)=bopp
                      else
                         bop(i)=bopp/num
                      endif

!                      write(unit=2,fmt=*)x(i),y(i),bop(i),num
           

                   enddo  ! i loop END 
                   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                 write(unit=3,fmt=*)pnum
!                  write(unit=3,fmt=*)
!	        do i=1,pnum
!                     !if(bop(i).le.bop_thresh)then
!                        write(unit=3,fmt=*)x(i),y(i),bop(i)
!                     !endif
!	        enddo
                  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                bopavg(k) = (sum(bop))/pnum
!                write(unit=4,fmt=*)bopavg	     

	enddo      ! k loop END
	!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~ Printing to File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   1000 format (4x,A,10x,A,5x,A)   
!   write(unit=4,fmt=1000)"External Force","Interparticle Friction","Diffusion-coefficient"
   write(unit=4,fmt=*)vc,stef,v0,sum(bopavg)/(nfiles-initial_nfiles)
   write(*,*)"Done"                    ! particle average 
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   

	end program bondorderparameter
	
