!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!                                                                                     !
!                          OUTPUT FILE EXTRACTION PROGRAM                             ! 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  This program reads the output file of the main program and converts it into user   !
!  defined format, so that it can be read in any visualization software such as ovito !
!                                                                                     !
!                          modified :: 21st August 2020                               !
!                                                                                     !
! Mohammad Samsuzzaman                                                                !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
              
	program readfile
	implicit none
	!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	   real :: z
           real*8 :: b,dt,F,boxlen
           character(len=30)::aa,bb,cc,dd,ee,ff,gg
	   integer::pnum,nfiles,i,j,time,interval,nfiles_expected
	   real*8, allocatable,dimension(:)::xold,yold,vxold,vyold,theta,Fx	   
	   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~		   
	   open(unit=1, file = 'output',form="unformatted", action='read')
	   open(unit=2,file="nfiles.dat",form="formatted", action='read')
	   open(unit=4,file="sam.dat")
           !~~~~~~~~~~~~~~~~~~~~~ Values of the variables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           read(unit=2,fmt=*)nfiles        ! Total number of configurations
           pnum=1000
           allocate(xold(1:pnum),yold(1:pnum),vxold(1:pnum),vyold(1:pnum))
           allocate(theta(1:pnum),Fx(1:pnum))
           z=1.0d0
           !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                   
            do i=1,nfiles
                             
	         write(unit=4,fmt=*)pnum
                 write(unit=4,fmt=*)
                         		        
		do j=1,pnum
		      read(unit=1)xold(j),yold(j),vxold(j),vyold(j)!,Fx(j)
		      write(unit=4,fmt=*)xold(j),yold(j),vxold(j),vyold(j)!,Fx(j)
	       end do

            end do 
             
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
           end program readfile
