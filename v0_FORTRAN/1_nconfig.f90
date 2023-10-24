!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!                                                                                     !
!                 PROGRAM TO DETERMINE NUMBER OF CONFIGURATIONS                       ! 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  This program reads the output file of the main program and tells us the number of  !
!  configuration files. This program requires program_parameters.dat, output file     !
!                                                                                     !
!                          modified :: 21st August 2020                               !
!                                                                                     !
! Mohammad Samsuzzaman                                                                !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  
        program readfile
	    implicit none
           !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	    real*8 :: b,dt,F,boxlen
	    integer :: nlines,io,pnum,nfiles,time,interval
            character(len=30)::aa,bb,cc,dd,ee,ff,gg
	   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	    open(unit=1, file = 'output', action='read', position='rewind',form="unformatted")		
	    open(unit=4,file="nfiles.dat")
           !~~~~~~~~~~~~~~~~~~~~~ Values of the variables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	    pnum=1000	
	    nlines = 0
	   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           do
	      read(unit=1,iostat=io)  !READ(1,*,iostat=io) use this for formatted type
	      if(io/=0) exit
	         nlines = nlines + 1
	   end do
		         
	   nfiles=nlines/pnum
	   close(unit=1)
!	   write(*,*)"*Particle_Number*","* Nconfig*"
!          write(*,*)pnum, nfiles
           write(unit=4,fmt=*)nfiles
	  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   end program readfile
