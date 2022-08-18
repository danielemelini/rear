!
!  Copyright Daniele Melini & Giorgio Spada, 2014-2022
!
!  This file is part of REAR.
!
!  REAR is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  REAR is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with REAR.  If not, see <http://www.gnu.org/licenses/>.
!
! ------------------------------------------------------------------
!
! File "hrm.f90" contains the following program units:
!
!       - Subroutine "FINDPX":
!	       Coordinates of the pixels of the Tegmark grid        		
!              
!       - Functions "SINDD" and "COSDD"
!              Real*8 cosine and sine (argument in degres)
!
!       - subroutine "count_rows":
!              Counts the number of rows in a text file
!
!       - subroutine "count_headers":
!              Counts the number of header lines in a text file
!              (i.e., those lines beginning with '!' or '#')
!
!       - subroutine "read_data_line":
!              Reads the next non-commented line from a text file
!              (i.e., the next line which do not begins with
!              '#' or '!')
!
! Last review by DM on Aug 17, 2022. 
!
! ------------------------------------------------------------------
!
!
!
!
!
!
	Subroutine FINDPX(RES,NP,LONP,LATP)
!
	implicit NONE 
!
! # Given resolution RES, this routine returns the number of pixels NP,
!   and the longitude and latitude arrays LONP, LATP (1:NP) with the pixels
!   coordinates. THis done through a call to "pixel2vector" by M. Tegmark- 
!   
!   Written by GS on October 2, 2008 for the implementation of the 
!   3D velocity maps for SELEN 2.7 - but perhaps useful also for other 
!   purposes... 
!
	INTEGER J, NP, RES
        REAL*4 R(0:19,3,3), V(0:11,3), VECT(3), X,Y,Z  
	REAL*4 LONP(1:NP), LATP(1:NP), LONPX, LATPX  
	REAL*4,  PARAMETER :: PI=3.14159265358979323840 
!
        call compute_matrices (r)
        call compute_corners  (v)
!
	np=40*res*(res-1)+12 
!
        do 1 j=0,np-1
	  	call pixel2vector (j,res,r,v,vect)
	        x=vect(1)
		y=vect(2)
		z=vect(3)
! 
! --- Polar pixel 
          	if(x==0..and.y==0.) then 
          		lonpx=0.
                	if(z>=0.) latpx = +90.
                	if(z<=0.) latpx = -90.
! --- Ordinary pixel 
	  		else                         	  
	  		lonpx = atan2(y,x)*180./pi 
          		if (lonpx < 0.) lonpx=360.+ lonpx
          		latpx = 90.-acos(z)*180./pi 
          	endif
	 lonp(j+1)=lonpx 
	 latp(j+1)=latpx 	
!
1       continue 
!
	end subroutine FINDPX 
!
!
!
!
!
!
!
!
	FUNCTION SINDD(ALFA)
	implicit NONE
!
! --- Assuming that ALFA is given in DEGREES, it returns the SINE of ALFA
!     *********************** GS and FC 11-02-2009 **********************
!	
	REAL*8,  PARAMETER :: PI=3.14159265358979323840 
	REAL*8 SINDD, ALFA, ALFARAD              
	ALFARAD=ALFA*PI/180.	
	SINDD=SIN(ALFARAD)
	
	END FUNCTION SINDD 
!
!
!
!
!
!
	FUNCTION COSDD(ALFA)
	implicit NONE
!
! --- Assuming that ALFA is given in DEGREES, it returns the COSINE of ALFA
!     ************************ GS and FC 11-02-2009 ***********************
!	
	REAL*8,  PARAMETER :: PI=3.14159265358979323840 
	REAL*8 COSDD, ALFA, ALFARAD               
	ALFARAD=ALFA*PI/180.	
	COSDD=COS(ALFARAD)	
!
	END FUNCTION COSDD
!
!
!
!
!
    subroutine count_rows(lun,n)
    implicit none
!
! --- Returns as N the number of rows in text file opened at LUN
!     DM 28.04.2014
!
    integer lun,n
!
    rewind(lun)
!    
    n=0
901 read(lun,*,end=902)
    n=n+1
    go to 901
    902 continue
!
    return
!
    end subroutine count_rows
!    
!
!
!
!
    subroutine count_headers(lun,n)
    implicit none
!
! --- Returns as N the number of header lines in text file opened at LUN
!     Header lines begin with '#' or '!'
!     DM 01.05.2022
!
    integer lun,n
    character :: c
!
    rewind(lun)
!    
    n=0
901 read(lun,*,end=902) c
    if ( ( c.ne.'#' ) .and. ( c.ne.'!' ) ) go to 902
    n=n+1
    go to 901
    902 continue
!
    return
!
    end subroutine count_headers
!    
!
!
!
    subroutine read_data_line(lun,buffer)
    implicit none
    character(*)   :: buffer
    character :: ch
    integer :: lun
    integer :: idx
!
910 read(lun,'(a)') buffer
    ch = buffer(1:1)
    if ((ch.eq.'#').or.(ch.eq.'!')) go to 910
!
    if( (index( buffer, '#' ).ne.0) .or. (index( buffer,'!' ).ne.0) ) then
! 
       if( index(buffer,'#').eq.0 ) then
           idx = index(buffer,'!')
       elseif( index(buffer,'!').eq.0 ) then
           idx = index(buffer,'#')
       else
           idx = min( index( buffer, '#' ), index( buffer, '!' ) )
       endif	   
!
       buffer = buffer( 1:(idx-1) )
!
    end if
!
    return
!
    end
!
!    
