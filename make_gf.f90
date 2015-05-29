!
!  Copyright Giorgio Spada & Daniele Melini, 2014
!
!  This file is part of REAR (v. 1.0)
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
! ------------------------------------------------------------------------
!                                                                         
! Computes the "Green functions" for a disc-shaped surface load  
!                    
! Notes:   
!    * Requires the load-deformation coefficients for a SNREI Earth model
!    * Dependencies are included from file "hrm.f90"
!    * The LDCs are computed on 1D grid with uniform or varying spacing 
!    * See file "input_data_gf.inc" for configuration and settings
!
! Output data:
!    - vertical and horizontal displacement
!    - geoid height displacement 
!
! Copyright (C) Giorgio Spada (UNIURB) and Daniele Melini (INGV), 2013-15
!
! History: 
! First written by GS in 2011 for the GJI paper on Greenland 
! Updated several times until July 2013 
! Updated GS August 2013 for the "Wang" 1D grid implementation 
! Updated GS August 20 2013 (meaning of the "h" constant)
! Updated DM April 28, 2014 (OpenMP)
! Exits with an error for inconsistent ID=2 grids - DM Feb 5, 2015
!
! ------------------------------------------------------------------------  
!
! Declarations 
!
 include "hrm.f90"
!
 PROGRAM GREEN_FUNCTION
 implicit NONE
!
 include "input_data_gf.inc"
! 
! --- Legendre Pols
 real*8 leg(0:lmax+1)   ! Legendre polynomial  
 real*8 dleg            ! Derivative of leg wrt colatitude
!
! --- Load function 
 real*8 sigma (0:lmax)  ! Spectral amplitude of the load 
!
! --- Love numbers 
 real*8 h_love(0:lmax)  ! "h" Load deformation coefficient
 real*8 l_love(0:lmax)  ! "l" Load deformation coefficient
 real*8 k_love(0:lmax)  ! "k" Load deformation coefficient
!
! --- Miscellanea variables and constants  
 character*10 header
 integer i, l, ij 
 real*8, parameter :: from_m_to_mm=1D3 
 real*8  cosdd, ampl, r, q, x
!
! --- Dynamic arrays
 real*8, allocatable :: theta(:)
 real*8, allocatable :: u(:), v(:), g(:)
 integer :: ngr
! 
! ----------------------------------------------------------  
!
!
!
 Write(*,*) "make_gf: Reading the Love numbers from file: ", trim(adjustl(file_love_numbers))
!
  	open(1,file=file_love_numbers,status='unknown') 
 	do i=1, nh
 		read(1,'(a10)')header 
 	enddo 
!if(lmax.le.1024) then 
 		do i=lmin, lmax
                l=i 
                read(1,*) ij, h_love(l), l_love(l), k_love(l)   
                enddo 	
!endif

!	if(lmax.le.1024) then 
! 		do i=0, lmax
!                l=i 
!                read(1,*) ij, h_love(l), l_love(l), k_love(l)   
!                enddo 	
!        endif
!
!	if(lmax.gt.1024) then 
! 		do i=0, 1024
! 		l=i 
! 		read(1,*) ij, h_love(l), l_love(l), k_love(l)   
! 	enddo 	
!        	do l=1025, lmax
!	        h_love(l)= -6.2129D0
!	        l_love(l)=  1.8889D0 / float(l)
!	        k_love(l)= -3.0564D0 / float(l)
!	enddo
!	endif
!			 
!
!
 Write(*,*) "make_gf: Computing the harmonic coefficients of the load ..."
!
   call  plegendre(leg, lmax+1, cosdd(alfa))
!
!   do l=lmin, lmax
! 	  if(l==0) sigma(l)= 1d0-cosdd(alfa)
! 	  if(l>=1) sigma(l)= leg(l-1)-leg(l+1)
!          sigma(l)=sigma(l)/2d0 
!   enddo
!
!
   do l=lmin, lmax
 	  if(l==0) sigma(l)= 0D0
 	  if(l>=1) sigma(l)= -(leg(l+1)-leg(l-1))/(1d0+cosdd(alfa))
   enddo
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
 Write(*,*) "make_gf: Computing displacements and geoid height at the grid points ... "
! 
!
     IF    (GRID_OPT==1) THEN 
    	 Write(*,*) "make_gf: The Green function grid has a uniform spacing"
     ELSEIF(GRID_OPT==2) THEN
    	 Write(*,*) "make_gf: The Green function grid has an increasing spacing"
  	     R = 2d0*(theta_max-theta_min)/float(ngrid-1)/spac_min - 1d0 
  	     Q = spac_min
         IF (R .LT. 1) THEN
             Write(*,*) "make_gf: ERROR: Inconsistent grid"  
             Write(*,*) "make_gf:    The grid must satisfy spac_min < (theta_max - theta_min)/(ngrid - 1)."
             Write(*,*) "make_gf:    Either set spac_min < ", &
                           (theta_max-theta_min)/float(ngrid-1)
             Write(*,*) "make_gf:    or ngrid < ", &
                           nint( (theta_max-theta_min)/spac_min ) + 1
             stop
         ENDIF
     ENDIF
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
     IF     (GRID_OPT==1) then 
            ngr = nint( (theta_max - theta_min) / theta_inc ) + 1
            allocate( theta(ngr), u(ngr), v(ngr), g(ngr) )
            do ij=1,ngr
               theta(ij) = theta_min + (ij-1)*theta_inc
            end do        
     ELSEIF (GRID_OPT==2) then 
            ngr = ngrid
            allocate( theta(ngr), u(ngr), v(ngr), g(ngr) )            
            theta(1)=theta_min
            do ij=2,ngr
               theta(ij) = theta(ij-1) +  Q*( 1D0 + (R-1D0)*float(ij-2)/float(ngrid-2) )
            end do
     ENDIF
!     
     u = 0d0 
     v = 0d0
     g = 0d0
!     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP    PRIVATE(IJ,LEG,L,X,DLEG,AMPL) &
!$OMP    SHARED(NGR,THETA,SIGMA,U,V,G,H_LOVE,K_LOVE,L_LOVE) &
!$OMP    SCHEDULE(GUIDED)
     do ij=1,ngr
!
     call  PLegendre(leg, lmax+1, cosdd(theta(ij)))  
!   
     do l=lmin, lmax 
!   
        x = cosdd(theta(ij))
!        
        if(abs(x)==1d0) dleg  = 0d0 
        if(abs(x)/=1d0) dleg= - (l+1)*(x*leg(l)-leg(l+1))/((1-x)*(1+x))*sqrt(1-x**2)            
!
        ampl = sigma(l)/(2d0*float(l)+1d0)
!
   		u(ij)  =  u(ij) +        h_love(l) * ampl * leg(l)
   		g(ij)  =  g(ij) +  (1d0+k_love(l)) * ampl * leg(l) 
   		v(ij)  =  v(ij) +        l_love(l) * ampl * dleg
!   
     enddo
!   
     u(ij) = u(ij) * (3d0*rhoice/rhoear) * 1D0 * from_m_to_mm
     v(ij) = v(ij) * (3d0*rhoice/rhoear) * 1D0 * from_m_to_mm 
     g(ij) = g(ij) * (3d0*rhoice/rhoear) * 1D0 * from_m_to_mm

!
     end do
!$OMP END PARALLEL DO
!     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
     open(10,file=file_gf, status='unknown', recl=132) 
!
!
    IF(GRID_OPT==1)    THEN  
!
   Write(10,*) trim(adjustl(" - Love numbers from file: ")), trim(file_love_numbers)
   Write(10,*) trim(adjustl(" - Green functions for MAX degree: ")), LMAX 
   Write(10,*) trim(adjustl(" - Ice density (kg/m^3): ")), RHOICE
   Write(10,*) trim(adjustl(" - Earth density (kg/m^3): ")), RHOEAR
   Write(10,*) trim(adjustl(" - Load half-amplitude (deg): ")), ALFA
   Write(10,*) trim(adjustl(" - Ice thickness variation: +1 m/year"))
   Write(10,*) trim(adjustl(" - ======= Grid properties (the grid has costant spacing) ======="))
   Write(10,*) trim(adjustl(" - Min and Max colatitudes (deg): ")), THETA_MIN, THETA_MAX
   Write(10,*) trim(adjustl(" - Colatitude increment(deg): ")), THETA_INC
   Write(10,*) trim(adjustl(" - Number of grid points: ")), INT((THETA_MAX-THETA_MIN)/THETA_INC)+1
   Write(10,'(a1)') "#"
   Write(10,'(a1, a92)') "#", "   Line        colat, deg        vert. vel., mm/yr    hor. vel., mm/yr     geoid vel., mm/yr"
   Write(10,'(a1)') "#"

!
   ELSEIF(GRID_OPT==2) THEN
!
   Write(10,*) trim(adjustl(" - Love numbers from file: ")), trim(adjustl(file_love_numbers))
   Write(10,*) trim(adjustl(" - Green functions for MAX degree: ")), LMAX 
   Write(10,*) trim(adjustl(" - Ice density (kg/m^3): ")), RHOICE
   Write(10,*) trim(adjustl(" - Earth density (kg/m^3): ")), RHOEAR
   Write(10,*) trim(adjustl(" - Load half-amplitude (deg): ")), ALFA
   Write(10,*) trim(adjustl(" - Ice thickness variation: +1 m/year"))
   Write(10,*) trim(adjustl(" - ======= Grid properties (the grid has an increasing spacing) ======="))
   Write(10,*) trim(adjustl(" - Min and Max colatitudes (deg): ")), THETA_MIN, THETA_MAX
   Write(10,*) trim(adjustl(" - Minimum spacing (deg): ")), SPAC_MIN
   Write(10,*) trim(adjustl(" - Number of grid points: ")), NGRID
   Write(10,'(a1)') "#"
   Write(10,'(a1, a92)') "#", "   Line        colat, deg        vert. vel., mm/yr    hor. vel., mm/yr     geoid vel., mm/yr"
   Write(10,'(a1)') "#"
   ENDIF
!
     do ij=1,ngr
        write(10,'(i7,1x,e20.8,3(1x,e20.8))') ij, theta(ij), u(ij), v(ij), g(ij) 
     end do
     close(10)
!
!   
!
   Write(*,*) "make_gf: The Green functions are written on file: ", file_gf  
!
!
   deallocate( theta, u, v, g )
!
   STOP
!
   END PROGRAM GREEN_FUNCTION
!
!
!
!
!
