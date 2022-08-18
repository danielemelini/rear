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
! Configuration is read from an input file - DM Mar 10, 2022
! Changed status of input files from UNKNOWN to OLD - DM Apr 30, 2022
! Updated DM May 1, 2022 - added auto-detection of header lines in 
!    input files (and removed the NH options in the config files)
! Updated DM May 5, 2022: Added the option for uncompensated loads
!    and the 'asymptotic extension' of LNs
! Updated DM Aug 11, 2022: Fixed a bug in the shape factors of the
!    uncompensated disc load 
! Updated DM Aug 12, 2022: Re-written the LDC input block and fixed a 
!    bug in the asymptotic expansion of the LNs
! Updated DM Aug 14, 2022: Added a check for the ROT by Bevis et al. (2016)
!
! ------------------------------------------------------------------------  
!
! Declarations 
!
 PROGRAM GREEN_FUNCTION
 implicit NONE
!
 real*8, parameter :: PI=3.14159265358979323840d0  ! Pi 
! 
!
! ##### Earth model parameters       
 real*8, parameter :: ggg=6.67384d-11                	   ! Newton's constant (SI units) 
 real*8, parameter :: radius=6371d0                        ! Radius of the Earth (km)
 real*8, parameter :: radiusm=radius*1d3                   ! Radius of the Earth (m)
 real*8, parameter :: grav=9.8046961d+00                   ! Surface gravity (m/s/s) 
 real*8, parameter :: rhoear=3d0*grav/4d0/ggg/pi/radiusm   ! Average Earth density (kg/m^3) 
!
! +++++++++++++++++++++++ CONFIGURATION PARAMETERS +++++++++++++++++++++++
!
! ##### Disc load parameters
! NOTE: The thickness variation is of +1 m/yr regardless of the load size 
!
 real*8 :: alfa                   ! Half-amplitude of the load (deg)		     
 real*8 :: rhoice                 ! Ice density (kg/m^3) 
 integer :: icomp                 ! Uncompensated (0) or compensated (1) load
!
!
! #### 1D Grid
!
 integer :: grid_opt              ! type of grid (1/2)
!
 real*8 :: theta_min              ! Min.colatitude (deg) both grids
 real*8 :: theta_max              ! Max colatitude (deg) both grids 
!
 real*8 ::  theta_inc      	      ! Grid increment (deg)  [GRIDTYPE #1]
!
!
 integer :: ngrid                 ! number of points         [GRIDTYPE #2]
 real*8 ::  spac_min              ! min. grid spacing (deg)  [GRIDTYPE #2]
!
! ##### Min/max harmonic degrees for synthesis 
 integer :: lmin, lmax
!
! ##### Asymptotic extension of LNs
 integer :: iasympt
 integer :: n_asympt
 integer :: lmax_asympt
 real*8  :: h_inf, l_inf, k_inf
 integer :: nn
!
! ##### File name for Load-deformation coefficients (input) 
 character*200 :: file_love_numbers   
 integer ::       nh                   ! Number of header lines in file_love_numbers
!
! ##### File name for gridded Green Function (output) 
 character*200 :: file_gf
!
! +++++++++++++++++++++++ END OF CONFIGURATION PARAMETERS +++++++++++++++++++++++

! 
! --- Legendre Pols
 real*8, allocatable :: leg(:)         ! Legendre polynomial  [0:lmax+1]
 real*8 dleg                           ! Derivative of leg wrt colatitude
!
! --- Load function 
 real*8, allocatable :: sigma (:)      ! Spectral amplitude of the load       [0:lmax]
!
! --- Love numbers 
 real*8, allocatable :: h_love(:)      ! "h" Load deformation coefficient     [0:lmax]
 real*8, allocatable :: l_love(:)      ! "l" Load deformation coefficient     [0:lmax]
 real*8, allocatable :: k_love(:)      ! "k" Load deformation coefficient     [0:lmax]
 real*8 ll, hh, kk
!
! --- Miscellanea variables and constants  
 character*10 header
 integer i, l, ij 
 real*8, parameter :: from_m_to_mm=1D3 
 real*8  cosdd, ampl, r, q, x
 character*100 :: cfg_f
 character*200 :: buffer
 logical :: flag
 integer :: ieof
 integer :: n_in
!
! --- Dynamic arrays
 real*8, allocatable :: theta(:)
 real*8, allocatable :: u(:), v(:), g(:)
 integer :: ngr
! 
! ----------------------------------------------------------  
!
!
 write(*,*) ""
 write(*,*) "<<<<<<<<  REAR: a Regional ElAstic Rebound calculator  >>>>>>>>"
 write(*,*) ""
!
 call getarg(1,cfg_f)
!
 inquire(file=cfg_f,exist=flag)
!
 if( .not. flag ) then
    write(*,*) "make_gf: Cannot find configuration file: '", trim(cfg_f), "'" 
	stop
 end if
! 
 open(99,file=trim(cfg_f),status='old')
!
 write(*,*) 'make_gf: Using configuration file: ', trim(cfg_f)
!
 call read_data_line(99,buffer)   ;   read(buffer,*) rhoice
 call read_data_line(99,buffer)   ;   read(buffer,*) alfa
 call read_data_line(99,buffer)   ;   read(buffer,*) icomp
 call read_data_line(99,buffer)   ;   read(buffer,*) grid_opt
 call read_data_line(99,buffer)   ;   read(buffer,*) theta_min, theta_max
 call read_data_line(99,buffer)   ;   read(buffer,*) theta_inc, ngrid 
 call read_data_line(99,buffer)   ;   read(buffer,*) lmin, lmax
 call read_data_line(99,buffer)   ;   read(buffer,*) iasympt
 call read_data_line(99,buffer)   ;   read(buffer,*) n_asympt
 call read_data_line(99,buffer)   ;   read(buffer,*) lmax_asympt
 call read_data_line(99,buffer)   ;   file_love_numbers = trim(adjustl(buffer))
!call read_data_line(99,buffer)   ;   read(buffer,*) nh
 call read_data_line(99,buffer)   ;   file_gf = trim(adjustl(buffer))
!
 close(99)
!
!
 if( theta_inc<0 ) then
     theta_inc = alfa / abs(theta_inc)
 end if
 spac_min = theta_inc
!
 if( iasympt.eq.0 ) then
    allocate( leg(0:lmax+1) )
    allocate( sigma(0:lmax) )
    allocate( h_love(0:lmax) )
    allocate( l_love(0:lmax) )
    allocate( k_love(0:lmax) )
 else
    allocate( leg(0:lmax_asympt+1) )
    allocate( sigma(0:lmax_asympt) )
    allocate( h_love(0:lmax_asympt) )
    allocate( l_love(0:lmax_asympt) )
    allocate( k_love(0:lmax_asympt) ) 
 end if
!
 Write(*,*) "make_gf: Reading the Love numbers from file: ", trim(adjustl(file_love_numbers))
!
 open(1,file=file_love_numbers,status='old')
!
 call count_headers(1, nh)
 rewind(1) 
 do i=1, nh
    read(1,*)
 end do
!
 ieof=0
 n_in=0
 do while( ieof==0 )
      read( 1, *, iostat=ieof ) l, hh, ll, kk
      if( ieof==0 ) then
         if( (l.ge.lmin) .and. (l.le.lmax) ) then
             h_love(l) = hh
             l_love(l) = ll
             k_love(l) = kk
             n_in = n_in + 1
         end if
      end if
end do
!
 if( n_in .ne. (lmax - lmin + 1) ) then
    write(*,*) "make_gf: ERROR: the LDC database does not contains all the harmonic degrees in the requested range"
    stop
 end if
!
! do i=lmin, lmax
!       l=i 
!       read(1,*) ij, h_love(l), l_love(l), k_love(l)   
! enddo
!
 close(1)
!
!
 if (iasympt.eq.0) then
       Write(*,*) "make_gf: No asymptotic extension of LNs has been requested"
 else
       Write(*,*) "make_gf: Computing asymptotic extension for LNs from L=",lmax+1," to L=", lmax_asympt
!
       if( (lmax-n_asympt+1) .lt. lmin ) then 
			write(*,*) "make_gf: ERROR: number of harmonic degrees for the evaluation of asymptotic limits is too large"
			stop
       end if	   
!
       h_inf = 0.d0
	   l_inf = 0.d0
	   k_inf = 0.d0
!
	   do l=(lmax-n_asympt+1),lmax
	      h_inf = h_inf + h_love(l)
		  l_inf = l_inf + l_love(l) * l
		  k_inf = k_inf + k_love(l) * l
	   end do
!
	   h_inf = h_inf / dble(n_asympt)
	   l_inf = h_inf / dble(n_asympt)
	   k_inf = h_inf / dble(n_asympt)
!	   
       do l=lmax+1, lmax_asympt
	       h_love(l) = h_inf
		   l_love(l) = l_inf / dble(l)
		   k_love(l) = k_inf / dble(l)
	   end do
!
       lmax = lmax_asympt
!
 end if
!
!
!
! Check that alpha and lmax satisfy the "rule of thumb"
 if( lmax .lt. nint( 360/alfa ) ) then
    Write(*,*) "make_gf: WARNING: A disc of half-amplitude ", alfa, " deg requires at least lmax=",nint(360/alfa)
	Write(*,*) "                  *** Results may be inaccurate. Consider increasing lmax ***"
 end if
!
!
!
 Write(*,*) "make_gf: Computing the harmonic coefficients of the load ..."
!
   call  plegendre(leg, lmax+1, cosdd(alfa))
!
   if (icomp==0) then         ! Uncompensated load
!
      Write(*,*) "make_gf: The load is NOT compensated"
!
      do l=lmin, lmax
 	     if(l==0) sigma(l)= 1d0-cosdd(alfa)
 	     if(l>=1) sigma(l)= leg(l-1)-leg(l+1)
      end do
!
      sigma = sigma * 0.5d0
!
   elseif (icomp==1) then     ! Compensated load
!
      Write(*,*) "make_gf: The load is compensated"
!
      do l=lmin, lmax
 	     if(l==0) sigma(l)= 0D0
 	     if(l>=1) sigma(l)= -(leg(l+1)-leg(l-1))/(1d0+cosdd(alfa))
      end do
!
   else
!
      Write(*,*) "make_gf: ERROR: unknown compensation option"
	  stop
!
   end if
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
!$OMP    SHARED(NGR,THETA,SIGMA,U,V,G,H_LOVE,K_LOVE,L_LOVE,LMIN,LMAX) &
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
!
     end do
!$OMP END PARALLEL DO
!
 u = u * (3d0*rhoice/rhoear) * 1D0 * from_m_to_mm
 v = v * (3d0*rhoice/rhoear) * 1D0 * from_m_to_mm 
 g = g * (3d0*rhoice/rhoear) * 1D0 * from_m_to_mm
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
   deallocate( leg )
   deallocate( sigma )
   deallocate( h_love, l_love, k_love )
!
   STOP
!
   END PROGRAM GREEN_FUNCTION
!
!
!
!
!
