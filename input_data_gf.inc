!
! ###########################################   
! Input parameters for program "make_gf.f90"
! ###########################################   
!
!  Copyright Giorgio Spada & Daniele Melini, 2014
!
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
!
 real*8, parameter :: PI=3.14159265358979323840d0  ! Pi 
! 
!
! ##### Earth model parameters      
 real*8, parameter ::                   & 
         ggg=6.67384d-11,                &	  ! Newton's constant (SI units) 
         radius=6371d0,                   &    ! Radius of the Earth (km)
         radiusm=radius*1d3,               &    ! Radius of the Earth (m)
         grav=9.8046961d+00,                &	 ! Surface gravity (m/s/s) 
         rhoear=3d0*grav/4d0/ggg/pi/radiusm       ! Average Earth density (kg/m^3) 

! ##### Disc load parameters
! NOTE: The thickness variation is of +1 m/yr regardless of the load size 
 real*8, parameter ::            & 
         rhoice = 917d0,   	  &	  ! Ice density (kg/m^3) 
         alfa=0.02537d0           ! Half-amplitude of the load (deg)		     
!
!
!
!
! #### 1D Grid (two options available)
!
 integer, parameter :: grid_opt = 1    ! type of grid (1/2)
!
 real*8, parameter :: theta_min = 0d0     ! Min.colatitude (deg) both grids
 real*8, parameter :: theta_max = 60d0    ! Max colatitude (deg) both grids 
!
! ##### grid_opt=1: constant stepsize 
 real*8, parameter ::  theta_inc = alfa/5d0 	    ! Increment (deg)
!
! ##### grid_opt=2: stepsize increases with colatitude 
 integer, parameter :: ngrid=  1001	            ! number of points 
 real*8, parameter ::   spac_min=    alfa/5d0       ! min. grid spacing (deg)  
!
! ##### Max harmonic degrees for synthesis 
 integer, parameter ::     & 
          lmin=1,         &          ! lmin
	  lmax=32768                   ! lmax 
!
! ##### File name for Load-deformation coefficients (input) 
 character*200, parameter :: & 
                file_love_numbers="./DATA/REF_6371_loading_love_numbers.txt"   
  integer, parameter :: nh=1 ! Number of header lines in file_love_numbers
!
! ##### File name for gridded Green Function (output) 
 character*200, parameter :: file_gf="gf_Greenland.dat" 
!
