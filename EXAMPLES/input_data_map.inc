!
! ###########################################   
! Input parameters for program "make_map.f90"
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
! >>>> INPUT PARAMETERS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!
 real*8, parameter :: RADIUS = 6371d0  ! Earth radius (km) 
 real*8, parameter :: PI=3.14159265358979323840d0 
 real*8, parameter :: D2R=PI/180.D0
 real*8, parameter :: R2D=180.D0/PI
! 
! >>>> ICE MODEL  
!
 character*200, parameter :: file_ice="./DATA/greeM3R.dat"   ! File of the ice model  
 integer, parameter :: NH_ICE=10                             ! Header lines in file_ice
!
! >>>> GREEN'S FUNCTION  
!
 character*200, parameter :: file_gf="gf_Greenland.dat"     ! gridded Green function 
 integer, parameter :: GRID_OPT=1                           ! grid style (1=uniform, 2=increasingly sparse)
 integer, parameter :: NH_GF=13                             ! Header lines in file_gf

!
!
!
! >>>> SELECTION of the OUTPUT <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!
! IWHERE=1 ---> Rates of displacement are on a 2D grid
! IWHERE=2 ---> Rates of displacement on isolated points
!
  integer, parameter :: iwhere = 1
!
! -------------------------------------------------------------
! >>>> PARAMETERS for IWHERE = 1 (2D grid) <<<<<<<<<<<<<<<<<<<<
! -------------------------------------------------------------
!
 REAL*8, PARAMETER :: GRD_TYPE=2             ! Grid type (1=Tegmark, 2=lon/lat)       
!
 REAL*8, PARAMETER :: PIXEL_AREA=100d0       ! Pixel area for GRD_TYPE=1 (km^2)    
 REAL*8, PARAMETER :: GRD_SPACING=0.5D0      ! Grid spacing for GRD_TYPE=2 (deg) 
!
 REAL*4, PARAMETER :: LONMIN_GR=285.0, LONMAX_GR=350.0   ! Region boundaries (lon) (for both types)
 REAL*4, PARAMETER :: LATMIN_GR= 58.0, LATMAX_GR= 86.0   ! Region boundaries (lat) 
! 
 CHARACTER*200, PARAMETER :: FILE_OUT_GRID="uvg.dat"    ! Output filename (for both types)
!
! -------------------------------------------------------------
! >>>> PARAMETERS for IWHERE = 2 (isolated points) <<<<<<<<<<<<
! -------------------------------------------------------------
!
 CHARACTER*200, PARAMETER :: FILE_IN_POINTS="GPS-points.dat"  ! Input Filename  (points)
 integer, parameter :: NH_IP=6                               ! Header lines 
!
 CHARACTER*200, PARAMETER :: FILE_OUT_POINTS="GPS-out.dat"    ! Output filename (points) 
!
! END of FILE
!
!
!
!
