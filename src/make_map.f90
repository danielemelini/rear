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
! Computes the response to a surface load on a 2D grid or on sparse points
!                    
! History: 
! First written by DM for the GJI paper on Greenland 
! Retouched by GS summer 2013 for layout ETC 
! A bug fixed by DM on Aug 7, 2013 
! Updated several times until July 2013 
! Updated GS August 2013 for the "Wang" 1D grid implementation 
! Updated DM Aug 18, 2013 with the optimized Wang 1D grid
! Updated DM Aug 19, 2013 -- NEW ice model format with height rate (m/yr)
! Seen by Giorgio Aug 20: Thanks Daniele! 
! Updated DM April 28, 2014 (OpenMP)
! Configuration is read from an input file - DM Mar 10, 2022
! Added uncertainties - DM Mar 22, 2022
! Updated DM Apr 27, 2022 - added some checks to avoid numerical
!    issues when the observer is on a disc center. Parameter EPS
!    controls the tolerance. If cos(distance)>(1-EPS), observer
!    is considered at the disc center. EPS=1e-11 as a default
! Added an option to include/not include uncertainties - DM Apr 30 2022
!    also changed status of input files from UNKNOWN to OLD
! Updated DM May 1, 2022 - added auto-detection of header lines in
!    input files (and removed the NH options in the config files)
!
! ------------------------------------------------------------------------  
! 
!
! 
 PROGRAM MAKE_MAP  
 IMPLICIT NONE
! 
!
! ///// Declarations
!
! --- Miscellanea 
 INTEGER I, J, K, IJ
 CHARACTER*30 JUNK
 CHARACTER*200 BUFFER  
 CHARACTER*100 CFG_F
 LOGICAL :: FLAG
!
!
 real*8, parameter :: RADIUS = 6371d0  ! Earth radius (km) 
 real*8, parameter :: PI=3.14159265358979323840d0 
 real*8, parameter :: D2R=PI/180.D0
 real*8, parameter :: R2D=180.D0/PI
 real*8, parameter :: EPS=1.D-11
! 
!
! ///// Configuration parameters
!
!
! 
! >>>> ICE MODEL  
!
 character*200 :: file_ice   ! File of the ice model  
 integer :: NH_ICE           ! Header lines in file_ice
 integer :: ERR_OPT          ! Uncertainties in ice model (0=no, 1=yes)
 logical :: ERR_FLAG
!
! >>>> GREEN'S FUNCTION  
!
 character*200 :: file_gf         ! gridded Green function 
 integer :: GRID_OPT              ! grid style (1=uniform, 2=increasingly sparse)
 integer, parameter :: NH_GF=13   ! Header lines in file_gf
!
! >>>> SELECTION of the OUTPUT <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!
 integer :: iwhere
!
! -------------------------------------------------------------
! >>>> PARAMETERS for IWHERE = 1 (2D grid) <<<<<<<<<<<<<<<<<<<<
! -------------------------------------------------------------
!
 INTEGER :: GRD_TYPE=2            ! Grid type (1=Tegmark, 2=lon/lat)       
 REAL*8 :: PIXEL_AREA             ! Pixel area for GRD_TYPE=1 (km^2)    
 REAL*8 :: GRD_SPACING            ! Grid spacing for GRD_TYPE=2 (deg) 
 REAL*4 :: LONMIN_GR, LONMAX_GR   ! Region boundaries (lon) (for both types)
 REAL*4 :: LATMIN_GR, LATMAX_GR   ! Region boundaries (lat) 
!
 CHARACTER*200 :: FILE_OUT_GRID   ! Output filename (for both types)
!
! -------------------------------------------------------------
! >>>> PARAMETERS for IWHERE = 2 (isolated points) <<<<<<<<<<<<
! -------------------------------------------------------------
!
 CHARACTER*200 :: FILE_IN_POINTS  ! Input Filename  (points)
 integer :: NH_IP                 ! Header lines 
 CHARACTER*200 :: FILE_OUT_POINTS ! Output filename (points) 
!
!
! ///// Green Function grid 
!
 INTEGER NGRID
 REAL*8 DELTA 
 REAL*8 COSTHETA, COSGAMMA, SINGAMMA
 REAL*8 THETA, THETA_STEP
 REAL*8 COSDD, SINDD 
 REAL*8, ALLOCATABLE :: THETA_GRID(:)
 REAL*8, ALLOCATABLE :: UR_GRID(:), & 
                        UH_GRID(:), &
                        N_GRID(:)
 REAL*8 T1(1), T2(1)  
 INTEGER NSEC(1)  
 REAL*8 XI, R, Q, SPAC_MIN
!
!
! ///// Observer (where the solution is computed) 
!
 INTEGER NP, NPP, IN_GR, RES_OBS, NLON, NLAT
 REAL*8 GRD_STEP
 LOGICAL COND_LAT, COND_LON  
 REAL*4, ALLOCATABLE :: LONP(:), LATP(:)
 INTEGER, ALLOCATABLE :: IDX(:)
 REAL*8 UR_D, UTH_D, UPHI_D, N_D, UH_D
 REAL*8, ALLOCATABLE :: U_OBS(:),    & 
                        UTH_OBS(:),  & 
			UPHI_OBS(:), & 
			N_OBS(:)
 REAL*8, ALLOCATABLE :: SIG_U_OBS(:),    & 
                        SIG_UTH_OBS(:),  & 
			SIG_UPHI_OBS(:), & 
			SIG_N_OBS(:)
 REAL*8, ALLOCATABLE :: LON_OBS(:), LAT_OBS(:) 
 REAL*8, ALLOCATABLE :: COS_LON_OBS(:), SIN_LON_OBS(:)
 REAL*8, ALLOCATABLE :: COS_LAT_OBS(:), SIN_LAT_OBS(:)
 CHARACTER*60 NAME_SITE
 CHARACTER*4, allocatable :: NAME_STAT(:)
!
!
! ///// Ice model information
!
 INTEGER NICE 
 REAL*8  ALFA_GF, ALFA
 REAL*8, PARAMETER :: ALFA_EPS = 1.D-8
 REAL*8, ALLOCATABLE :: ALFA_ICE(:)
 REAL*8, ALLOCATABLE :: LON_ICE(:)
 REAL*8, ALLOCATABLE :: LAT_ICE(:)
 REAL*8, ALLOCATABLE :: H1(:)
 REAL*8, ALLOCATABLE :: SIG_H1(:)
 REAL*8, ALLOCATABLE :: COS_LON_ICE(:)
 REAL*8, ALLOCATABLE :: SIN_LON_ICE(:)
 REAL*8, ALLOCATABLE :: COS_LAT_ICE(:)
 REAL*8, ALLOCATABLE :: SIN_LAT_ICE(:)
!
!
! 
! !!!!!!!!!!!!!!!!!!!!!!!!!!!! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Execution starts here     !    Execution starts here     !    Execution starts here
! !!!!!!!!!!!!!!!!!!!!!!!!!!!! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
    write(*,*) "make_map: Cannot find configuration file: '", trim(cfg_f), "'" 
	stop
 end if
! 
 open(99,file=trim(cfg_f),status='old')
!
 write(*,*) 'make_map: Using configuration file: ', trim(cfg_f)
!
 call read_data_line(99,buffer)   ;   file_ice = trim(adjustl(buffer))
 call read_data_line(99,buffer)   ;   read(buffer,*) err_opt 
 call read_data_line(99,buffer)   ;   file_gf = trim(adjustl(buffer))
 call read_data_line(99,buffer)   ;   read(buffer,*) grid_opt
 call read_data_line(99,buffer)   ;   read(buffer,*) iwhere
 call read_data_line(99,buffer)   ;   read(buffer,*) grd_type
 call read_data_line(99,buffer)   ;   read(buffer,*) pixel_area
 call read_data_line(99,buffer)   ;   read(buffer,*) grd_spacing
 call read_data_line(99,buffer)   ;   read(buffer,*) lonmin_gr, lonmax_gr
 call read_data_line(99,buffer)   ;   read(buffer,*) latmin_gr, latmax_gr
 call read_data_line(99,buffer)   ;   file_out_grid = trim(adjustl(buffer))
 call read_data_line(99,buffer)   ;   file_in_points = trim(adjustl(buffer))
 call read_data_line(99,buffer)   ;   file_out_points = trim(adjustl(buffer))
!
 close(99)
!
!
!
! -------------------------------------------------------------
 Write(*,*) "make_map: Reading ice model information "
! ------------------------------------------------------------- 
!
 if( err_opt.eq.0 ) then
    write(*,*) "make_map: Assuming NO uncertainties on ice model"
    err_flag=.false.
 elseif( err_opt.eq.1 ) then
    write(*,*) "make_map: Expecting uncertainties in the ice model file"
    err_flag=.true.
 else
    write(*,*) "make_map: ERROR: unknown option for ice model uncertainties"
    stop
 end if
!
 open(20,file=file_ice,status='old')
!
! Reading the number of lines in the ice model file
!
 call count_rows(20,nice)
 rewind(20)
 call count_headers(20,nh_ice)
 rewind(20)
 nice = nice-NH_ICE
 write(*,*) "make_map: Found ", nice, " ice elements in file ", trim(adjustl(file_ice))
!
 allocate( lon_ice(nice), lat_ice(nice), h1(nice), sig_h1(nice), alfa_ice(nice) )  
!
! Reading the header of the ice model file
 do i=1, NH_ICE 
	read(20,'(a30)')junk 
 enddo
!
! Reading the the ice model information
 if( err_flag ) then
    do i=1, nice  
       read(20,*) k, lon_ice(i), lat_ice(i), alfa_ice(i), h1(i), sig_h1(i) 
    end do
 else
    do i=1, nice  
       read(20,*) k, lon_ice(i), lat_ice(i), alfa_ice(i), h1(i)
    end do
    sig_h1 = 0.d0
 end if
! 
 close(20) 
!
!
! 
 if( abs(maxval(alfa_ice) - minval(alfa_ice)) .gt. alfa_eps ) then
    write(*,*) "make_map: Error: ice elements must have constant size."
    stop
 else
    alfa = alfa_ice(1)
 endif    
!
!
!
!--------------------------------------------------------------------
 Write(*,*) "make_map: Reading the Green Function file: ", trim(adjustl(file_gf))
!--------------------------------------------------------------------
!
 Open(30,file=file_gf,status='old') 
!
! Reading the disc half-amplitude from the header lines
!
 do i=1,4 ; read(30,*) ; end do
 read(30,'(a)')  buffer
 read(buffer(31:sizeof(buffer)),*) alfa_gf
 if( abs((alfa-alfa_gf)) .gt. alfa_eps ) then
    write(*,*) "make_map: Error: Disc element size mismatch between ice model and Green Function files"
    write(*,*) "make_map: Disc half-amplitude in ice model file is ", alfa 
    write(*,*) "make_map: Disc half-amplitude used for Green Functions is ", alfa_gf
    write(*,*) "make_map: The program will stop."
    stop
 end if
 rewind(30)
!
! Reading the number of lines
!
 call count_rows(30,ngrid)
 rewind(30)
 ngrid = ngrid-NH_GF
 write(*,*) "make_map: Reading ", ngrid, " points from the Green Function file " 
!
 allocate( theta_grid(ngrid), ur_grid(ngrid), uh_grid(ngrid), n_grid(ngrid) )
!
! Reading the header of the GF file
 do i=1, NH_GF 
	read(30,'(a30)') JUNK 
 enddo 
!
! Reading the GF information 
 do i=1, ngrid
 	read(30,*) k, theta_grid(i), ur_grid(i), uh_grid(i), n_grid(i)
 enddo
!
 close(30) 
!
 NSEC(1)=NGRID 
 T1(1)=theta_grid(1)
 T2(1)=theta_grid(ngrid)
! 
 write(*,*) "make_map: Grid data: min colat = ", T1
 write(*,*) "make_map: Grid data: max colat = ", T2
 write(*,*) "make_map: Number of points: = ", NGRID
! 
!
!
 IF(IWHERE == 1) THEN    
! --------------------------------------------------------------------------
   write(*,*) "make_map: Option 1: Gridded output "   
! --------------------------------------------------------------------------
!  
 IF( GRD_TYPE == 1) THEN
! --------------------------------------------------------------------------
   write(*,*) "make_map: Using the Tegmark gridding"   
! --------------------------------------------------------------------------
!
! New as of Aug 8, 2013 (GS)
! ---------------------------
! NP is ~ the ratio (area of the surface of the sphere)/(pixel area):  
  NPP=INT(4D0*PI*RADIUS*RADIUS/PIXEL_AREA)+1
!
! The Tegmark resolution is given by the inverse formula (see my notes): 
  RES_OBS=INT(0.5D0 + SQRT(10d0*NPP-20)/20D0)+1
!
! From the Tegmark formula, the effective number of pixels is:
  NP=2*RES_OBS*(RES_OBS-1)*20+12
!
   ALLOCATE( LONP(NP) )
   ALLOCATE( LATP(NP) )
   ALLOCATE( IDX(NP) )
!
! Finding the pixels coordinates globally 
   call FINDPX(res_obs, np, lonp, latp)
!
   write(*,*) "make_map: User supplied pixel area (km^2): ", pixel_area
   write(*,*) "make_map: Approximate global number of pixels: ", npp
   write(*,*) "make_map: Tegmark resolution: ", RES_OBS
   write(*,*) "make_map: Effective global number of pixels: ", np
   write(*,*) "make_map: Approximate spacing of the GRID (km): ", 2D0*SQRT(PIXEL_AREA/PI)
!
   in_gr=0
!
   do i=1, NP 
!
        cond_lat=.FALSE.
        cond_lon=.FALSE.
!
 	    cond_lat=( (latp(i).ge.latmin_gr).and.(latp(i).le.latmax_gr) )
            cond_lon=( (lonp(i).ge.lonmin_gr).and.(lonp(i).le.lonmax_gr) )
!
	    if(cond_lat.and.cond_lon) THEN   
		  in_gr= in_gr+1
                  idx(in_gr) = i
	    endif  
   enddo
!
   allocate( lon_obs(in_gr) )
   allocate( lat_obs(in_gr) )
!
   lon_obs = lonp( idx(1:in_gr) )
   lat_obs = latp( idx(1:in_gr) )
!
   deallocate( lonp )
   deallocate( latp )
   deallocate( idx )
!
 ELSEIF( GRD_TYPE == 2 ) THEN
! -------------------------------------------------------------
   write(*,*) "make_map: Using a uniform lon/lat grid"   
! -------------------------------------------------------------
!
! Revised as of Aug 8, 2013 (GS)
! ------------------------------
   NLON = NINT( (LONMAX_GR - LONMIN_GR) / GRD_SPACING ) + 1
   NLAT = NINT( (LATMAX_GR - LATMIN_GR) / GRD_SPACING ) + 1
   NP=NLON*NLAT
!   
   allocate( lon_obs(NP) )
   allocate( lat_obs(NP) )
!
   write(*,*) "make_map: Spacing of the 2D GRID of observers (deg): ", GRD_SPACING
   write(*,*) "make_map: Number of GRID points within the region: ", np
!
   in_gr=0
!
   do i=1, nlat
      do j=1, nlon
!
        in_gr = in_gr + 1
!
        lon_obs(in_gr) = lonmin_gr + (j-1) * GRD_SPACING
        lat_obs(in_gr) = latmin_gr + (i-1) * GRD_SPACING
!
	enddo
! 
  enddo
!
  in_gr = np
!
 ENDIF  ! On the type of grid 
!
 Write(*,*) "make_map: Number of pixels (observers) within the region: ", in_gr
!
!
!

 ELSEIF(IWHERE==2) THEN 
! --------------------------------------------------------------------------
   write(*,*) "make_map: Option 2: Sparse points "   
! --------------------------------------------------------------------------
!
 open(30,file=file_in_points,status='old') 
!
!
! Finding the number of observer points
 call count_rows(30,in_gr)
 rewind(30)
 call count_headers(30,nh_ip)
 rewind(30)
 in_gr = in_gr-NH_IP
 write(*,*) "make_map: Reading ", in_gr, " sparse points from file ", trim(adjustl(file_in_points))
!
 ALLOCATE( lon_obs(in_gr) )
 ALLOCATE( lat_obs(in_gr) )
 allocate( name_stat(in_gr) )
!
!
! Reading the header of the sparse points file
 do i=1, nh_ip
	read(30,'(a30)')junk  
 enddo
!
! Reading the sparse points information
 do i=1, in_gr 
	read(30,*) k, lon_obs(i), lat_obs(i), name_site
	NAME_STAT(I)=trim(adjustl(NAME_SITE))
	if(lon_obs(i).le.0d0) lon_obs(i)=360D0+lon_obs(i)
 enddo
!
 close(30) 
!
!
 ENDIF   ! on "IWHERE" (where the solutions are computed) 
!
!
 ALLOCATE( U_OBS(IN_GR) )
 ALLOCATE( UTH_OBS(IN_GR) )
 ALLOCATE( UPHI_OBS(IN_GR) )
 ALLOCATE( N_OBS(IN_GR) )
!
 ALLOCATE( SIG_U_OBS(IN_GR) )
 ALLOCATE( SIG_UTH_OBS(IN_GR) )
 ALLOCATE( SIG_UPHI_OBS(IN_GR) )
 ALLOCATE( SIG_N_OBS(IN_GR) )
!
!
 Write(*,*) "make_map: Pre-computing trigonometric functions"
!
 allocate( sin_lon_obs( in_gr ) )
 allocate( cos_lon_obs( in_gr ) )
 allocate( sin_lat_obs( in_gr ) )
 allocate( cos_lat_obs( in_gr ) )
!
 allocate( sin_lon_ice( nice ) )
 allocate( cos_lon_ice( nice ) )
 allocate( sin_lat_ice( nice ) )
 allocate( cos_lat_ice( nice ) )
!
 sin_lon_obs = sin( lon_obs * d2r )
 cos_lon_obs = cos( lon_obs * d2r )
 sin_lat_obs = sin( lat_obs * d2r )
 cos_lat_obs = cos( lat_obs * d2r )
!
 sin_lon_ice = sin( lon_ice(1:nice) * d2r )
 cos_lon_ice = cos( lon_ice(1:nice) * d2r )
 sin_lat_ice = sin( lat_ice(1:nice) * d2r )
 cos_lat_ice = cos( lat_ice(1:nice) * d2r )
!
! --- Computing some variables for the "increasingly sparse" grid
!
 if (grid_opt.eq.2) then
    spac_min = theta_grid(2) - theta_grid(1)
    R = 2d0*(t2(1)-t1(1))/float(ngrid-1)/spac_min - 1d0
    Q = spac_min
endif
!
 Write(*,*) "make_map: Computing displacements at the requested points"
!
 u_obs    = 0d0
 uth_obs  = 0d0
 uphi_obs = 0d0
 n_obs    = 0d0
! 
 sig_u_obs    = 0d0
 sig_uth_obs  = 0d0
 sig_uphi_obs = 0d0
 sig_n_obs    = 0d0
!
!$OMP PARALLEL DO &
!$OMP    DEFAULT(NONE) &
!$OMP    PRIVATE(I,J,COSTHETA,K,IJ,THETA,XI,DELTA,THETA_STEP,UR_D,UH_D,N_D,COSGAMMA,SINGAMMA,UTH_D,UPHI_D,FLAG) &
!$OMP   SHARED(IN_GR,NICE,SIN_LAT_ICE,COS_LAT_ICE,SIN_LON_ICE,COS_LON_ICE, &
!$OMP          SIN_LAT_OBS,SIN_LON_OBS,COS_LAT_OBS,COS_LON_OBS, &
!$OMP          T1,T2,NSEC,NGRID,THETA_GRID,UR_GRID,N_GRID,UH_GRID,Q,R,H1,SIG_H1, &
!$OMP          U_OBS,UTH_OBS,UPHI_OBS,N_OBS, &
!$OMP          SIG_U_OBS,SIG_UTH_OBS,SIG_UPHI_OBS,SIG_N_OBS,GRID_OPT,ERR_FLAG)
!
 loop_obs: do i = 1,in_gr
!
! if(mod(i,1000)==0)write(*,*) " Observer: ", i
!		
	loop_ice: do j = 1, nice 
!
!       Angular separation between the i-th observed and the j-th ice element
!  
        costheta = sin_lat_ice(j) * sin_lat_obs(i) + &
                   cos_lat_ice(j) * cos_lat_obs(i) * &
                   (cos_lon_ice(j) * cos_lon_obs(i) + sin_lon_ice(j) * sin_lon_obs(i))
!
        if (costheta > (1.d0 - eps)) then
           costheta = 0.d0
           theta = 0.d0
           flag  = .true.
        else
           theta = acos(costheta)*r2d
           flag  = .false.
        end if
!
! --- Finding the position of the observer in the 1D grid          
!
        if (grid_opt .eq. 1) then                   ! A uniform grid    
!	
            k = -1
!
            loop_grid: do ij=1,1
!
            if( (theta.ge.t1(ij)) .and. (theta.le.t2(ij)) ) then  
!
                 k = ceiling( (theta - t1(ij)) / (t2(ij)-t1(ij)) * (nsec(ij)-1) )
!
                 k = k + sum( nsec(1:(ij-1)) )
!
            end if
!
            end do loop_grid            
!
            if (k .eq. -1)  cycle  loop_ice
!
        elseif (grid_opt .eq. 2) then                   ! Variable-step grid
!	    
            delta = (theta - t1(1)) / q
!
            xi   = (1 - 2*ngrid +3*r + &
               sqrt( dble(9 + 16*delta -12*ngrid - 8*delta*ngrid + &
                     4*ngrid**2 + 6*r - 16*delta*r - 4*ngrid*r + &
                      8*delta*ngrid*r + r**2) ) ) / &
	              ( 2*(r-1) )
!		      
            k    = floor(xi)
!
            if ( (k.lt.1) .or. (k.gt.ngrid) )  cycle loop_ice
!
        end if
!               
    	theta_step= theta_grid(k+1) - theta_grid(k)			
!
! --- Vertical deformation ---
!
        delta = ur_grid(k+1)-ur_grid(k)
        ur_d = ur_grid(k) + ((theta - theta_grid(k))/theta_step)*delta  	
!       ur_d = ur_d * h1(j) 
!
! --- Geoid height ---
!
	delta = n_grid(k+1)-n_grid(k)
        n_d  = n_grid(k) + ((theta - theta_grid(k))/theta_step)*delta  	
!       n_d  = n_d * h1(j)
!
! --- Horizontal velocity ---
!
	delta = uh_grid(k+1)-uh_grid(k)
        uh_d = uh_grid(k) + ((theta - theta_grid(k))/theta_step)*delta 	
!       uh_d = uh_d* h1(j)
!               
        cosgamma = ( sin_lat_ice(j) - sin_lat_obs(i) * costheta ) / &
           ( cos_lat_obs(i) * sqrt( 1 - costheta ** 2 ) )
!
        singamma = ( cos_lon_ice(j) * sin_lon_obs(i) - cos_lon_obs(i) * sin_lon_ice(j) ) * &
           cos_lat_ice(j)   / sqrt( 1 - costheta ** 2 )  
!       
        if( flag ) then
           uphi_d = 0.d0
           uth_d  = 0.d0
        else
           uphi_d = uh_d * singamma
           uth_d  = uh_d * cosgamma
        end if
!
!
! --- Updating the vertical displacement at the observer 
!
        u_obs(i) =    u_obs(i)    + h1(j) * ur_d 
        uth_obs(i)  = uth_obs(i)  + h1(j) * uth_d 
        uphi_obs(i) = uphi_obs(i) + h1(j) * uphi_d 
        n_obs(i)    = n_obs(i)    + h1(j) * n_d 
!
        if( err_flag ) then
           sig_u_obs(i)    = sig_u_obs(i)    + ( sig_h1(j) * ur_d    )**2  
           sig_uth_obs(i)  = sig_uth_obs(i)  + ( sig_h1(j) * uth_d   )**2
           sig_uphi_obs(i) = sig_uphi_obs(i) + ( sig_h1(j) * uphi_d  )**2
           sig_n_obs(i)    = sig_n_obs(i)    + ( sig_h1(j) * n_d     )**2
        end if
!
     end do loop_ice 
! 
end do loop_obs
!
 Write(*,*) "make_map: END of computation"
!
 if( err_flag ) then
    sig_u_obs    = sqrt( sig_u_obs    )
    sig_uth_obs  = sqrt( sig_uth_obs  )
    sig_uphi_obs = sqrt( sig_uphi_obs )
    sig_n_obs    = sqrt( sig_n_obs    )
 end if
!
 Write(*,*) "make_map: Writing the output"
! 
! 
    if(IWHERE==1) Open(44,file=file_out_grid,  status='unknown')
    if(IWHERE==2) Open(44,file=file_out_points,status='unknown')
!
    if(IWHERE==1)  Write(44,*) "#", " Lon, Lat, UrDOT, UthetaDOT, UlambdaDOT, NDOT (mm/yr) "  
    if(IWHERE==2)  Write(44,*) "#", " Lon, Lat, Up, South, East, Geoid (mm/yr), Site"  
    Write(44,*) "#", " Green function from file: ", trim(adjustl(file_gf)) 
! 
    if( err_flag ) then
!            
    do i=1,in_gr
	IF    (IWHERE==1) THEN 
	   write(44,'(2(f10.4,2x),8(e12.6,2x))')  lon_obs(i),  &
	                                          lat_obs(i),  & 
						  u_obs(i),    sig_u_obs(i),    &
	                                          uth_obs(i),  sig_uth_obs(i),  & 
						  uphi_obs(i), sig_uphi_obs(i), & 
						  n_obs(i),    sig_n_obs(i)
        ELSEIF(IWHERE==2) THEN
	   write(44,'(2(f15.9,2x),8(e12.6,2x),A4)')  lon_obs(i),  &
	                                          lat_obs(i),  & 
						  u_obs(i),    sig_u_obs(i),    &
	                                          uth_obs(i),  sig_uth_obs(i),  & 
						  uphi_obs(i), sig_uphi_obs(i), & 
						  n_obs(i),    sig_n_obs(i),    &
						  NAME_STAT(i)  
	ENDIF
!
    enddo
!
!
    else
!
!
    do i=1,in_gr
	IF    (IWHERE==1) THEN 
	   write(44,'(2(f10.4,2x),4(e12.6,2x))')  lon_obs(i),  &
	                                          lat_obs(i),  & 
						  u_obs(i),    &
	                                          uth_obs(i),  & 
						  uphi_obs(i), & 
						  n_obs(i)
        ELSEIF(IWHERE==2) THEN
	   write(44,'(2(f15.9,2x),4(e12.6,2x),A4)')  lon_obs(i),  &
	                                          lat_obs(i),  & 
						  u_obs(i),    &
	                                          uth_obs(i),  & 
						  uphi_obs(i), & 
						  n_obs(i),    &
						  NAME_STAT(i)  
	ENDIF
!
    enddo
!
    end if
!
!
    close(44)
!
    if(IWHERE==1) Write(*,*) "make_map: The output is reported on file: ", trim(adjustl(file_out_grid))   
    if(IWHERE==2) Write(*,*) "make_map: The output is reported on file: ", trim(adjustl(file_out_points))  
!
 deallocate( lon_obs )
 deallocate( lat_obs )
 deallocate( u_obs )
 deallocate( uth_obs )
 deallocate( uphi_obs )
 deallocate( n_obs )
!
 deallocate( sig_u_obs )
 deallocate( sig_uth_obs )
 deallocate( sig_uphi_obs )
 deallocate( sig_n_obs )
!
 deallocate( h1 )
 deallocate( sig_h1 )
!
 deallocate( sin_lon_obs )
 deallocate( cos_lon_obs )
 deallocate( sin_lat_obs )
 deallocate( cos_lat_obs )
!
 deallocate( sin_lon_ice )
 deallocate( cos_lon_ice )
 deallocate( sin_lat_ice )
 deallocate( cos_lat_ice )
!
 end program MAKE_MAP
!
!
!
!
!
