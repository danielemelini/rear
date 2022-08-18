!
!
!
subroutine PLegendre(p, lmax, z)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This subroutine evalutates all of the unnormalized Legendre polynomials 
!	up to degree lmax. 
!
!	Calling Parameters:
!		Out
!			p:	A vector of all unnormalized Legendgre polynomials evaluated at 
!				z up to lmax. The lenght must by greater or equal to (lmax+1).
!		IN
!			lmax:	Maximum degree to compute.
!			z:	Value within [-1, 1], cos(colatitude) or sin(latitude).
!
!	Notes:
!	
!	1.	The integral of Pl**2 over (-1,1) is 2/(2l+1).
!	2.	Values are calculated accoring to the following recursion scheme:
!			P_0(z) = 1.0, P_1(z) = z, and 
!			P_l(z) = (2l-1) * z * P_{l-1}(z) / l - (l-1) * P_{l-2}(z) / l
!
!	Dependencies:	None
!
!	Written by Mark Wieczorek June 2004
!
! ----> Modified to SINGLE PRECISION by Giorgio Spada 2007 
! ----> Also modified for the management of ERROR conditions 
!
!	Copyright (c) 2005, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
	implicit none
	integer, intent(in) ::	lmax
	real*8 paux(lmax+1)	
	real*8, intent(out) ::	p(0:lmax+1)		
       	real*8, intent(in) ::	z
       	real*8 :: pm2, pm1, pl
      	integer :: l, j

	if(size(p) < lmax+1) then
 	        Write(* ,*) "Error --- PlegendreL" 
		Write(* ,*) "P must be dimensioned as (LMAX+1) where LMAX is ", lmax
		Write(* ,*) "Input array is dimensioned ", size(p)
		Write(* ,*) "The program will STOP ----------------"	   
	        Stop
     	elseif (lmax < 0) then 
 	        Write(* ,*) "Error --- PlegendreL" 
		Write(* ,*) "LMAX must be greater than or equal to 0."
		Write(* ,*) "Input value is ", lmax
		Write(* ,*) "The program will STOP ----------------"	   
     	elseif(abs(z) > 1.) then
 	        Write(* ,*) "Error --- PlegendreL" 
		Write(* ,*) "ABS(Z) must be less than or equal to 1."
		Write(* ,*) "Input value is ", z
		Write(* ,*) "The program will STOP ----------------"	   
	        Stop
     	endif
      	
   	pm2  = 1.
      	paux(1) = 1.
      	
      	pm1  = z
      	paux(2) = pm1
      	
      	do l = 2, lmax
         	pl = (float(2*l-1) * z * pm1 - float(l-1) * pm2)/float(l)
         	paux(l+1) = pl
         	pm2  = pm1
         	pm1  = pl
      	enddo
!
        do j=0, lmax
             p(j)=paux(j+1)
	enddo
!
end subroutine PLegendre
!
!
!
!
!
!
