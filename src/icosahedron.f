!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!       THE ICOSAHEDRON PACKAGE FOR PIXELIZING THE SPHERE        !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
!    Written by Max Tegmark, Max-Planck-Institut fuer Physik, Munich
!    April 1996
!    Currently I'm at MIT, tegmark@mit.edu
!    Bug in vector2pixel fixed by Robert Berry in May 13
!
!  WHAT IS IT?
!    This FORTRAN package lets the user pixelize the sphere at 
!    a wide range of resolutions. It was written primarily for 
!    map-making in astronomy and cosmology. It is also useful 
!    for doing integrals over the sphere when the integrand is 
!    expensive to evaluate, so that one wishes to minimize the 
!    number of points used.
!
!  DOCUMENTATION:
!    The package and its purpose is described in detail in 
!    a postscript file available from
!      http://www.mpa-garching.mpg.de/~max/icosahedron.html
!    (faster from Europe) and from from 
!      http://sns.ias.edu.edu/~max/icosahedron.html
!    (faster from the US). This site also contains the latest 
!    version of the source code.
!
!  RULES:
!    The package is public domain, which means that you are
!    allowed to use it for any non-commercial purpose whatsoever 
!    free of charge. The only requirement is that you include an 
!    appropriate acknowledgement in any publications based on
!    work where the package has been used. Also, if you 
!    redistribute the package, you must not remove this text.
!
!  HOW IT WORKS:
!    As a supplement to the above-mentioned postscript file, 
!    here is a brief summary of the nitty-gritty details. 
!    To use the package, first call the subroutines 
!    compute_matrices and compute_corners once and for all, 
!    as in the demo routine below. This precomputes some 
!    geometric stuff to save time later.
!    The RESOLUTION is a positive integer 1, 2, 3, ... that you
!    can choose freely. It determines the number of pixels, which
!    is given by 
!       N = 40*resolution*(resolution-1)+12.
!    For instance, resolution=13 gives N=6252 pixels.
!    The only subroutines you ever need to use are these two:
!    * vector2pixel takes a point on the sphere (specified by a 
!      unit vector) and returns the number of the pixel to which
!      this point belongs, i.e., an integer between 0 and N-1.
!    * pixel2vector does the opposite: it takes a pixel number and
!      computes the corresponding unit vector.
!    The subroutine "demo" below illustrates how the routines are used.
!    It produces a text file with the coordinates of all pixels
!    together with their pixel numbers. It also outputs the pixel 
!    number reconstructed with vector2pixel, so you can verify that 
!    it all works as it should by checking that the first two columns
!    in the file are identical.
!
!  YOU DON'T NEED TO KNOW THIS:      
!    The resolution is defined so that the number of pixels along 
!    the side of a triangular face is 2*resolution, so there are 
!    resolution*(2*resolution+1) pixels on each of the 20 faces of 
!    the icosahedron. To avoid counting pixels on edges more than 
!    once, the edge pixels are split up half-and-half between the
!    two faces to which the edge belongs, so each face in fact
!    only contains 2*resolution*(resolution-1) pixels if you ignore the corners. 
!    The 12 corner pixels aren't lumped in with any face at all,
!    so you can see them listed separately as the last 12 pixels
!    in test.dat if you run the demo. 
!    This makes 40*resolution*(resolution-1) + 12 pixels all in all.
!    Thanks to Christopher Weth for catching typos in an earlier version
!    of this documentation!
!
!  FEEDBACK:
!    If you have any questions, comments or suggestions regarding 
!    this package, please email them to me at max@ias.edu.
!    Thanks,
!    ;-)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!     THESE SUBROUTINES ARE ALL YOU NEED TO CALL FROM OUTSIDE    !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
!!!	These subroutines convert between unit vectors and         !!!
!!!     pixel numbers, and are the only ones that the user of this !!!
!!!     package calls repeatedly:				   !!!
!!!	  subroutine vector2pixel(vector,resolution,R,v,pixel)     !!!
!!!	  subroutine pixel2vector(pixel,resolution,R,v,vector)     !!!
!!!                                                                !!!
!!!	These subroutines are called only once, in the beginning,  !!!
!!!     and compute the necessary rotation matrices and corner     !!!
!!!     vectors once and for all:                                  !!!
!!!	  subroutine compute_matrices(R)                           !!!
!!!	  subroutine compute_corners(v)                            !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	

	subroutine vector2pixel(vector,resolution,R,v,pixel)
	implicit none
	real vector(3), R(0:19,3,3), v(0:11,3)
	real A(3,3), vec(3), x, y
	integer resolution, pixel
	integer pix, face, pixperface, ifail
	if (resolut ion.lt.1) pause 'Resolution must exceed 0'
	pixperface = 2*resolution*(resolution-1)
	call find_face(vector,R,face)
	call getmatrix(face,R,A)
	call vecmatmul2(A,vector,vec)	
	x	= vec(1)/vec(3)
	y	= vec(2)/vec(3)
	call adjust(x,y)
	call tangentplanepixel(resolution,x,y,pix,ifail)
	if (ifail.gt.0) then 
	  ! Try the runner-up face:
	  call find_another_face(vector,R,face)
	  call getmatrix(face,R,A)
	  call vecmatmul2(A,vector,vec)	
	  x	= vec(1)/vec(3)
	  y	= vec(2)/vec(3)
	  call adjust(x,y)
	  call tangentplanepixel(resolution,x,y,pix,ifail)
	end if	
	pixel = face*pixperface + pix
	if (ifail.gt.0) then
	  ! The pixel wasn't in any of those two faces,
	  ! so it must be a corner pixel.
	  call find_corner(vector,v,pix)
	  pixel = 20*pixperface + pix
	end if
	return
	end

	subroutine pixel2vector(pixel,resolution,R,v,vector)
	! Returns a unit vector pointing towards pixel.
	! Resolution must be an even, positive integer.
	implicit none
	real vector(3), R(0:19,3,3), v(0:11,3)
	real A(3,3), x, y, norm
	integer resolution, pixel
	integer pix, face, pixperface
	if (resolution.lt.1) pause 'Resolution must exceed 0'
	pixperface = 2*resolution*(resolution-1)
	if (pixel.lt.0) pause 'Error: negative pixel number'	
	if (pixel.ge.20*pixperface+12) 
     &	  pause 'Error: pixel number too large'
	if (pixperface.gt.0) then 
	  face = pixel/pixperface
	  if (face.gt.20) face = 20
	else ! There are no pixels at all on the faces - it's all just corners.
	  face = 20
	end if	
	pix = pixel - face*pixperface
	if (face.lt.20) then
	  ! The pixel is on one of the 20 faces:
	  call tangentplanevector(pix,resolution,x,y)
	  call unadjust(x,y)
	  norm 	= sqrt(x*x+y*y+1)
	  vector(1)	= x/norm
	  vector(2)	= y/norm
	  vector(3)	= 1./norm
	  call getmatrix(face,R,A)
	  call vecmatmul1(A,vector,vector)
	else
	  ! This is a corner pixel:
	  if (pix.gt.11) pause 'Error: pixel number too big'
	  vector(1) = v(pix,1)
	  vector(2) = v(pix,2)
	  vector(3) = v(pix,3)
	end if
	return
	end

	subroutine compute_matrices(R)
	! On exit, R will contain the 20 rotation matrices
	! that rotate the 20 icosahedron faces 
	! into the tangent plane
	! (the horizontal plane with z=1). 
	! Only called once, so speed is irrelevant.
	implicit none
	real R(0:19,3,3)
	real A(3,3), B(3,3), C(3,3), D(3,3), E(3,3)
	real pi, sn, cs, ct, x
	integer i,j,n
	do i=1,3
	  do j=1,3
	    A(i,j) = 0.
	    B(i,j) = 0.
	    C(i,j) = 0.
	    D(i,j) = 0.
	  end do
	end do
	pi 	= 4.*atan(1.)
	x 	= 2.*pi/5.
	cs 	= cos(x)
	sn	= sin(x)
	A(1,1)	= cs
	A(1,2)	= -sn
	A(2,1)	= sn
	A(2,2)	= cs
	A(3,3) 	= 1.
	! A rotates by 72 degrees around the z-axis.
	x 		= pi/5.
	ct 		= cos(x)/sin(x)
	cs		= ct/sqrt(3.)
	sn		= sqrt(1-ct*ct/3.)
	C(1,1)	= 1
	C(2,2)	= cs
	C(2,3)	= -sn
	C(3,2)	= sn
	C(3,3) 	= cs
	! C rotates around the x-axis so that the north pole 	
	! ends up at the center of face 1.
	cs		= -0.5
	sn		= sqrt(3.)/2
	D(1,1)	= cs
	D(1,2)	= -sn
	D(2,1)	= sn
	D(2,2)	= cs
	D(3,3)	= 1.
	! D rotates by 120 degrees around z-axis.
	call matmul1(C,D,E)
	call matmul2(E,C,B)	! B = CDC^t
	! B rotates face 1 by 120 degrees. 
	do i=1,3
	  do j=1,3
	    E(i,j) = 0.
	  end do
	  E(i,i) = 1.
	end do	! Now E is the identity matrix.
	call putmatrix(0,R,E)	
	call matmul1(B,A,E)
	call matmul1(B,E,E)
	call putmatrix(5,R,E)	
	call matmul1(E,A,E)
	call putmatrix(10,R,E)
	call matmul1(E,B,E)
	call matmul1(E,B,E)
	call matmul1(E,A,E)
	call putmatrix(15,R,E)
	do n=0,15,5	
	  call getmatrix(n,R,E)
	  do i=1,4
	    call matmul1(A,E,E)
	    call putmatrix(n+i,R,E)
	  end do
	end do
	! Now the nth matrix in R will rotate 
	! face 1 into face n. 
	! Multiply by C so that they will rotate
	! the tangent plane into face n instead:
	do n=0,19
	  call getmatrix(n,R,E)
	  call matmul1(E,C,E)
	  call putmatrix(n,R,E)
	end do
	return
	end

	subroutine compute_corners(v)
	! On exit, v will contain unit vectors pointing toward 
	! the 12 icoshedron corners. 
	implicit none
	real v(0:11,3)
	real pi, z, rho, dphi
	integer i
	pi = 4.*atan(1.)
	dphi = 2.*pi/5.
	! First corner is at north pole:
	v(0,1) = 0.
	v(0,2) = 0.
	v(0,3) = 1.
	! The next five lie on a circle, with one on the y-axis:
	z = 0.447213595	! This is 1/(2 sin^2(pi/5)) - 1
	rho = sqrt(1.-z*z)
	do i=0,4
	  v(1+i,1) = -rho*sin(i*dphi)
	  v(1+i,2) =  rho*cos(i*dphi)
	  v(1+i,3) = z 
	end do
	! The 2nd half are simply opposite the first half:
	do i=0,5
	  v(6+i,1) = -v(i,1)
	  v(6+i,2) = -v(i,2)
	  v(6+i,3) = -v(i,3)
	end do
	return
	end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!     THE SUBROUTINES BELOW ARE SUBORDINATE TO THOSE ABOVE, AND  !!!
!!!     CAN BE SAFELY IGNORED BY THE GENERAL USER OF THE PACKAGE.  !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
!!!	These subroutines perform some standard linear algebra:    !!!
!!!	  subroutine matmul1(A,B,C)                                !!!
!!!	  subroutine matmul2(A,B,C)                                !!!	
!!!	  subroutine matmul3(A,B,C)                                !!!	
!!!	  subroutine vecmatmul1(A,b,c)	                           !!!	
!!!	  subroutine vecmatmul2(A,b,c)	                           !!!	
!!!                                                                !!!
!!!	These subroutines copy matrices in and out of storage:     !!!
!!!	  subroutine getmatrix(n,R,A)                              !!!
!!!	  subroutine putmatrix(n,R,A)                              !!!
!!!                                                                !!!
!!!     These subroutines help vector2pixel reduce the 3D sphere   !!!
!!!     problem to a problem on an equilateral triangle in the     !!!
!!!     z=1 tangent plane (an icosahedron face):                   !!!		
!!!	  subroutine find_face(vector,R,face)                      !!!
!!!	  subroutine find_another_face(vector,R,face)              !!!
!!!	  subroutine find_corner(vector,v,corner)                  !!!
!!!                                                                !!!
!!!     These subroutines pixelize this triangle with a regular    !!!
!!!     triangular grid:                                           !!!
!!!	  subroutine find_mn(pixel,resolution,m,n)                 !!!
!!!	  subroutine tangentplanepixel(resolution,x,y,pix,ifail)   !!!
!!!	  subroutine tangentplanevector(pix,resolution,x,y)        !!!
!!!                                                                !!!
!!!     These subroutines reduce the area equalization problem to  !!!
!!!     one on the right triangle in the lower right corner:       !!!
!!!	  subroutine find_sixth(x,y,rot,flip)                      !!!
!!!	  subroutine rotate_and_flip(rot,flip,x,y)	           !!!
!!!	  subroutine adjust(x,y)                                   !!!
!!!	  subroutine unadjust(x,y)                                 !!!
!!!                                                                !!!
!!!     These subroutines perform the area equalization mappings   !!!
!!!     on this right triangle:                                    !!!
!!!	  subroutine adjust_sixth(x,y)                             !!!
!!!	  subroutine unadjust_sixth(x,y)                           !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	

	subroutine matmul1(A,B,C)	
	! Matrix multiplication C = AB.
	! A, B and C are allowed to be physiclly the same.
	implicit none
	real A(3,3), B(3,3), C(3,3), D(3,3), sum
	integer i,j,k
	sum = 0.
	do i=1,3
	  do j=1,3
	    sum = 0.
	    do k=1,3
	      sum = sum + A(i,k)*B(k,j)
	    end do
	    D(i,j) = sum
	  end do
	end do
	call copymatrix(D,C)
	return
	end

	subroutine matmul2(A,B,C)	
	! Matrix multiplication C = AB^t
	! A, B and C are allowed to be physically the same.
	implicit none
	real A(3,3), B(3,3), C(3,3), D(3,3), sum
	integer i,j,k
	sum = 0.
	do i=1,3
	  do j=1,3
	    sum = 0.
	    do k=1,3
	      sum = sum + A(i,k)*B(j,k)
	    end do
	    D(i,j) = sum
	  end do
	end do
	call copymatrix(D,C)
	return
	end

	subroutine matmul3(A,B,C)	
	! Matrix multiplication C = A^t B
	! A, B and C are allowed to be physically the same.
	implicit none
	real A(3,3), B(3,3), C(3,3), D(3,3), sum
	integer i,j,k
	sum = 0.
	do i=1,3
	  do j=1,3
	    sum = 0.
	    do k=1,3
	      sum = sum + A(k,i)*B(k,j)
	    end do
	    D(i,j) = sum
	  end do
	end do
	call copymatrix(D,C)
	return
	end

	subroutine vecmatmul1(A,b,c)	
	! Matrix multiplication c = Ab
	! b and c are allowed to be physically the same.
	implicit none
	real A(3,3), b(3), c(3), d(3), sum
	integer i,j
	sum = 0.
	do i=1,3
	  sum = 0.
	  do j=1,3
	    sum = sum + A(i,j)*b(j)
	  end do
	  d(i) = sum
	end do
	call copyvector(d,c)
	return
	end

	subroutine vecmatmul2(A,b,c)	
	! Matrix multiplication c = A^tb
	! b and c are allowed to be physiclly the same.
	implicit none
	real A(3,3), b(3), c(3), d(3), sum
	integer i,j
	sum = 0.
	do i=1,3
	  sum = 0.
	  do j=1,3
	    sum = sum + A(j,i)*b(j)
	  end do
	  d(i) = sum
	end do
	call copyvector(d,c)
	return
	end

	subroutine copymatrix(A,B)	
	! B = A
	implicit none
	real A(3,3), B(3,3)
	integer i,j
	do i=1,3
	  do j=1,3
	    B(i,j) = A(i,j)
	  end do
	end do
	return
	end

	subroutine copyvector(a,b)	
	! b = a
	implicit none
	real a(3), b(3)
	integer i
	do i=1,3
	  b(i) = a(i)
	end do
	return
	end

	subroutine getmatrix(n,R,A)	
	! A = the nth matrix in R
	implicit none
	real R(0:19,3,3), A(3,3)
	integer i,j,n
	do i=1,3
	  do j=1,3
	    A(i,j) = R(n,i,j)
	  end do
	end do
	return
	end

	subroutine putmatrix(n,R,A)	
	! the nth matrix in R = A
	implicit none
	real R(0:19,3,3), A(3,3)
	integer i,j,n
	do i=1,3
	  do j=1,3
	    R(n,i,j) = A(i,j)
	  end do
	end do
	return
	end

	subroutine find_face(vector,R,face)
	! Locates the face to which vector points.
	! Computes the dot product with the vectors
	! pointing to the center of each face and picks the
	! largest one. 
	! This simple routine can be substantially accelerated
	! by adding a bunch of if-statements, to avoid looping
	! over more than a few faces.
	implicit none
	real vector(3), R(0:19,3,3), dot, max
	integer n,face,i
	max = -17.
	do n=0,19
	  dot = 0.
	  do i=1,3
	    dot = dot + R(n,i,3)*vector(i)
	  end do
	  if (dot.gt.max) then
	    face = n
	    max = dot
	  end if
	end do
	return
	end

	subroutine find_another_face(vector,R,face)
	! Computes the dot product with the vectors
	! pointing to the center of each face and picks the
	! largest one other than face.
	! This simple routine can be substantially accelerated
	! by adding a bunch of if-statements, to avoid looping
	! over more than a few faces.
	implicit none
	real vector(3), R(0:19,3,3), dot, max
	integer n,face,facetoavoid,i
	facetoavoid = face
	max = -17.
	do n=0,19
	  if (n.ne.facetoavoid) then 
	    dot = 0.
	    do i=1,3
	      dot = dot + R(n,i,3)*vector(i)
	    end do
	    if (dot.gt.max) then
	      face = n
	      max = dot
	    end if
	  end if
	end do
	return
	end

	subroutine find_corner(vector,v,corner)
	! Locates the corner to which vector points.
	! Computes the dot product with the vectors
	! pointing to each corner and picks the
	! largest one. 
	! This simple routine can be substantially accelerated
	! by adding a bunch of if-statements, but that's pretty
	! pointless since it gets called so rarely.
	implicit none
	real vector(3), v(0:11,3), dot, max
	integer corner,n,i
	max = -17.
	do n=0,11
	  dot = 0.
	  do i=1,3
	    dot = dot + v(n,i)*vector(i)
	  end do
	  if (dot.gt.max) then
	    corner = n
	    max = dot
	  end if
	end do
	return
	end

	subroutine find_mn(pixel,resolution,m,n)
	! Computes the integer coordinates (m,n) of the pixel 
	! numbered pix on the basic triangle.
	implicit none
	integer pixel, pix, resolution, m, n
	integer interiorpix , pixperedge
	pix 	    = pixel
	interiorpix = (2*resolution-3)*(resolution-1)
	pixperedge  = (resolution)-1
	if (pix.lt.interiorpix) then 
	  ! The pixel lies in the interior of the triangle.
	  m = (sqrt(1.+8.*pix)-1.)/2. + 0.5/resolution
	  ! 0.5/resolution was added to avoid problems with
	  ! rounding errors for the case when n=0.
	  ! As long as you don't add more than 2/m, you're OK.
	  n = pix - m*(m+1)/2
	  m = m + 2
	  n = n + 1
	  goto 555
	end if
	pix = pix - interiorpix 
	if (pix.lt.pixperedge) then
	  ! The pixel lies on the bottom edge.
	  m = 2*resolution-1
	  n = pix+1
	  goto 555
	end if
	pix = pix - pixperedge
	if (pix.lt.pixperedge) then
	  ! The pixel lies on the right edge.
	  m = 2*resolution-(pix+2)
	  n = m
	  goto 555
	end if
	pix = pix - pixperedge
	! The pixel lies on the left edge.
	m = pix+1
	n = 0
555	return
	end

	subroutine tangentplanepixel(resolution,x,y,pix,ifail)
	! Finds the hexagon in which the point (x,y) lies
	! and computes the corresponding pixel number pix.
	! Returns ifail=0 if (x,y) lies on the face, 
	! otherwise returns ifail=1.
	implicit none
	real x, y, a, b, c, d, edgelength
	parameter (c=0.866025404)	! sqrt(3)/2
	parameter(edgelength=1.3231690765)
	! The edge length of the icosahedron is 
	! sqrt(9 tan^2(pi/5) - 3) when scaled so that
	! it circumscribes the unit sphere. 
	integer resolution, pix, ifail, i, j, k, m, n, r2
	r2	= 2*resolution
	a 	= 0.5*x
	b 	= c*y
	!d 	= 0.5*edgelength/r2
	d 	= 0.5*edgelength/(r2-1)	! Fixed by Robert Berry 5/2-13
	i 	= x/d 	  + r2
	j 	= (a+b)/d + r2
	k 	= (a-b)/d + r2
	m 	= (r2+r2-j+k-1)/3
	n 	= (i+k+1-r2)/3
	pix 	= (m-2)*(m-1)/2 + (n-1)
	ifail = 0
	if (m.eq.r2-1) then		! On bottom row
	  if ((n.le.0).or.(n.ge.resolution)) then
	    ifail=1
	  end if
	  goto 666				! Pix already correct
	end if
	if (n.eq.m) then			! On right edge
	  k = (r2-1) - m
	  if ((k.le.0).or.(k.ge.resolution)) then
	    ifail = 1
	  else
	    pix = (r2-2)*(resolution-1) + k - 1
 	  end if
	  goto 666
	end if
	if (n.eq.0) then			! On left edge
	  if ((m.le.0).or.(m.ge.resolution)) then
	    ifail = 1
	  else
	    pix = (r2-1)*(resolution-1) + m - 1
	  end if
	end if
666	return
	end

	subroutine tangentplanevector(pix,resolution,x,y)
	! Computes the coordinates (x,y) of the pixel 
	! numbered pix on the basic triangle.
	implicit none
	real x, y, c1, c2, edgelength
	parameter(c1=0.577350269)	! 1/sqrt(3)
	parameter(c2=0.866025404)	! sqrt(3)/2
	parameter(edgelength=1.3231690765)
	! The edge length of the icosahedron is 
	! sqrt(9 tan^2(pi/5) - 3) when scaled so that
	! it circumscribes the unit sphere. 
	integer pix, resolution, m, n
	call find_mn(pix,resolution,m,n)
	x	= edgelength*(n-0.5*m)/(2*resolution-1)
	y 	= edgelength*(c1-(c2/(2*resolution-1))*m)
	return
	end

	subroutine find_sixth(x,y,rot,flip)
	! Find out in which sixth of the basic triangle
	! the point (x,y) lies, identified by the
	! two integers rot (=0, 1 or 2) and flip = 0 or 1).
	! rot and flip are defined such that the sixth is 
	! mapped onto the one at the bottom right by
	! these two steps:
	! 1. Rotate by 120 degrees anti-clockwise, rot times.
	! 2. Flip the sign of x if flip = 1, not if flip=0.
	! The if-statements below go through the six cases 
	! anti-clockwise, starting at the bottom right.
	implicit none
	real x, y, c, d
	parameter(c=1.73205081)		! sqrt(3)
	integer rot, flip
	d = c*y
	if (x.ge.0) then
	  if (x.le.-d) then
	    rot  = 0
	    flip = 0
	  else
	    if (x.ge.d) then
	      rot  = 2
	      flip = 1
	    else 
	      rot  = 2
	      flip = 0
	    end if
	  end if
	else
	  if (x.ge.-d) then
	    rot  = 1
	    flip = 1
	  else
	    if (x.le.d) then
	      rot  = 1
	      flip = 0
	    else 
	      rot  = 0
	      flip = 1
	    end if
	  end if
	end if
	return
	end

	subroutine rotate_and_flip(rot,flip,x,y)
	implicit none
	real x, y, x1, cs, sn, c
	integer rot, flip
	parameter(cs=-0.5)
	parameter(c=0.866025404)	! sqrt(3)/2
	if (rot.gt.0) then
	  if (rot.eq.1) then
	    sn = c	! Rotate 120 degrees anti-clockwise 
	  else 
	    sn = -c	! Rotate 120 degrees anti-clockwise 
	  end if
	  x1 = x
	  x  = cs*x1 - sn*y
	  y  = sn*x1 + cs*y
	end if
	if (flip.gt.0) x = -x
	return
	end

	subroutine adjust(x,y)
	! Maps the basic triangle onto itself in such a way
	! that pixels will have equal area when mapped onto
	! the sphere. 
	implicit none
	real x, y
	integer rot, flip
	call find_sixth(x,y,rot,flip)
	call rotate_and_flip(rot,flip,x,y)
	call adjust_sixth(x,y)
	! Now rotate & flip the sixth back into its
	! original position:
	if ((flip.eq.0).and.(rot.gt.0)) then  
	  call rotate_and_flip(3-rot,flip,x,y)
	else
	  call rotate_and_flip(rot,flip,x,y)
	end if
	return
	end
	
	subroutine unadjust(x,y)
	! Performs the inverse of what adjust does. 
	implicit none
	real x, y
	integer rot, flip
	call find_sixth(x,y,rot,flip)
	call rotate_and_flip(rot,flip,x,y)
	call unadjust_sixth(x,y)
	! Now rotate & flip the sixth back into its
	! original position:
	if ((flip.eq.0).and.(rot.gt.0)) then  
	  call rotate_and_flip(3-rot,flip,x,y)
	else
	  call rotate_and_flip(rot,flip,x,y)
	end if
	return
	end
	
	subroutine adjust_sixth(x,y)
	! Maps the basic right triangle (the sixth of the face that
	! is in the lower right corner) onto itself in such a way
	! that pixels will have equal area when mapped onto the sphere. 
	implicit none
	real x, y, u, v, g, v2, root, trig, eps, scale
	parameter(eps=1.e-14, scale=1.09844)
	parameter(g=1.7320508075689)		! sqrt(3)
	u 		= x  + eps
	v		= -y + eps
	v2		= v*v
	root		= sqrt(1.+4.*v2)
	trig		= atan((g*root-g)/(root+3.))
	y		= sqrt(trig*2./g)
	x		= sqrt((1.+4.*v2)/(1.+u*u+v2))*u*y/v
	x		= scale*x
	y 		= -scale*y	
	return
	end
	
	subroutine unadjust_sixth(x,y)
	! Performs the inverse of what adjust_sixth does.
	implicit none
	real x, y, u, v, g, v2, y2, tmp, trig, eps, scale
	parameter(eps=1.e-14, scale=1.09844)
	parameter(g=1.7320508075689)		! sqrt(3)
	u 		=  x/scale + eps
	v		= -y/scale + eps
	v2		= v*v
	trig		= tan(g*v2/2.)
	tmp		= (g+3.*trig)/(g-trig)
	y2		= (tmp*tmp-1.)/4.
	y		= sqrt(y2)
	tmp		= v2*(1.+4.*y2) - u*u*y2
	x	 	= u*y*sqrt((1.+y2)/tmp) 
	y 		= -y	
	return
	end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!     END OF THE ICOSAHEDRON PACKAGE FOR PIXELIZING THE SPHERE   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	

