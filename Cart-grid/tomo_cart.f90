!---------------------------------------------------------------------------------------------
! <Ne3dTomo> provides the tomographic reconstruction of 3D electron density of the corona
! with an option by applying radial weighting factor to the regularization term to provide 
! balanced smoothing at different heights and help stablizing the solution. The detailed 
! method and applications are referred to the publication "Wang, T., Arge, C.N, & Jones, I.J.
! (2025), Sol.Phys., in Press. Copyright (C) 2025 by Tongjiang Wang (wangt@cua.edu)
!
! First version released on March 28, 2025
! (download available from https://stereo-ssc.nascom.nasa.gov/data/ins_data/secchi/N3D_COR1B)
!
! This program is free software: you can redistribute it and/or modify it under the terms of the 
!    GNU General Public License as published by the Free Software Foundation, either version 3 of the 
!    License, or (at your option) any later version.
! This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
!    without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
!    See the GNU General Public License for more details.
!-------------------------------------------------------------------------------------------
MODULE InternPro 
  IMPLICIT NONE
 
   INTEGER, PARAMETER :: double=SELECTED_REAL_KIND(15)   
   REAL (kind=double), PARAMETER:: PI=3.14159265358979D0
   REAL, PARAMETER:: radeg=180.D0/PI

   INTERFACE img2car_coord
      MODULE PROCEDURE img2car_coord_1d,img2car_coord_scl
   END INTERFACE

   INTERFACE looktable
      MODULE PROCEDURE looktable_R, looktable_I
   END INTERFACE

   interface findmax
       module procedure  findmax_I, findmax_R, findmax_D
   end interface

   interface findmin
       module procedure  findmin_I, findmin_R, findmin_D
   end interface


   CONTAINS


! DTIME in standard Fortran.
Function dtime(time)
  Real, intent(out) ::  time(2)
  Real :: dtime
  Double Precision,Save :: last_time = 0
  Double Precision :: this_time
  Intrinsic Cpu_Time

  Call Cpu_Time(this_time)
  time(1) = this_time - last_time
  time(2) = 0
  dtime = time(1)
  last_time = this_time
End Function dtime


SUBROUTINE locate(xx,x,j)
   IMPLICIT NONE
   REAL (kind=double), INTENT(IN) :: xx(:),x
   INTEGER, INTENT(OUT) :: j
   INTEGER :: n, jl, ju, jm
! Given an array xx(1:n), and a given value x, return a value j such
!that x is between xx(j) and xx(j+1). xx(1:n) must be monotonic, either
!increasing or decreasing. j=0 or j=n is return to indicate that x is out of
!range
n=size(xx)
jl=0
ju=n+1
DO WHILE (ju-jl >1)
   jm=(ju+jl)/2
   IF (xx(n) >= xx(1) .eqv. x >= xx(jm)) THEN
      jl=jm
   ELSE
      ju=jm
   ENDIF
ENDDO
IF (x == xx(1)) THEN
  j=1
ELSE IF (x == xx(n)) THEN
  j=n-1
ELSE
  j=jl
ENDIF
END SUBROUTINE locate

SUBROUTINE interp(xData, yData, xVal, yVal )
! Inputs: xData = a vector of the x-values of the data to be interpolated
! (monotonically increasing)
!         yData = a vector of the y-values of the data to be interpolated
!         xVal  = a vector of the x-values where interpolation should be
!         performed (monotonic is not necessary)
! Return yVal  = a vector of the resulting interpolated values

  implicit none

  real (kind=double), intent(in), dimension(:) :: xData, yData
  real (kind=double), intent(in), dimension(:) :: xVal
  real (kind=double), intent(out), dimension(size(xVal)) :: yVal

  integer :: inputIndex, dataIndex
  real (kind=double) ::  weight

  ! Possible checks on inputs could go here
  ! Things you may want to check:  monotonically increasing xData
  DO dataIndex=1, size(xData)-1
     IF (xData(dataIndex+1) < xData(dataIndex)) THEN
        print *,"xData must increase monotonically for INTERP()!"
        STOP
     ENDIF
  ENDDO

  IF  (size(xData) /= size(yData)) THEN
     print *, "xData and yData must have the same size! "
     STOP
  ENDIF
  !   size(xVal) == size(yVal)

  do inputIndex = 1, size(xVal)
       call locate(xData, xVal(inputIndex), dataIndex)
       weight = (xVal(inputIndex) - xData(dataIndex))/(xData(dataIndex+1)-xData(dataIndex))
       yVal(inputIndex) = (1.0-weight)*yData(dataIndex) + weight*yData(dataIndex+1)
  end do
END SUBROUTINE interp

FUNCTION  looktable_R (xData, yData, xVal, state) RESULT(yVal)
! purpose: search index of given x-values in xData, and return the
! values in yData for the corresponding index. return yVal(index)=0 and set
! state(0)=number of elements not found; state(index)=index of x-Val not found
! Inputs: xData = a vector of the x-values of the data to be searched
! (monotonically increasing)
!         yData = a vector of the y-values of the data to be interpolated
!         xVal  = a vector of the x-values where interpolation should be
!         performed (monotonic is not necessary)
! Return yVal  = a vector of the resulting searched values
  implicit none

  real, parameter :: eps=1.d-5
  real (kind=double), intent(in), dimension(:) :: xData, yData, xVal
  real (kind=double), dimension(size(xVal)) :: yVal
  integer, intent(out), allocatable :: state(:)
  integer, dimension(0:size(xVal)) :: stat0
  integer :: inputIndex, DataIndex, id, nnf, i,k

  ! Possible checks on inputs could go here
  ! Things you may want to check:  monotonically increasing xData
  DO dataIndex=1, size(xData)-1
     IF (xData(dataIndex+1) < xData(dataIndex)) THEN
        print *,"xData must increase monotonically for LOOKTABLE_R()!"
        STOP
     ENDIF
  ENDDO

  IF  (size(xData) /= size(yData)) THEN
     print *, "xData and yData must have the same size! "
     STOP
  ENDIF
  nnf=0
  do inputIndex = 1, size(xVal)
       call locate(xData, xVal(inputIndex), id)
       if (abs(xData(id)-xVal(inputIndex)) < eps) then
             yVal(inputIndex)=yData(id)
             stat0(inputIndex)=0
       else if (abs(xData(id+1)-xVal(inputIndex)) < eps) then
             yVal(inputIndex)=yData(id+1)
             stat0(inputIndex)=0
       else
             yVal(inputIndex)=0
             stat0(inputIndex)=inputIndex
             nnf=nnf+1
!        print*,'--- x-Value=',xVal(inputIndex),' is not found!'
       endif
  end do
  stat0(0)=nnf

if(nnf >0) then
    allocate(state(0:nnf))
    k=1
    do i=1,size(xVal)
       if (stat0(i)/=0) then
          stat0(k)=stat0(i)
          k=k+1
       endif
    enddo
    if (nnf /= k-1) stop "There is an error in looktable"
    state=stat0(0:nnf)
    print*,'Number of non-found elements with looktable_R=',state(0)
else
   allocate(state(0:0))
   state(0)=0
endif
return
END FUNCTION looktable_R

FUNCTION looktable_I (xData, yData, xVal, state) RESULT(yVal)
! input and output variables are integer types
  implicit none

  integer, intent(in), dimension(:) :: xData, yData, xVal
  integer, dimension(size(xVal)) :: yVal
  integer, intent(out), allocatable :: state(:)
  integer, dimension(0:size(xVal)) :: stat0
  integer :: inputIndex, DataIndex, id, nnf, i,k

  ! Possible checks on inputs could go here
  ! Things you may want to check:  monotonically increasing xData
  DO dataIndex=1, size(xData)-1
     IF (xData(dataIndex+1) < xData(dataIndex)) THEN
        print *,"xData must increase monotonically for LOOKTABLE_I!"
        print*,'xData=',xData
        STOP
     ENDIF
  ENDDO

  IF  (size(xData) /= size(yData)) THEN
     print *, "xData and yData must have the same size! "
     STOP
  ENDIF
  nnf=0
  do inputIndex = 1, size(xVal)
       call locate(dble(xData), dble(xVal(inputIndex)), id)
       if (xData(id) == xVal(inputIndex)) then
             yVal(inputIndex)=yData(id)
             stat0(inputIndex)=0
       else if (xData(id+1) == xVal(inputIndex)) then
             yVal(inputIndex)=yData(id+1)
             stat0(inputIndex)=0
       else
             yVal(inputIndex)=0
             stat0(inputIndex)=inputIndex
             nnf=nnf+1
       endif
  end do
  stat0(0)=nnf

if(nnf >0) then
    allocate(state(0:nnf))
    k=1
    do i=1,size(xVal)
       if (stat0(i)/=0) then
          stat0(k)=stat0(i)
          k=k+1
       endif
    enddo
    if (nnf /= k-1) stop "There is an error in looktable"
    state=stat0(0:nnf)
    print*,'Number of un-found elements with looktable_I=',state(0)
else
   allocate(state(0:0))
   state(0)=0
endif
return
END FUNCTION looktable_I


SUBROUTINE xy2sp(x,y,z,Bd0,Pd,Ld,Bd,R)
!+ Purpose
! transform from helicentric-cartesian coord. 
!  [(x,y) is image plane, z-toward to observer
! into Stonyhurst heliographic coord
! (same as Carr.coord. but with Lon. at meridian=0, L in [-180, 180]
!Input:  x,y,z - cartesian coord
!Output: Ld,Bd,R - Lon, Lat (deg), and Rsun (arcsec)
! Keyword
! Bd0 - Carrington Latitude of solar disk center
! Pd  - position angle of solar north, positive when 

   IMPLICIT NONE 
   REAL, INTENT(IN), DIMENSION(:) :: x,y,z
   REAL, INTENT(IN) :: Bd0,Pd
   REAL, INTENT(OUT),DIMENSION(SIZE(x)) :: Ld,Bd,R

   REAL :: B0,P
   REAL, DIMENSION(SIZE(x)) :: x1,y1,z1,x2,y2,z2

 B0=Bd0/radeg
 P=Pd/radeg
 R=sqrt(x*x+y*y+z*z)
 x1=x*cos(P)+y*sin(P)
 y1=y*cos(P)-x*sin(P)
 z1=z
 x2=x1
 y2=y1*cos(B0)+z1*sin(B0)
 z2=z1*cos(B0)-y1*sin(B0)
 Bd=radeg*asin(y2/R)
 Ld=radeg*atan2(x2,z2)

END SUBROUTINE xy2sp

SUBROUTINE img2car_coord_1d (cobs, xp, yp, zp, x, y, z)
!+purpose: Conversion of image coordinates (xp,yp,zp)
! into the Carrington coordinates (x,y,z). The image coordinate
! system has y-axis upward to the solar North and z-axis
! toward the observer (opposite to LOS). The Carrington coordinate
! system has z-axis along solar axis toward North, and xy-plane in
! solar equatoral plane. Two systems have the same origin.
!  Carrington longitude (lon) starts from
! x-axis, positive for counter-clockwise. Carrington latitude
! is relative to xy-plane, positive toward z-axis.
! Assuming the Carrington coordinate of LOS is (phi, 90-theta),
! where theta is the angle of LOS to z-axis, then
! the Carrington coordinate of observer (or satellite) (lon, lat)
! is: lon=phi-180, lat=theta-90
!Input:
! cobs[2] - including Carrington coordinate of observer, [lon, lat]
!           in degree
! xp,yp,zp -  array for image coord.
!Output:
! x,y,z -  array for Carrington coord.
IMPLICIT NONE
REAL (kind=double), INTENT(IN), DIMENSION(:) :: xp, yp, zp
REAL (kind=double), INTENT(IN), DIMENSION(2) :: cobs
REAL (kind=double), INTENT(OUT), DIMENSION(SIZE(xp)) :: x, y, z
REAL (kind=double) :: t11,t12,t13,t21,t22,t23,t31,t32,t33
REAL (kind=double) :: phi, theta

!position of LOS
  phi  =(cobs(1)+180.d0)/radeg
  theta=(cobs(2)+90.d0)/radeg
! transform maxtrix
 t11=  sin(phi)
 t12= -cos(theta)*cos(phi)
 t13= -sin(theta)*cos(phi)
 t21= -cos(phi)
 t22= -cos(theta)*sin(phi)
 t23= -sin(theta)*sin(phi)
 t31= 0.d0
 t32=  sin(theta)
 t33= -cos(theta)
 !from [xp,yp,zp] to [x,y,z]
 x=t11*xp + t12*yp +t13*zp
 y=t21*xp + t22*yp +t23*zp
 z=t31*xp + t32*yp +t33*zp
 END SUBROUTINE img2car_coord_1d


SUBROUTINE img2car_coord_scl (cobs, xp, yp, zp, x, y, z)
! xp,yp,zp - scale  for image coord.
!Output:
! x,y,z - scale for Carrington coord.
IMPLICIT NONE
REAL (kind=double), INTENT(IN) :: xp, yp, zp
REAL (kind=double), INTENT(IN), DIMENSION(2) :: cobs
REAL (kind=double), INTENT(OUT) :: x, y, z
REAL (kind=double) :: t11,t12,t13,t21,t22,t23,t31,t32,t33
REAL (kind=double) :: phi, theta

!position of LOS
  phi  =(cobs(1)+180.d0)/radeg
  theta=(cobs(2)+90.d0)/radeg
! transform maxtrix
 t11=  sin(phi)
 t12= -cos(theta)*cos(phi)
 t13= -sin(theta)*cos(phi)
 t21= -cos(phi)
 t22= -cos(theta)*sin(phi)
 t23= -sin(theta)*sin(phi)
 t31= 0.d0
 t32=  sin(theta)
 t33= -cos(theta)
 !from [xp,yp,zp] to [x,y,z]
 x=t11*xp + t12*yp +t13*zp
 y=t21*xp + t22*yp +t23*zp
 z=t31*xp + t32*yp +t33*zp
 END SUBROUTINE img2car_coord_scl


SUBROUTINE car2sp (x,y,z,lon,lat,r)
!+ Purpose
! transform Carrington coord. from Cartesian sys. to Spheric sys.
!Input:  x,y,z - cartesian coord. z defined as solar axis toward N-pole
!      [x,y] in equatial plane
!Output:
! lon - longitude in rad, range in [-pi, pi], measured from x-axis;
!       positive in anti-clockwise
! lat - latitude in rad, range in [-pi/2, pi/2]
! r   - heliocentric distance in unit same as x,y,or z
!-
IMPLICIT NONE
REAL (kind=double), INTENT(IN), DIMENSION(:) :: x,y,z
REAL (kind=double), INTENT(OUT), DIMENSION(size(x)) :: lon,lat,r

 r=sqrt(x*x+y*y+z*z)
 lat=asin(z/r)
 lon=atan2(y,x)
END SUBROUTINE car2sp

  function findmax_D (a,w,n) result(amax)
 !   a(:) is 1D vector
 !   w(:) 1D vector including indice of maximum elements
 !   n is the number of maximum elements in a(:)
      implicit none
      real (kind=8), intent(in), dimension (:) :: a
      integer, intent(out), allocatable :: w(:)
      integer :: n, mask(size(a)), i,k
      real (kind=8) :: amax

      amax=maxval(a)
      where (a == amax)
         mask=1
      elsewhere
         mask=0
      endwhere
      n=sum(mask)
      allocate(w(n))
      k=1
      do  i=1,size(a)
        if (mask(i)==1) then
          w(k)=i
          k=k+1
        endif
      enddo
!      print*,'max value=',amax
!      print*,'number of max elements=',n
!      print*,'max location=',w
!     if (n/=(k-1)) STOP "--error!"
  return
  end function findmax_D


  function findmax_R (a,w,n) result(amax)
 !   a(:) is 1D vector
 !   w(:) 1D vector including indice of maximum elements
 !   n is the number of maximum elements in a(:)
      implicit none
      real (kind=4), intent(in), dimension (:) :: a
      integer, intent(out), allocatable :: w(:)
      integer :: n, mask(size(a)), i,k
      real (kind=4) :: amax

      amax=maxval(a)
      where (a == amax)
         mask=1
      elsewhere
         mask=0
      endwhere
      n=sum(mask)
      allocate(w(n))
      k=1
      do  i=1,size(a)
        if (mask(i)==1) then
          w(k)=i
          k=k+1
        endif
      enddo
  return
  end function findmax_R

  function findmax_I (a,w,n) result(amax)
 !   a(:) is 1D vector
 !   w(:) 1D vector including indice of maximum elements
 !   n is the number of maximum elements in a(:)
      implicit none
      integer, intent(in), dimension (:) :: a
      integer, intent(out), allocatable :: w(:)
      integer :: n, mask(size(a)), i,k
      integer :: amax

      amax=maxval(a)
      where (a == amax)
         mask=1
      elsewhere
         mask=0
      endwhere
      n=sum(mask)
      allocate(w(n))
      k=1
      do  i=1,size(a)
        if (mask(i)==1) then
          w(k)=i
          k=k+1
        endif
      enddo
  return
  end function findmax_I

  function findmin_D (a,w,n) result(amin)
 !   a(:) is 1D vector
 !   w(:) 1D vector including indice of minimum elements
 !   n is the number of minimum elements in a(:)
      implicit none
      real (kind=8), intent(in), dimension (:) :: a
      integer, intent(out), allocatable :: w(:)
      integer :: n, mask(size(a)), i,k
      real (kind=8) :: amin

      amin=minval(a)
      where (a == amin)
         mask=1
      elsewhere
         mask=0
      endwhere
      n=sum(mask)
      allocate(w(n))
      k=1
      do  i=1,size(a)
        if (mask(i)==1) then
          w(k)=i
          k=k+1
        endif
      enddo
  return
  end function findmin_D

  function findmin_R (a,w,n) result(amin)
 !   a(:) is 1D vector
 !   w(:) 1D vector including indice of minimum elements
 !   n is the number of minimum elements in a(:)
      implicit none
      real (kind=4), intent(in), dimension (:) :: a
      integer, intent(out), allocatable :: w(:)
      integer :: n, mask(size(a)), i,k
      real (kind=4) :: amin

      amin=minval(a)
      where (a == amin)
         mask=1
      elsewhere
         mask=0
      endwhere
      n=sum(mask)
      allocate(w(n))
      k=1
      do  i=1,size(a)
        if (mask(i)==1) then
          w(k)=i
          k=k+1
        endif
      enddo
  return
  end function findmin_R

  function findmin_I (a,w,n) result(amin)
 !   a(:) is 1D vector
 !   w(:) 1D vector including indice of minimum elements
 !   n is the number of minimum elements in a(:)
      implicit none
      integer, intent(in), dimension (:) :: a
      integer, intent(out), allocatable :: w(:)
      integer :: n, mask(size(a)), i,k
      integer :: amin

      amin=minval(a)
      where (a == amin)
         mask=1
      elsewhere
         mask=0
      endwhere
      n=sum(mask)
      allocate(w(n))
      k=1
      do  i=1,size(a)
        if (mask(i)==1) then
          w(k)=i
          k=k+1
        endif
      enddo
  return
  end function findmin_I


END MODULE InternPro


MODULE TomoPro
  USE InternPro, only: double
  IMPLICIT NONE

!INTEGER, PARAMETER :: double=SELECTED_REAL_KIND(15)
REAL (kind=double), PARAMETER :: Rsun=6.96d10 !cm  solar radius
REAL (kind=double), PARAMETER :: grid_Rout=4.0_double
REAL (kind=double), PARAMETER :: grid_Rin=1.5_double
INTEGER, PARAMETER :: grid_sz=128

TYPE SparseMat  !Compressed Sparse Row (CSR) for m-row x n-col Matrix
 !ifst is allocated to the number of rows m plus 1, ie. m+1
 !icol and A are allocated to asiz
 !the entries in ifst[1:m] must be in range [1,size], ifst[m+1] should always be size+1
 !the entries in icol must be in range of column numbers [1,n]
  INTEGER :: m,n  ! number of rwos and cols
  INTEGER :: asiz ! number of non-zero elements
  INTEGER, ALLOCATABLE :: ifst(:)  ! index in A(:) of first non-0 element for given row
  INTEGER, ALLOCATABLE :: icol(:)  ! column number of non-0 elements stored
  REAL, ALLOCATABLE :: A(:)   ! element value
END TYPE SparseMat 

TYPE dataRys
   INTEGER :: n
   INTEGER, ALLOCATABLE :: ja(:)
   REAL, ALLOCATABLE :: aa(:)
END TYPE dataRys

TYPE Mapstr
      REAL(kind=8) :: dx,dy,xc,yc,rsun
      REAL, DIMENSION(:,:), ALLOCATABLE :: data
END TYPE Mapstr


  INTERFACE thomsonscatter
    MODULE PROCEDURE thomsonscatter_1d, thomsonscatter_scl
  END INTERFACE

  INTERFACE index_to_ijk
     MODULE PROCEDURE index_to_ijk_1d, index_to_ijk_scl
  END INTERFACE

  INTERFACE coord2pix
    MODULE PROCEDURE  coord2pix_1d,  coord2pix_scl
  END INTERFACE

  INTERFACE getmap_coord
    MODULE PROCEDURE getmap_coord_scl, getmap_coord_1d, getmap_coord_2d
  END INTERFACE

  CONTAINS


SUBROUTINE Allocate_Rys (row, n)
  IMPLICIT NONE
  TYPE (dataRys), INTENT(INOUT) :: row
  INTEGER, INTENT(IN) :: n
  INTEGER :: AllocateStatus
  IF (ALLOCATED(row%aa)) THEN
     IF (row%n /= n) THEN
         call DeAllocate_Rys (row)
     ENDIF   
  ENDIF
  IF (row%n /= n) THEN
     row%n = n
     ALLOCATE(row%ja(n), row%aa(n), STAT=AllocateStatus)
     IF (AllocateStatus /= 0) STOP " ***Not enough MEM for allocating type (dataRys) ***" 
  ENDIF
  RETURN
END SUBROUTINE Allocate_Rys

SUBROUTINE DeAllocate_Rys (row)
  IMPLICIT NONE
  TYPE (dataRys), INTENT(INOUT) :: row
  DEALLOCATE(row%ja, row%aa)
  row%n = 0
  RETURN
END SUBROUTINE DeAllocate_Rys

FUNCTION thomsonscatter_1d (Rin, Tht, limb, is_rad) RESULT(pB)
!+calculate Thomson Scatter light
!Input:
! Rin - impact distance in unit of Rsun 
!    (nearest distance from Sun center to LOS)
! Tht - angle in degree of the radial direction through an electron
!     with respect to the plane of sky
!keyword: limb=0.63 default (520nm in green), limb darkening coef.
! output in unit of mean disk Brightness
! similar to the routine ELTHEORY in SSW
!  rad  -if set, Theta in radian, otherwise in degree
! Return: polarized brightness (pB) from single electron
!-
 USE InternPro, only : radeg
 IMPLICIT NONE

  REAL, INTENT(IN), DIMENSION(:) :: Rin, Tht
  REAL, INTENT(IN) :: limb
  LOGICAL, INTENT(IN) :: is_rad
  REAL(kind=double), DIMENSION(size(Rin)) :: pB
  REAL (kind=double), PARAMETER :: sigmat=6.65245D-25, Rs=1.0_double 
 !sigmat is Thomson scatter cross-section. Rs is Solar radius in unit of Rsun
  REAL (kind=double), DIMENSION(size(Rin)) :: theta,sinchi,sinw,cosw,A,B
  REAL (kind=double) :: u,Icen

 IF (is_rad) THEN
    theta=dble(Tht)
 ELSE
    theta=dble(Tht)/radeg
 ENDIF

 IF (limb <0) THEN
    u=0.63d0 ! set the default value if limb is not given
 ELSE
    u=limb
 ENDIF
 Icen=1.D0/(1.d0-u/3.d0) !Sun center brightness in unit of mean brightness
 sinchi=cos(theta)
! R=Rin/sinchi
 sinw=Rs*sinchi/Rin
 cosw=sqrt(1.d0-sinw**2)  !w <90. deg
!calculate Billing's coef.
 A=cosw*sinw**2
 B=-1.d0/8.d0*(1.d0-3.d0*sinw**2-cosw**2/sinw*(1.d0 &
      +3.d0*sinw**2)*log((1.d0+sinw)/cosw))
!C=4./3.-cosw-cosw^3/3.
!D=1./8.*(5.+sinw^2-cosw^2/sinw*(5.-sinw^2)*alog((1.+sinw)/cosw))
!polarized brightness pB =Bt - Br
pB=3.d0/16.d0*Icen*sigmat*sinchi**2*((1.d0-u)*A+u*B)
RETURN
END FUNCTION thomsonscatter_1d 

FUNCTION thomsonscatter_scl (Rin, Tht, limb, is_rad) RESULT(pB)
! input and output arguments are scalar
 USE InternPro, only : radeg
 IMPLICIT NONE
  REAL, INTENT(IN) :: Rin, Tht
  REAL, INTENT(IN) :: limb
  LOGICAL, INTENT(IN) :: is_rad
  REAL (kind=double) :: pB
  REAL (kind=double), PARAMETER :: sigmat=6.65245D-25, Rs=1.0_double
 !sigmat is Thomson scatter cross-section. Rs is Solar radius in unit of Rsun
  REAL (kind=double) :: theta,sinchi,sinw,cosw,A,B
  REAL (kind=double) :: u,Icen

 IF (is_rad) THEN
    theta=dble(Tht)
 ELSE
    theta=dble(Tht)/radeg
 ENDIF

IF (limb <0) THEN
    u=0.63d0 ! set the default value if limb is not given
 ELSE
    u=limb
 ENDIF
 Icen=1.D0/(1.d0-u/3.d0) !Sun center brightness in unit of mean brightness
 sinchi=cos(theta)
! R=Rin/sinchi
 sinw=Rs*sinchi/Rin
 cosw=sqrt(1.d0-sinw**2)  !w <90. deg
!calculate Billing's coef.
 A=cosw*sinw**2
 B=-1.d0/8.d0*(1.d0-3.d0*sinw**2-cosw**2/sinw*(1.d0 &
      +3.d0*sinw**2)*log((1.d0+sinw)/cosw))
!C=4./3.-cosw-cosw^3/3.
!D=1./8.*(5.+sinw^2-cosw^2/sinw*(5.-sinw^2)*alog((1.+sinw)/cosw))
!polarized brightness pB =Bt - Br
pB=3.d0/16.d0*Icen*sigmat*sinchi**2*((1.d0-u)*A+u*B)
RETURN
END FUNCTION thomsonscatter_scl


FUNCTION coord2pix_1d (x,y,z) RESULT (f)
!+ transform coord. [x,y,z] of grid points into pixel coord.
! x,y,z-  1d vector
! return - [3,n] where n is number of total pixels
IMPLICIT NONE
REAL (kind=double), INTENT(IN), DIMENSION(:) :: x,y,z
REAL (kind=double), DIMENSION(3,SIZE(x)) :: f
REAL (kind=double) :: ds
 !define a fixed 3d grid:
 ds=2.0_double*grid_Rout/dble(grid_sz-1) !size of the cell in the domain
 f(1,:)=(x+grid_Rout)/ds
 f(2,:)=(y+grid_Rout)/ds
 f(3,:)=(z+grid_Rout)/ds
RETURN
END FUNCTION coord2pix_1d

FUNCTION coord2pix_scl (x,y,z) RESULT (f)
!+ transform coord. [x,y,z] of one grid point into pixel coord.
! x,y,z- scalar
! return -  a vector f[3] including pixel coord (xp,yp,zp)
! where xp in [0: M-1], yp in [0, N-1], zp in [0, K-1] assuming
! th grid size as MxNxK 
IMPLICIT NONE
REAL (kind=double), INTENT(IN) :: x,y,z
REAL (kind=double), DIMENSION(3) :: f
REAL (kind=double) :: ds
 !define a fixed 3d grid:
 ds=2.0_double*grid_Rout/dble(grid_sz-1) !size of the cell in the domain
 f(1)=(x+grid_Rout)/ds
 f(2)=(y+grid_Rout)/ds
 f(3)=(z+grid_Rout)/ds
RETURN
END FUNCTION coord2pix_scl

! read map-str from an unformatted data file created by IDL
FUNCTION read_mapstr(filename) RESULT(map)
  IMPLICIT NONE
  CHARACTER(LEN=50), INTENT(IN) :: filename
  TYPE (Mapstr) :: map
  INTEGER :: isz
  OPEN(1, FILE=filename, FORM="UNFORMATTED")
  READ(1) isz
  ALLOCATE(map%data(isz,isz))
  READ(1) map%xc, map%yc, map%dx, map%dy, map%rsun
  READ(1) map%data
  CLOSE(1)
RETURN
END FUNCTION read_mapstr

SUBROUTINE getmap_coord_scl(map,xp,yp,xcoord,ycoord)
  IMPLICIT NONE
  TYPE (Mapstr), INTENT(IN) :: map
  REAL (kind=double), INTENT(IN) :: xp,yp
  REAL (kind=double), INTENT(OUT) :: xcoord,ycoord
  REAL :: xdim, ydim
  xdim=size(map%data,1)
  ydim=size(map%data,2)
  xcoord=(xp-(xdim+1.d0)/2.d0)*map%dx +map%xc
  ycoord=(yp-(ydim+1.d0)/2.d0)*map%dy +map%yc
RETURN
END SUBROUTINE getmap_coord_scl

SUBROUTINE getmap_coord_1d(map,xp,yp,xcoord,ycoord)
  IMPLICIT NONE
  TYPE (Mapstr), INTENT(IN) :: map
  REAL (kind=double), INTENT(IN), DIMENSION(:) :: xp,yp
  REAL (kind=double), INTENT(OUT),DIMENSION(:) :: xcoord,ycoord
  REAL :: xdim, ydim
  xdim=size(map%data,1)
  ydim=size(map%data,2)
  xcoord=(xp-(xdim+1.d0)/2.d0)*map%dx +map%xc
  ycoord=(yp-(ydim+1.d0)/2.d0)*map%dy +map%yc
RETURN
END SUBROUTINE getmap_coord_1d

SUBROUTINE getmap_coord_2d(map,xp,yp,xcoord,ycoord)
  IMPLICIT NONE
  TYPE (Mapstr), INTENT(IN) :: map
  REAL(kind=double), INTENT(IN), DIMENSION(:,:) :: xp,yp
  REAL (kind=double), INTENT(OUT),DIMENSION(:,:) :: xcoord,ycoord
  REAL :: xdim, ydim
  xdim=size(map%data,1)
  ydim=size(map%data,2)
  xcoord=(xp-(xdim+1.d0)/2.d0)*map%dx +map%xc
  ycoord=(yp-(ydim+1.d0)/2.d0)*map%dy +map%yc
RETURN
END SUBROUTINE getmap_coord_2d


SUBROUTINE index_to_ijk_1d(index,iw,jw,kw)
 ! convert 1D-index of array(I,J,K) in range=[0:I*J*K-1] into 3 indice [i,j,k]
 IMPLICIT NONE
 INTEGER, INTENT(IN), DIMENSION(:) :: index
 INTEGER, INTENT(OUT), DIMENSION(size(index)) :: iw,jw,kw
 INTEGER :: tmp, i
   DO i=1, size(index)
      kw(i)=index(i)/grid_sz**2+1
      tmp=mod(index(i),grid_sz**2)
      jw(i)=tmp/grid_sz+1
      iw(i)=mod(tmp,grid_sz)+1
   ENDDO
END SUBROUTINE index_to_ijk_1d

SUBROUTINE index_to_ijk_scl(index,iw,jw,kw)
! convert one index of array(I,J,K) with value in range=[0:I*J*K-1]
!   into 3 indice [i,j,k] where i=[1,I], j=[1,J], k=[1,K]
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: index
 INTEGER, INTENT(OUT) :: iw,jw,kw
 INTEGER :: tmp
      kw=index/grid_sz**2+1
      tmp=mod(index,grid_sz**2)
      jw=tmp/grid_sz+1
      iw=mod(tmp,grid_sz)+1
END SUBROUTINE index_to_ijk_scl

FUNCTION ijk_to_index (i,j,k) RESULT(index)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: i,j,k
  INTEGER :: index
  index = (i-1) + (j-1)*grid_sz + (k-1)*grid_sz**2
RETURN
END FUNCTION ijk_to_index


FUNCTION point_weight (x,y,z,dx,dy,dz) RESULT(fwt)
!+ purpose: calculate the weight (line integral) of density for
! segment of a ray in a unit cell with N=1 at the vertice
! (0,0,0) and N=0 for other vertices.
! Input: (x,y,z) entry point;
!       (dx,dy,dz) vector from entry point to exit point
IMPLICIT NONE
REAL (kind=double), INTENT(IN) :: x,y,z,dx,dy,dz
REAL (kind=double) :: d,fwt,f1,f2,f3,f4

 d=sqrt(dx*dx+dy*dy+dz*dz)
 f1= d*(1.0-x-y-z +x*y+x*z+y*z -x*y*z)
 f2= d/2.0*(-dx-dy-dz +dx*z+dx*y+dy*x+dy*z+dz*x+dz*y -dx*y*z-dy*x*z-dz*x*y)
 f3= d/3.0*(dx*dy+dx*dz+dy*dz -dx*dy*z-dx*dz*y-dy*dz*x)
 f4=-d/4.0*dx*dy*dz
 fwt=f1+f2+f3+f4
 return
END FUNCTION point_weight


FUNCTION cell_weigths (index,xt,yt,zt) RESULT(wt)
!+purpose: return weights calculated for a cell with given index (i,j,k)
! input:
!  index[3] - pixel index of a given cell, corresponding to [ix,jy,kz]
!        where ix=int(xp), jy=int(yp), kz=int(zp) where (xp,yp,zp) is
!        pixel coord., so ix in [0, N-1], jy in [0, M-1], kz in [0, K-1]
!        for grid size NxMxK
!  xt[2],yt[2],zt[2] - entry[i=1] and exit[i=2] coords. for a ray crossing the cell
! output:
!   return weights wt[2,2,2] for 8 vertices with i=0,1;j=0,1;k=0,1
!-
IMPLICIT NONE
INTEGER, INTENT(IN), DIMENSION(3) :: index
REAL (kind=double), INTENT(IN), DIMENSION(2) :: xt,yt,zt
REAL (kind=double), DIMENSION(2,2,2) :: wt
REAL (kind=double), DIMENSION(3,2) :: f
REAL (kind=double) :: x0,y0,z0,x1,y1,z1,x,y,z,dx,dy,dz
INTEGER :: i,j,k
 
  f=coord2pix(xt,yt,zt)

  x0=f(1,1)-real(index(1))
  y0=f(2,1)-real(index(2))
  z0=f(3,1)-real(index(3)) !entry point coord in unit cell

  x1=f(1,2)-real(index(1))
  y1=f(2,2)-real(index(2))
  z1=f(3,2)-real(index(3)) !exit point coord in unit cell

  DO 100 k=0,1 
  DO 100 j=0,1 
  DO 100 i=0,1
     x=(1.0-i)*x0 +i*(1.0-x0)
     y=(1.0-j)*y0 +j*(1.0-y0)
     z=(1.0-k)*z0 +k*(1.0-z0)
     dx=(-1)**i*(x1-x0)
     dy=(-1)**j*(y1-y0)
     dz=(-1)**k*(z1-z0)
   wt(i+1,j+1,k+1)=point_weight(x,y,z,dx,dy,dz)
  100 CONTINUE
return
END FUNCTION cell_weigths


FUNCTION ray_weights (xe,ye,ze,ccell) RESULT(wt)
!+purpose: calculate weights (contribution func to line integral) for a ray
!passing through a grid box
! input:
!  xe,ye,ze[N] - Carrington coords. of N crossing points for a ray
! ccell[3,N-1] - indice of N-1 cells which are passed through
! return: weights wt[SZ,SZ,SZ] for grid contributed by the ray
!   where SZ=grid_size
!-
IMPLICIT NONE
REAL (kind=double), INTENT(IN), DIMENSION(:) :: xe,ye,ze
INTEGER, INTENT(IN), DIMENSION(3,size(xe)-1) :: ccell
REAL (kind=double), DIMENSION(0:grid_sz-1,0:grid_sz-1,0:grid_sz-1) :: wt
INTEGER :: i,j,k,i0,j0,k0, m, index(3)
REAL (kind=double), DIMENSION(2,2,2) :: f
REAL (kind=double), DIMENSION(2) :: xt,yt,zt


! n=n_elements(xe)-1 ; total number of crossed cells
! wt=dblarr(sz,sz,sz)

  wt=0. 
 !note: initial values must be set because wt=wt+f is used!
 DO m=1, size(xe)-1 
    index=ccell(:,m)
    xt=xe(m:m+1)
    yt=ye(m:m+1)
    zt=ze(m:m+1)
    f=cell_weigths(index,xt,yt,zt)
    DO 100 k=0,1 
    DO 100 j=0,1 
    DO 100 i=0,1 
        i0=index(1)+i
        j0=index(2)+j
        k0=index(3)+k
        wt(i0,j0,k0)=wt(i0,j0,k0)+f(i+1,j+1,k+1)
    100 CONTINUE
 ENDDO

return
END FUNCTION ray_weights


SUBROUTINE get_cross_points (cobs,xi,yi,zi,rdis,is_spb, xe_cr,ye_cr,ze_cr,ccell_cr, &
              cross, ni, xt,yt,zt)
!;+ for a given point with image coord. [xi,yi,zi] and an observer with Carrington coord
!; cobs[2]=[lon, lat], determine the entry point of LOS (defined from the object to
!; observer) at the boundary of grid box
!; input: cobs[2] -observer's position in degrees
!;        xi,yi,zi - image coord. of a point
!;        rdis[n,n,n]- radial distance of grids
!; output:
!;  xe[N],ye[N],ze[N] - Carrington coords. of N cross-points
!;  ccell[3,N-1] - index for N-1 crossing cells
!:  xt[2],yt[2],zt[2] - coord. of extry amd exit points
!; keyword: cross =1, return entry point; cross=0, not crossing.
!;          if cross =1, calculate Carrington coord of ientry and exit points
!;              at the grid boundaries
!;       ci[3] - if set, using a point in Carringon coord as a input to
!;            replace [xi,yi,zi] in image coord.
!;      ni[3] - direction vector of the LOS
!;v1. improve the agorithm for calculation of entry point
!!;    add keyword ci[3] for test point
 USE InternPro
 IMPLICIT NONE

 INTEGER, PARAMETER :: maxnum=grid_sz**2
 REAL (kind=double), INTENT(IN), DIMENSION(2) :: cobs
 REAL (kind=double), INTENT(IN) :: xi,yi,zi
 REAL, INTENT(IN), DIMENSION(:,:,:) :: rdis
 LOGICAL, INTENT(IN) :: is_spb !if true, set spheric boundary for grid box
 REAL (kind=double), INTENT(OUT), DIMENSION(:), allocatable :: xe_cr,ye_cr,ze_cr
 INTEGER, INTENT(OUT), DIMENSION(:,:), allocatable :: ccell_cr
 INTEGER, INTENT(OUT) :: cross
 REAL (kind=double), INTENT(OUT), DIMENSION(3) :: ni
 REAL, INTENT(OUT), DIMENSION(2) :: xt,yt,zt
 REAL (kind=double) :: ds
 INTEGER :: num, i0,j0,k0, AllocateStatus
 REAL (kind=double) , DIMENSION(maxnum) :: xe,ye,ze
 INTEGER, DIMENSION(3,maxnum) :: ccell
 REAL (kind=double), DIMENSION(3) :: dli,dle,dlni,dlne
 REAL (kind=double) :: phi,tht, disi_max,dise_min 
 REAL (kind=double) :: x,y,z, x1,y1,z1, x2,y2,z2, xbi,ybi,zbi,xbe,ybe,zbe
 REAL (kind=double) :: tlen,r12, xe0,ye0,ze0,xe1,ye1,ze1, f(3)
 LOGICAL :: flag_inbox 
 INTEGER, DIMENSION(1) :: wi
 INTEGER :: i,j,k, ix0,iy0,iz0, ixv(8),iyv(8),izv(8), nvt, wpb(maxnum),npb, nmw
 REAL (kind=double) :: fmx
 INTEGER, ALLOCATABLE :: windx(:)
 INTEGER :: nfm
 flag_inbox=.TRUE.  
!note: initial value of flag_inbox must be set here. 
!   If its value is set during declaration, it becomes a static variable. 
!   It will keep its value as called in the last time
 cross=1
 ds=2.d0*grid_Rout/dble(grid_sz-1)

!the boundaries of grid
 x1=-grid_Rout 
 y1=-grid_Rout 
 z1=-grid_Rout 
 x2=grid_Rout
 y2=grid_Rout
 z2=grid_Rout

 !calculate the direction vector of LOS
 phi=cobs(1)/radeg
 tht=cobs(2)/radeg
 ni=(/cos(phi)*cos(tht), sin(phi)*cos(tht), sin(tht)/)
 IF (ni(1)==0 .or. ni(2)==0 .or. ni(3)==0) STOP '* One of direction-COS is 0 in SUB get_cross_points()'
 ! ;transform image coord of the point to Carrington coord.
  call img2car_coord (cobs, xi, yi, zi, x, y, z)
 
 !determine the entry and exit points of LOS-ray with the grid box
   if (ni(1) > 0) then 
      xbi=x1
      xbe=x2 
   else 
      xbi=x2
      xbe=x1
   endif
   if (ni(2) > 0) then
      ybi=y1
      ybe=y2
   else
      ybi=y2
      ybe=y1
   endif
   if (ni(3) > 0) then
      zbi=z1
      zbe=z2
   else
      zbi=z2
      zbe=z1
   endif

 dli=(/xbi-x,ybi-y,zbi-z/)
 dlni=dli/ni
 disi_max=maxval(dlni)

 dle=(/xbe-x,ybe-y,zbe-z/)
 dlne=dle/ni
 dise_min=minval(dlne)

  IF (dise_min < disi_max) THEN ! no crossing between ray and grid box
     cross=0
     print *,'Not crossing!'
     print *,'Carr.coord. of given point: [x,y,z]=',x,y,z,' [nx,ny,nz]=',ni
     xt=0.
     yt=0.
     zt=0.
     RETURN
  ELSE
!    print *,'11111 --the given ray crossing grid-box!'
     cross=1
     xt(1)=x + disi_max*ni(1) ! for entry point
     yt(1)=y + disi_max*ni(2)
     zt(1)=z + disi_max*ni(3)
     xt(2)=x + dise_min*ni(1) ! for exit point
     yt(2)=y + dise_min*ni(2)
     zt(2)=z + dise_min*ni(3)
  ENDIF

! set as the first entry point
 xe0=xt(1)
 ye0=yt(1)
 ze0=zt(1)

 tlen=0.  !;calculate total length of ray in grid box for test
 num=1 !;order of cross-cell
 xe(num)=xe0
 ye(num)=ye0
 ze(num)=ze0

!print*,'flag_inbox=',flag_inbox,' num=',num
 DO WHILE (flag_inbox .and. num<maxnum) 
   if (num ==1) then !for entry point
       wi=maxloc(dlni)
       fmx=findmax(dlni,windx,nfm)
   else if (num >1) then !for exit points
       wi=minloc(dlne)
       fmx=findmin(dlne,windx,nfm)
   endif 
   if (nfm >1) then
       print*,'number of max or min elements=',nfm
       STOP "---Error due to a ray crossing side or vertice of grid cell! Special treatment is required!"
   endif

   f=coord2pix(xe0,ye0,ze0)
   IF (wi(1) == 1) THEN !cross x-side x=xbi
      if (ni(1)>0) then  
         xbe=xe0+ds
         i0=nint(f(1))
      else
         xbe=xe0-ds
         i0=nint(f(1))-1
      endif
      j0=int(f(2))
      k0=int(f(3))
      if (ni(2)>0) then 
         ybe=y1+(j0+1)*ds 
      else 
         ybe=y1+j0*ds
      endif
      if (ni(3)>0) then 
         zbe=z1+(k0+1)*ds 
      else 
         zbe=z1+k0*ds
      endif
   ENDIF
       !cross y-side, y=ybi
   IF (wi(1) == 2) THEN 
       if (ni(2)>0)  then 
          ybe=ye0+ds 
          j0=nint(f(2))
       else 
          ybe=ye0-ds
          j0=nint(f(2))-1
       endif
       i0=int(f(1))
       k0=int(f(3))
      if (ni(1)>0) then 
         xbe=x1+(i0+1)*ds 
      else 
         xbe=x1+i0*ds
      endif
      if (ni(3)>0) then 
         zbe=z1+(k0+1)*ds 
      else 
         zbe=z1+k0*ds
      endif
   ENDIF 
   !cross z-side,z=zbi
   IF (wi(1) == 3) THEN
       if (ni(3)>0)  then
          zbe=ze0+ds
          k0=nint(f(3))
       else
          zbe=ze0-ds
          k0=nint(f(3))-1
       endif
       i0=int(f(1))
       j0=int(f(2))
      if (ni(1)>0) then
         xbe=x1+(i0+1)*ds
      else
         xbe=x1+i0*ds
      endif
      if (ni(2)>0) then
         ybe=y1+(j0+1)*ds
      else
         ybe=y1+j0*ds
      endif
   ENDIF
     !determine index of the cell to which an entry point belongs
      ccell(1,num)=i0
      ccell(2,num)=j0
      ccell(3,num)=k0
!print *, 'num=',num, ' ',xe0,ye0,ze0,ccell(:,num)
 !calculate exit point for current cell 
   dle=(/xbe-xe(1),ybe-ye(1),zbe-ze(1)/)
   dlne=dle/ni  !dlne will be used to find the next w(1)
   dise_min=minval(dlne) !distance of the ray in cell
   xe1=xe(1)+dise_min*ni(1)
   ye1=ye(1)+dise_min*ni(2)
   ze1=ze(1)+dise_min*ni(3)
! print,'exit point of cell(m)=',xe1,ye1,ze1
   r12=sqrt((xe1-xe0)**2+(ye1-ye0)**2+(ze1-ze0)**2)
   tlen=tlen+r12

  IF ((xe1 < x1).or.(xe1 > x2).or.(ye1 < y1).or.(ye1 > y2) &
      .or.(ze1 < z1).or.(ze1 > z2)) THEN
     flag_inbox=.FALSE.  ! the exit point is outside the grid box
     tlen=tlen-r12
!    print *,'--the exit point is out of grid box!---1'
!    print *,'total cross-point num=',num
!    print *,'total segment length by suming num-1 cells=',tlen
  ELSE 
     xe(num+1)=xe1
     ye(num+1)=ye1
     ze(num+1)=ze1
     xe0=xe1 ! set exit point of current cell as entry point for next cell
     ye0=ye1
     ze0=ze1
     num=num+1
  ENDIF
 END DO

 IF (num == maxnum .or. flag_inbox) THEN
     print*,' flag_inbox=',flag_inbox
!    print*,'current exit point=',xe1,ye1,ze1
!     print*,'entry point=',xt(1),yt(1),zt(1)
!     print*,'exit point=',xt(2),yt(2),zt(2)
     STOP "loop number for getting crossing points reaches Maximum--Error!"
 ENDIF
!print *,'total cross-point num=',num

 ! for spheric boundary if is_spb=true
 set_spb: IF (is_spb) THEN
!    print *,'---- trim crossed points and cells within Spheric Domain!!!!'
     k=1
   loop_spb: DO i=1, num-1 
     ix0=ccell(1,i)
     iy0=ccell(2,i)
     iz0=ccell(3,i)
     ixv=(/ix0,ix0+1,ix0+1,ix0, ix0,ix0+1,ix0+1,ix0/) !indice of 8 vertices for cell-i
     iyv=(/iy0,iy0,iy0+1,iy0+1, iy0,iy0,iy0+1,iy0+1/)
     izv=(/iz0,iz0,iz0,iz0, iz0+1,iz0+1,iz0+1,iz0+1/)
     nvt=0
     DO j=1,8
        if (rdis(ixv(j)+1,iyv(j)+1,izv(j)+1) > grid_Rout) nvt=nvt+1
     ENDDO
 !boundary cell taken as the one with at least
 !one vertice meeting r>Rout and at least one vertice meeating r<=Rout
     IF (nvt >0 .and. nvt<8) THEN
        wpb(k)=i
        k=k+1
     ENDIF
   END DO loop_spb
   npb=k-1
   !  print *,'boundary cell num=',npb
   if (npb ==0) then
      print*,'entry point=',xt(1),yt(1),zt(1) 
      print*,'exit point=',xt(2),yt(2),zt(2)
      STOP "there is an error in finding boundary cells!"
   endif
   ! select inbox-cells as the one with at least one vertice within grid-box
    nmw=wpb(npb)-wpb(1)+2
    ccell(:,1:nmw-1)=ccell(:,wpb(1):wpb(npb))
    xe(1:nmw)=xe(wpb(1):wpb(npb)+1)
    ye(1:nmw)=ye(wpb(1):wpb(npb)+1)
    ze(1:nmw)=ze(wpb(1):wpb(npb)+1)

  ENDIF set_spb

  allocate(xe_cr(nmw),ye_cr(nmw),ze_cr(nmw),ccell_cr(3,nmw-1),STAT = AllocateStatus)
  IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
  xe_cr=xe(1:nmw)
  ye_cr=ye(1:nmw)
  ze_cr=ze(1:nmw)
  ccell_cr=ccell(:,1:nmw-1)
END SUBROUTINE get_cross_points


SUBROUTINE Allocate_SP (A,m,n,asiz)
IMPLICIT NONE
TYPE(SparseMat), INTENT(inout) :: A
INTEGER, INTENT(in) :: m,n,asiz
INTEGER :: AllocateStatus
 IF (allocated(A%A)) THEN
    IF (A%n/=n .or. A%m/=m .or. A%asiz/=asiz) THEN
       call DeAllocate_SP(A)
    ENDIF
 ENDIF
 IF (A%n/=n .or. A%m/=m .or. A%asiz/=asiz) THEN
    allocate(A%ifst(m+1))
    allocate(A%icol(asiz), A%A(asiz), STAT=AllocateStatus)
    IF (AllocateStatus /= 0) STOP " ***Not enough MEM for allocating type (SparseMat) ***"
    A%n=n; A%m=m; A%asiz=asiz
 ENDIF
RETURN
END SUBROUTINE Allocate_SP

SUBROUTINE DeAllocate_SP(A)
  IMPLICIT NONE
  TYPE(SparseMat), INTENT(inout) :: A
  deallocate(A%ifst, A%icol, A%A)
  A%n=0; A%m=0; A%asiz=0
RETURN
END SUBROUTINE DeAllocate_SP

FUNCTION Sparse2Mat (A) RESULT(Mat)
!+ purpose: transform the compressed sparse Matrix into normal format
! Input:
! A - struct, CRS format matrix
! Return: normal m x n Matrix, M(n,m)
  IMPLICIT NONE
  TYPE(SparseMat), INTENT(in) :: A
  REAL, DIMENSION(A%m,A%n) :: Mat
  INTEGER :: i,j,istt,iend
  Mat=0. ! initial values must be given
  DO i=1,A%m
    istt=A%ifst(i)
    iend=A%ifst(i+1)-1
    DO j=istt,iend
       Mat(i,A%icol(j))=A%A(j)
    ENDDO
  ENDDO
 RETURN
END FUNCTION Sparse2Mat

FUNCTION Mat2Sparse (Mat) RESULT(B)
!  transform the normal-form Matrix into the CRS sparseMatrix
  IMPLICIT NONE
  REAL, INTENT(IN), DIMENSION(:,:) :: Mat !m-row x n-col
  TYPE(SparseMat) :: A,B
  INTEGER :: m,n, i,j,k, inz, ifst, asiz

  m=size(Mat,1); n=size(Mat,2)
  call Allocate_SP(A, m,n,m*n)
  ifst=1
  k=1 !index of non-zero value in A%A
  DO i=1,m
     inz=0  ! number of non-zero elements in i-row
     DO j=1,n
        IF (abs(Mat(i,j)) > 0) THEN
           inz=inz+1
           A%icol(k)=j
           A%A(k)=Mat(i,j)
           k=k+1
        ENDIF
     ENDDO
     A%ifst(i)=ifst
     ifst=ifst + inz
  ENDDO
  asiz= k-1  !total number of non-zero elements
  A%ifst(m+1)=ifst !This also gives total number of non-zero elements plus 1, ie.= asiz+1
  IF (A%ifst(m+1)== asiz+1) THEN !check error
     print*,'total non-zero elements number=',asiz
  ELSE
     STOP "A%ifst(m+1) not = asiz --ERROR!"
  ENDIF
  call Allocate_SP(B, m,n,asiz)
  B%ifst=A%ifst
  B%icol=A%icol(1:asiz)
  B%A   =A%A(1:asiz)
  print *,'Compressed Sparse A is created!'
  call DeAllocate_SP(A)
RETURN
END FUNCTION Mat2Sparse

FUNCTION Index_SP(A,i,j)
! returns the element index (in SP%A) of element (i,j) of the SparseMAt matrix A
! If element (i,j) cannot be found (ie.is a zero-element), set index=9
IMPLICIT NONE
INTEGER :: Index_SP
TYPE (SparseMat), INTENT(IN) :: A
INTEGER, INTENT(IN) :: i,j
INTEGER :: ii
  ii=A%ifst(i)
  DO WHILE(A%icol(ii)/=j.AND.ii<A%ifst(i+1)-1)
     ii=ii+1
  ENDDO
  IF (A%icol(ii)/=j) THEN
    Index_SP=0
  ELSE
    Index_SP=ii
  ENDIF
RETURN
END FUNCTION Index_SP

SUBROUTINE Aij_SP(A,index, i,j)
!return the index (i,j) of a non-zero element with indexi_sp (in SP%A) for
!Sparse matrix A
IMPLICIT NONE
TYPE (SparseMat), INTENT(IN) :: A
INTEGER, INTENT(IN) :: index
INTEGER, INTENT(OUT):: i,j
INTEGER :: ii
  IF (Index>A%asiz) THEN
      i=0; j=0
      print*,'ERROR! index_sp is out of index range for SP%A!'
  ENDIF
  DO ii=1,A%m
     IF (index >=A%ifst(ii).AND.index<A%ifst(ii+1)) THEN
        i=ii
        EXIT
     ENDIF
  ENDDO
  j=A%icol(index)
END SUBROUTINE Aij_SP


FUNCTION RMatMul_SP (A,x,n) RESULT(y)
!+purpose: right multiply SparseMat A with vector x: RMatMul_SP=A*x
!Input:
!  A-struct -Compressed Row Sparse (CRS) format Matrix
! x(n) - Vector
!  n  - number of elements in x
! Return:  vector RMatMul_SP=y
! if x(n) is shorter than A.n, x is implicitly zero filled
! if x(n) is longer than A.n, x is truncated to A.n
IMPLICIT NONE
TYPE (SparseMat), INTENT(IN) :: A
INTEGER, INTENT(IN) :: n
REAL, INTENT(IN), DIMENSION(n) :: x
REAL, DIMENSION(A%m) :: y
INTEGER :: irow,i, istt,iend

DO irow=1, A%m
   y(irow)=0.
   istt=A%ifst(irow)
   iend=A%ifst(irow+1)-1
   DO i=istt, iend
      IF (A%icol(i)<=n) THEN
         y(irow)=y(irow)+A%A(i)*x(A%icol(i))
      ENDIF
   ENDDO
ENDDO
RETURN
END FUNCTION RMatMul_SP

FUNCTION LMatMul_SP (x,m,A) RESULT(y)
!+purpose: left multiply SparseMat A with vector x: LMatMul_SP=x*A=A^T*x
!Input:
!  A-struct -Compressed Row Sparse (CRS) format Matrix
! x(m) - Vector
!  m  - number of elements in x
! Return:  vector RMatMul_SP=y
! if x(m) is shorter than A.m, x is implicitly zero filled
! if x(m) is longer than A.m,  x is truncated to A.m
IMPLICIT NONE
INTEGER, INTENT(IN) :: m
REAL, INTENT(IN), DIMENSION(m):: x
TYPE (SparseMat), INTENT(IN) :: A
REAL, DIMENSION(A%n) :: y
INTEGER :: i,j, istt, iend, irow

y=0.
irow=min(A%m,m)
DO i=1, irow
   istt=A%ifst(i)
   iend=A%ifst(i+1)-1
   DO j=istt, iend
      y(A%icol(j))=y(A%icol(j)) + x(i)*A%A(j)
   ENDDO
ENDDO
END FUNCTION LMatMul_SP


FUNCTION ConjugateGrad_SP (A,ydata,mu,R, tol,sav_res,sav_iter,resume)  RESULT(xsolu)
!+purpose, solving the regularized minimization,
!    x_mu=argmin(|A*x-y|^2+mu*|R*x|^2 by solving
! M*x=b, where M=A^T*A + mu*R^T*R and b=A^T*y
! using the conjugate gradient method, where matrix M is a
! real, symmetric, positive-definite matrix
! Input:
! A: struct, m-row X n-col compressed row sparse (CRS) matrix for los integral
! y(m)   -a vector including observed data
! R: struct, 3n-row X n-col CRS matrix for regularization
! mu: scalar, regularization factor. If mu=0, no regularization
!    tol: the desired convergence tolerance
!       define tol=|A*x-y|/|y|
!     if sav_res=1: save the middle results for resume in a file [resume_iter.dat]. 
!     if sav_iter=1: save data and regularization errors for each step if se
!     if resum=1: to resume the calculation (e.g. for tol set to be smaller value)
!-
IMPLICIT NONE
INTEGER, PARAMETER ::  itmax =20000 
TYPE (SparseMat), INTENT(IN):: A, R
REAL, INTENT(IN), DIMENSION(A%m) :: ydata
REAL,INTENT(IN) :: mu, tol
INTEGER, INTENT(IN) :: sav_res, sav_iter, resume
REAL, DIMENSION(A%n) :: bvec,xk,xk1,rk,rk1,pk,pk1,Apk,Rpk, xsolu, xknn
REAL :: Ap(A%m),Rp(R%m),btol, rrk,rrk1,alphak,betak,rchi
REAL :: resA(A%m),resR(R%m),Data_err,Reg_err, itol,rerr 
INTEGER :: k
LOGICAL :: exist

 print*,'Solving linear system by conjugate method -----'
 IF (mu <= 0.) print*,'---Regularization not included!'
 IF (tol <=0.) THEN
     itol=1.d-3  !default value of iteration tolerance
 ELSE
     itol=tol
 ENDIF
 bvec=LMatMul_SP(ydata,size(ydata),A)
 btol=sqrt(sum(bvec**2)) 
 xk=0. !set initial solution as zero
 rk=bvec !r0=b-A*x0, set initial residual
 pk=rk
 k=0
 rrk=dot_product(rk,rk)
 ! resume iterating calculation by reading the saved middle parameters
 !   from [resume_iter.dat]
 IF (resume ==1) THEN
    OPEN(2, FILE='resume_iter.dat',FORM='unformatted',ACTION='read')
    READ(2) k,rrk
    READ(2) rk,pk,xk
    CLOSE(2) 
    print*,'Resume the calculation by conjugate gradient---'
 ENDIF

 DO WHILE (k < itmax)
    Ap=RMatMul_SP(A,pk,A%n)  ! Ap[m] =A*pk
    Apk=LMatMul_SP(Ap,A%m,A)  !Apk[n]= A^T*Ap=Ap*A
    IF (mu > 0) THEN   !include regularization
       Rp=RMatMul_SP(R,pk,A%n)   !Rp[m1]=R*pk where m1=R.m
       Rpk=LMatMul_SP(Rp,R%m,R)  !Rpk[n]=R^T*Rp
       Apk=Apk+mu*Rpk            !Apk=A^T*A*pk + mu*R^T*R*pk
    ENDIF
    alphak=rrk/dot_product(pk,Apk)
    xk1=xk +alphak*pk
    rk1=rk -alphak*Apk
    rerr=sqrt(sum(rk1**2))/btol  !relative error
  print*,'k=',k,' res=',rerr, ' (tol=',tol,')'
    IF (rerr <= tol) EXIT
    rrk1=dot_product(rk1,rk1)
    betak=rrk1/rrk
    pk1=rk1 +betak*pk
    k=k + 1
    rk=rk1
    pk=pk1
    xk=xk1
    rrk=rrk1
    IF (sav_res==1) THEN !save the middle result for resume
        OPEN(1,FILE='resume_iter.dat',FORM='unformatted',ACTION='write',STATUS='replace')
        WRITE(1) k,rrk
        WRITE(1) rk,pk,xk
        CLOSE(1)
    ENDIF
    IF (sav_iter==1) THEN !save data and reg. errors for each iteration step
    !  calculate data and regularization errors
       WHERE(xk1 >0) !set xknn including non-negative elements in xk1(:) 
          xknn=xk1
       ELSEWHERE
          xknn=0.
       END WHERE
       resA=RMatMul_SP(A,xknn,A%n)-ydata
       Data_err=sum(resA**2)
       resR=RMatMul_SP(R,xknn,A%n)
       Reg_err=sum(resR**2)
       print*,'data err=||Ax-y||^2=',Data_err,' reg err=||Rx||^2=',Reg_err

       INQUIRE(FILE='iter_errs.tex', exist=exist)
          if (.not. exist) then
               OPEN(3, FILE='iter_errs.tex', ACTION='write',STATUS='new')
          else
               OPEN(3, FILE='iter_errs.tex', ACTION='write',POSITION='append',STATUS='old')
          endif
       WRITE(3,*) k,Data_err,Reg_err, rerr
       CLOSE(3)
  ENDIF
 ENDDO
  xsolu=xk1
RETURN
END FUNCTION ConjugateGrad_SP


FUNCTION Reg2ndsm2cdf_SP (wdis, rdis, nr) RESULT(R)
!+ purpose: Second-order smoothing using central differetiate of 2nd-order ccuracy
! Input: n - dimension of grid
!       wdis -index for grid points with r>Rin
!       rdis - distance of grids to Sun center for grid [n,n,n]
! return: R -mm x nn sparse matrix for regularization, where
!          nn=ndis and mm=3*nn
! add weight=1/nr for each grid point remove r-dependence effect of regularization
USE InternPro
IMPLICIT NONE
INTEGER, INTENT(IN), DIMENSION(:) :: wdis
REAL, INTENT(IN), DIMENSION(:,:,:) :: rdis
REAL (kind=double), INTENT(IN), DIMENSION(:,:) :: nr
TYPE (SparseMat) :: R
TYPE (dataRys) :: Nx, Ny, Nz, mRys(3*size(wdis))
REAL (kind=double) :: Rin, Rout, ds
INTEGER :: nn, mm, nnz, i,j,k, iw1,iw2,iw3, ix2(2),ix3(3)
INTEGER, DIMENSION(size(wdis)) :: iwd, iw,jw,kw
REAL (kind=double), DIMENSION(size(wdis)) :: xnr, ynr
INTEGER, DIMENSION(3*size(wdis)+1) :: IA
INTEGER, ALLOCATABLE :: state(:)
INTEGER :: im, ii, nja

print*,'---create regularized R_sprs using Second-order differetiate, |d^2f/dx^2|=|(f(i+1)-2f(i)+f(i-1))/dx^2|'

ds=2.0_double*grid_Rout/dble(grid_sz-1)
Rin=grid_Rin - 2*ds
Rout=grid_Rout + 2.d0*ds
nn=size(wdis)
mm=3*nn
iwd=(/(i,i=1,nn)/)  !column number of R_sparse matrix

call Allocate_Rys(Nx, 3)
call Allocate_Rys(Ny, 3)
call Allocate_Rys(Nz, 3)

IF (nr(1,1)>0) THEN
   print*,'====> Radial weight Nr is applied to R_sparse!'
   call index_to_ijk(wdis,iw,jw,kw)
   DO i=1,nn
      xnr(i)=rdis(iw(i),jw(i),kw(i))
   ENDDO
   call interp(nr(:,1),nr(:,2),xnr,ynr)
ENDIF
  ii=1  ! row number of R_sparse matrix
  IA(1)=1  ! index in R%aa for 1st non-0 element in 1st row of R_sparse
row_loop: DO im=1, nn
   call index_to_ijk(wdis(im),i,j,k)
 !for df/dx=N(i-1)-2*N(i)+N(i+1)
    IF (i==1) THEN
       Nx%n =2 
       iw1=ijk_to_index(i,j,k)
       iw2=ijk_to_index(i+1,j,k)
       ix2=(/iw1,iw2/)
       Nx%ja(1:2)= looktable(wdis,iwd,ix2, state)
       Nx%aa(1:2)= (/-1, 1/)
    ELSE IF (i==grid_sz) THEN
       Nx%n =2
       iw1=ijk_to_index(i-1,j,k)
       iw2=ijk_to_index(i,j,k)
       ix2=(/iw1,iw2/)
       Nx%ja(1:2)= looktable(wdis,iwd,ix2, state)
       Nx%aa(1:2)= (/-1, 1/)
    ELSE
       Nx%n =3 
       iw1=ijk_to_index(i-1,j,k)
       iw2=ijk_to_index(i,j,k)
       iw3=ijk_to_index(i+1,j,k)
       if ((rdis(i-1,j,k)<Rin .and. rdis(i,j,k)>Rin) .or. &
           (rdis(i-1,j,k)>Rout .and. rdis(i,j,k)<Rout)) then
          Nx%n =2
          iw1=ijk_to_index(i,j,k)
          iw2=ijk_to_index(i+1,j,k)
          ix2=(/iw1,iw2/)
          Nx%ja(1:2)= looktable(wdis,iwd,ix2, state)
          Nx%aa(1:2)= (/-1, 1/)
       else if ((rdis(i+1,j,k)<Rin .and. rdis(i,j,k)>Rin) .or. &
                (rdis(i+1,j,k)>Rout .and. rdis(i,j,k)<Rout)) then
          Nx%n =2
          iw1=ijk_to_index(i-1,j,k)
          iw2=ijk_to_index(i,j,k)
          ix2=(/iw1,iw2/)
          Nx%ja(1:2)= looktable(wdis,iwd,ix2, state)
          Nx%aa(1:2)= (/-1, 1/)
       else
          ix3=(/iw1,iw2,iw3/)
          Nx%ja(:)= looktable(wdis,iwd,ix3, state)
          Nx%aa(:)= (/1, -2, 1/)
       endif
    ENDIF

 !for df/dy=N(j-1)-2*N(j)+N(j+1)
    IF (j==1) THEN
       Ny%n =2
       iw1=ijk_to_index(i,j,k)
       iw2=ijk_to_index(i,j+1,k)
       ix2=(/iw1,iw2/)
       Ny%ja(1:2)= looktable(wdis,iwd,ix2, state)
       Ny%aa(1:2)= (/-1, 1/)
    ELSE IF (j==grid_sz) THEN
       Ny%n =2
       iw1=ijk_to_index(i,j-1,k)
       iw2=ijk_to_index(i,j,k)
       ix2=(/iw1,iw2/)
       Ny%ja(1:2)= looktable(wdis,iwd,ix2, state)
       Ny%aa(1:2)= (/-1, 1/)
    ELSE
       Ny%n =3
       iw1=ijk_to_index(i,j-1,k)
       iw2=ijk_to_index(i,j,k)
       iw3=ijk_to_index(i,j+1,k)
       if ((rdis(i,j-1,k)<Rin .and. rdis(i,j,k)>Rin) .or. &
           (rdis(i,j-1,k)>Rout .and. rdis(i,j,k)<Rout)) then
          Ny%n =2
          iw1=ijk_to_index(i,j,k)
          iw2=ijk_to_index(i,j+1,k)
          ix2=(/iw1,iw2/)
          Ny%ja(1:2)= looktable(wdis,iwd,ix2, state)
          Ny%aa(1:2)= (/-1, 1/)
       else if ((rdis(i,j+1,k)<Rin .and. rdis(i,j,k)>Rin) .or. &
                (rdis(i,j+1,k)>Rout .and. rdis(i,j,k)<Rout)) then
          Ny%n =2
          iw1=ijk_to_index(i,j-1,k)
          iw2=ijk_to_index(i,j,k)
          ix2=(/iw1,iw2/)
          Ny%ja(1:2)= looktable(wdis,iwd,ix2, state)
          Ny%aa(1:2)= (/-1, 1/)
       else
          ix3=(/iw1,iw2,iw3/)
          Ny%ja(:)= looktable(wdis,iwd,ix3, state)
          Ny%aa(:)= (/1, -2, 1/)
       endif
    ENDIF

 !for df/dz=N(k-1)-2*N(k)+N(k+1)
    IF (k==1) THEN
       Nz%n =2
       iw1=ijk_to_index(i,j,k)
       iw2=ijk_to_index(i,j,k+1)
       ix2=(/iw1,iw2/)
       Nz%ja(1:2)= looktable(wdis,iwd,ix2, state)
       Nz%aa(1:2)= (/-1, 1/)
    ELSE IF (k==grid_sz) THEN
       Nz%n =2
       iw1=ijk_to_index(i,j,k-1)
       iw2=ijk_to_index(i,j,k)
       ix2=(/iw1,iw2/)
       Nz%ja(1:2)= looktable(wdis,iwd,ix2, state)
       Nz%aa(1:2)= (/-1, 1/)
    ELSE
       Nz%n =3
       iw1=ijk_to_index(i,j,k-1)
       iw2=ijk_to_index(i,j,k)
       iw3=ijk_to_index(i,j,k+1)
       if ((rdis(i,j,k-1)<Rin .and. rdis(i,j,k)>Rin) .or. &
           (rdis(i,j,k-1)>Rout .and. rdis(i,j,k)<Rout)) then
          Nz%n =2
          iw1=ijk_to_index(i,j,k)
          iw2=ijk_to_index(i,j,k+1)
          ix2=(/iw1,iw2/)
          Nz%ja(1:2)= looktable(wdis,iwd,ix2, state)
          Nz%aa(1:2)= (/-1, 1/)
       else if ((rdis(i,j,k+1)<Rin .and. rdis(i,j,k)>Rin) .or. &
                (rdis(i,j,k+1)>Rout .and. rdis(i,j,k)<Rout)) then
          Nz%n =2
          iw1=ijk_to_index(i,j,k-1)
          iw2=ijk_to_index(i,j,k)
          ix2=(/iw1,iw2/)
          Nz%ja(1:2)= looktable(wdis,iwd,ix2, state)
          Nz%aa(1:2)= (/-1, 1/)
       else
          ix3=(/iw1,iw2,iw3/)
          Nz%ja(:)= looktable(wdis,iwd,ix3, state)
          Nz%aa(:)= (/1, -2, 1/)
       endif
    ENDIF

    IF (nr(1,1) >0) THEN !add the radial weights
       Nx%aa=Nx%aa/ynr(im)
       Ny%aa=Ny%aa/ynr(im)
       Nz%aa=Nz%aa/ynr(im)
    ENDIF

    nja=Nx%n
    call Allocate_Rys(mRys(ii),nja)
    mRys(ii)%ja = Nx%ja(1:nja)
    mRys(ii)%aa = Nx%aa(1:nja)
    IA(ii+1) = IA(ii)+nja
    ii=ii+1

    nja=Ny%n
    call Allocate_Rys(mRys(ii),nja)
    mRys(ii)%ja = Ny%ja(1:nja)
    mRys(ii)%aa = Ny%aa(1:nja)
    IA(ii+1) = IA(ii)+nja
    ii=ii+1

    nja=Nz%n
    call Allocate_Rys(mRys(ii),nja)
    mRys(ii)%ja = Nz%ja(1:nja)
    mRys(ii)%aa = Nz%aa(1:nja)
    IA(ii+1) = IA(ii)+nja
    ii=ii+1
ENDDO row_loop
! print*,'ii=',ii,' mm=',mm, '  Should be: ii=mm+1'

 nnz = IA(mm+1)-1
 call Allocate_SP(R, mm, nn, nnz)
 R%ifst = IA
 DO i=1,mm
    R%icol(IA(i):IA(i+1)-1) = mRys(i)%ja
    R%A(IA(i):IA(i+1)-1) = mRys(i)%aa
    call DeAllocate_Rys(mRys(i))
 ENDDO
 call DeAllocate_Rys(Nx)
 call DeAllocate_Rys(Ny)
 call DeAllocate_Rys(Nz)
RETURN
END FUNCTION Reg2ndsm2cdf_SP


FUNCTION Reg2ndsm4cdf_SP (wdis, rdis, nr) RESULT(R)
!+ purpose: Second-order smoothing using central differetiate of 4nd-order ccuracy
! Input: n - dimension of grid
!       wdis -index for grid points with r>Rin
!       rdis - distance of grids to Sun center for grid [n,n,n]
! return: R -mm x nn sparse matrix for regularization, where
!          nn=ndis and mm=3*nn
! add weight=1/nr for each grid point remove r-dependence effect of
! regularization
USE InternPro
IMPLICIT NONE
INTEGER, INTENT(IN), DIMENSION(:) :: wdis
REAL, INTENT(IN), DIMENSION(:,:,:) :: rdis
REAL (kind=double), INTENT(IN), DIMENSION(:,:) :: nr
TYPE (SparseMat) :: R
TYPE (dataRys) :: Nx, Ny, Nz, mRys(3*size(wdis))
REAL (kind=double) :: Rin, Rout, ds
INTEGER :: nn, mm, nnz, i,j,k, iw1,iw2,iw3,iw4,iw5, ix2(2),ix5(5)
INTEGER, DIMENSION(size(wdis)) :: iwd, iw,jw,kw
REAL (kind=double), DIMENSION(size(wdis)) :: xnr, ynr
INTEGER, DIMENSION(3*size(wdis)+1) :: IA
INTEGER, ALLOCATABLE :: state(:)
INTEGER :: im, ii, nja
REAL :: c1,c2,c3,c4
integer :: i0,j0,k0

 print*,'---create regularized R_sprs using Second-order differetiate,&
     &|d^2f/dx^2|=|(-f(i+2)+16*f(i+1)-30*f(i)+16*f(i-1)-f(i-2))/12*dx^2|'

  ds=2.0_double*grid_Rout/dble(grid_sz-1)
  Rin=grid_Rin - 2*ds
  Rout=grid_Rout + 2.d0*ds
  nn=size(wdis)
  mm=3*nn
  iwd=(/(i,i=1,nn)/)  !column number of R_sparse matrix

 c1=-1.
 c2= 16.
 c3=-30.
 c4=12.  

call Allocate_Rys(Nx, 5)
call Allocate_Rys(Ny, 5)
call Allocate_Rys(Nz, 5)

IF (nr(1,1)>0) THEN
   print*,'====> Radial weight Nr is applied to R_sparse!'
   call index_to_ijk(wdis,iw,jw,kw)
   DO i=1,nn
      xnr(i)=rdis(iw(i),jw(i),kw(i))
   ENDDO
   call interp(nr(:,1),nr(:,2),xnr,ynr)
ENDIF

  ii=1  ! row number of R_sparse matrix
  IA(1)=1  ! index in R%aa for 1st non-0 element in 1st row of R_sparse
row_loop: DO im=1, nn
   call index_to_ijk(wdis(im),i,j,k)
! 2nd-order smoothing
!for d2f/dx2=c1*N(i+2)+c2*N(i+1)+c3*N(i)+c2*N(i-1)+c1*N(i-2)
    IF (i==1) THEN  !right diff.= c4*Rdf where Rfg=f(i+1)-f(I)
       Nx%n =2
       iw1=ijk_to_index(i,j,k)
       iw2=ijk_to_index(i+1,j,k)
       ix2=(/iw1,iw2/)
       Nx%ja(1:2)= looktable(wdis,iwd,ix2, state)
       Nx%aa(1:2)= (/-1, 1/)*c4
    ELSE IF (i==2) THEN  !central diff.= c4*(cdf/2) where cdf=f(i+1)-f(i-1)
       Nx%n =2
       iw1=ijk_to_index(i-1,j,k)
       iw2=ijk_to_index(i+1,j,k)
       if (rdis(i-1,j,k)>Rout) then
           iw1=ijk_to_index(i,j,k) !right diff
           Nx%aa(1:2)= (/-1, 1/)*c4
       else
           Nx%aa(1:2)= (/-1, 1/)*c4/2.
       endif
       ix2=(/iw1,iw2/)
       Nx%ja(1:2)= looktable(wdis,iwd,ix2, state)
    ELSE IF (i==grid_sz) THEN !left diff.
       Nx%n =2
       iw1=ijk_to_index(i-1,j,k)
       iw2=ijk_to_index(i,j,k)
       ix2=(/iw1,iw2/)
       Nx%ja(1:2)= looktable(wdis,iwd,ix2, state)
       Nx%aa(1:2)= (/-1, 1/)*c4
    ELSE IF (i==grid_sz-1) THEN !central diff.= c4*(cdf/2)
       Nx%n =2
       iw1=ijk_to_index(i-1,j,k)
       iw2=ijk_to_index(i+1,j,k)
       if (rdis(i+1,j,k)>Rout) then
           iw2=ijk_to_index(i,j,k) !left diff
           Nx%aa(1:2)= (/-1, 1/)*c4
       else
           Nx%aa(1:2)= (/-1, 1/)*c4/2.
       endif
       ix2=(/iw1,iw2/)
       Nx%ja(1:2)= looktable(wdis,iwd,ix2, state)
    ELSE
       Nx%n =5
       iw1=ijk_to_index(i-2,j,k)
       iw2=ijk_to_index(i-1,j,k)
       iw3=ijk_to_index(i,j,k)
       iw4=ijk_to_index(i+1,j,k)
       iw5=ijk_to_index(i+2,j,k)
       if ((rdis(i-2,j,k)<Rin.and.rdis(i-1,j,k)>Rin) .OR.&  !central diff
           (rdis(i+2,j,k)<Rin.and.rdis(i+1,j,k)>Rin) .OR.&
           (rdis(i-2,j,k)>Rout.and.rdis(i-1,j,k)<Rout) .OR.&
           (rdis(i+2,j,k)>Rout.and.rdis(i+1,j,k)<Rout)) then
           Nx%n =2
           iw1=ijk_to_index(i-1,j,k)
           iw2=ijk_to_index(i+1,j,k)
           ix2=(/iw1,iw2/)
           Nx%ja(1:2)= looktable(wdis,iwd,ix2, state)
           if (state(0) > 0) then
              print*,'---abandon this peculiar point with rdis=',rdis(i,j,k)
              cycle
           endif 
           Nx%aa(1:2)= (/-1, 1/)*c4/2.
       else if ((rdis(i-1,j,k)<Rin.and.rdis(i,j,k)>Rin) .OR.&  !right diff. 
                (rdis(i-1,j,k)>Rout.and.rdis(i,j,k)<Rout)) then
           Nx%n =2
           iw1=ijk_to_index(i,j,k)
           iw2=ijk_to_index(i+1,j,k)
           ix2=(/iw1,iw2/)
           Nx%ja(1:2)= looktable(wdis,iwd,ix2, state)
           Nx%aa(1:2)= (/-1, 1/)*c4
        else if ((rdis(i+1,j,k)<Rin.and.rdis(i,j,k)>Rin) .OR.&  !left diff.
                (rdis(i+1,j,k)>Rout.and.rdis(i,j,k)<Rout)) then
           Nx%n =2
           iw1=ijk_to_index(i-1,j,k)
           iw2=ijk_to_index(i,j,k)
           ix2=(/iw1,iw2/)
           Nx%ja(1:2)= looktable(wdis,iwd,ix2, state)
           Nx%aa(1:2)= (/-1, 1/)*c4 
        else
           ix5=(/iw1,iw2,iw3,iw4,iw5/)
           Nx%ja(:)= looktable(wdis,iwd,ix5, state)
           Nx%aa(:)= (/c1,c2,c3,c2,c1/)
           if (state(0) > 0) then
              print*,'---abandon this peculiar point with rdis=',rdis(i,j,k)
              cycle
           endif
        endif
    ENDIF 

!  ;for d2f/dy2=c1*N(j+2)+c2*N(j+1)+c3*N(j)+c2*N(j-1)+c1*N(j-2)
    IF (j==1) THEN  !right diff.= c4*Rdf where Rfg=f(i+1)-f(I)
       Ny%n =2
       iw1=ijk_to_index(i,j,k)
       iw2=ijk_to_index(i,j+1,k)
       ix2=(/iw1,iw2/)
       Ny%ja(1:2)= looktable(wdis,iwd,ix2, state)
       Ny%aa(1:2)= (/-1, 1/)*c4
    ELSE IF (j==2) THEN  !central diff.= c4*(cdf/2) where cdf=f(i+1)-f(i-1)
       Ny%n =2
       iw1=ijk_to_index(i,j-1,k)
       iw2=ijk_to_index(i,j+1,k)
       if (rdis(i,j-1,k)>Rout) then
           iw1=ijk_to_index(i,j,k) !right diff
           Ny%aa(1:2)= (/-1, 1/)*c4
       else
           Ny%aa(1:2)= (/-1, 1/)*c4/2.
       endif
       ix2=(/iw1,iw2/)
       Ny%ja(1:2)= looktable(wdis,iwd,ix2, state)
    ELSE IF (j==grid_sz) THEN !left diff.
       Ny%n =2
       iw1=ijk_to_index(i,j-1,k)
       iw2=ijk_to_index(i,j,k)
       ix2=(/iw1,iw2/)
       Ny%ja(1:2)= looktable(wdis,iwd,ix2, state)
       Ny%aa(1:2)= (/-1, 1/)*c4
    ELSE IF (j==grid_sz-1) THEN !central diff.= c4*(cdf/2)
       Ny%n =2
       iw1=ijk_to_index(i,j-1,k)
       iw2=ijk_to_index(i,j+1,k)
       if (rdis(i,j+1,k)>Rout) then
           iw2=ijk_to_index(i,j,k) !left diff
           Ny%aa(1:2)= (/-1, 1/)*c4
       else
           Ny%aa(1:2)= (/-1, 1/)*c4/2.
       endif
       ix2=(/iw1,iw2/)
       Ny%ja(1:2)= looktable(wdis,iwd,ix2, state)
    ELSE
       Ny%n=5
       iw1=ijk_to_index(i,j-2,k)
       iw2=ijk_to_index(i,j-1,k)
       iw3=ijk_to_index(i,j,k)
       iw4=ijk_to_index(i,j+1,k)
       iw5=ijk_to_index(i,j+2,k)
       if ((rdis(i,j-2,k)<Rin.and.rdis(i,j-1,k)>Rin) .OR.&  !central diff
           (rdis(i,j+2,k)<Rin.and.rdis(i,j+1,k)>Rin) .OR.&
           (rdis(i,j-2,k)>Rout.and.rdis(i,j-1,k)<Rout) .OR.&
           (rdis(i,j+2,k)>Rout.and.rdis(i,j+1,k)<Rout)) then
           Ny%n =2
           iw1=ijk_to_index(i,j-1,k)
           iw2=ijk_to_index(i,j+1,k)
           ix2=(/iw1,iw2/)
           Ny%ja(1:2)= looktable(wdis,iwd,ix2, state)
           if (state(0) > 0) then
              print*,'---abandon this peculiar point with rdis=',rdis(i,j,k)
              cycle
           endif
           Ny%aa(1:2)= (/-1, 1/)*c4/2.
       else if ((rdis(i,j-1,k)<Rin.and.rdis(i,j,k)>Rin) .OR.&  !right diff.
                (rdis(i,j-1,k)>Rout.and.rdis(i,j,k)<Rout)) then
           Ny%n =2
           iw1=ijk_to_index(i,j,k)
           iw2=ijk_to_index(i,j+1,k)
           ix2=(/iw1,iw2/)
           Ny%ja(1:2)= looktable(wdis,iwd,ix2, state)
           Ny%aa(1:2)= (/-1, 1/)*c4
        else if ((rdis(i,j+1,k)<Rin.and.rdis(i,j,k)>Rin) .OR.&  !left diff.
                (rdis(i,j+1,k)>Rout.and.rdis(i,j,k)<Rout)) then
           Ny%n =2
           iw1=ijk_to_index(i,j-1,k)
           iw2=ijk_to_index(i,j,k)
           ix2=(/iw1,iw2/)
           Ny%ja(1:2)= looktable(wdis,iwd,ix2, state)
           Ny%aa(1:2)= (/-1, 1/)*c4
        else
           ix5=(/iw1,iw2,iw3,iw4,iw5/)
           Ny%ja(:)= looktable(wdis,iwd,ix5, state)
           Ny%aa(:)= (/c1,c2,c3,c2,c1/)
           if (state(0) > 0) then
              print*,'---abandon this peculiar point with rdis=',rdis(i,j,k)
              cycle
           endif
        endif
    ENDIF

!;for d2f/dz2=c1*N(k+2)+c2*N(k+1)+c3*N(k)+c2*N(k-1)+c1*N(k-2)
    IF (k==1) THEN  !right diff.= c4*Rdf where Rfg=f(i+1)-f(I)
       Nz%n =2
       iw1=ijk_to_index(i,j,k)
       iw2=ijk_to_index(i,j,k+1)
       ix2=(/iw1,iw2/)
       Nz%ja(1:2)= looktable(wdis,iwd,ix2, state)
       Nz%aa(1:2)= (/-1, 1/)*c4
    ELSE IF (k==2) THEN  !central diff.= c4*(cdf/2) where cdf=f(i+1)-f(i-1)
       Nz%n =2
       iw1=ijk_to_index(i,j,k-1)
       iw2=ijk_to_index(i,j,k+1)
       if (rdis(i,j,k-1)>Rout) then
           iw1=ijk_to_index(i,j,k) !right diff
           Nz%aa(1:2)= (/-1, 1/)*c4
       else
           Nz%aa(1:2)= (/-1, 1/)*c4/2.
       endif
       ix2=(/iw1,iw2/)
       Nz%ja(1:2)= looktable(wdis,iwd,ix2, state)
    ELSE IF (k==grid_sz) THEN !left diff.
       Nz%n =2
       iw1=ijk_to_index(i,j,k-1)
       iw2=ijk_to_index(i,j,k)
       ix2=(/iw1,iw2/)
       Nz%ja(1:2)= looktable(wdis,iwd,ix2, state)
       Nz%aa(1:2)= (/-1, 1/)*c4
    ELSE IF (k==grid_sz-1) THEN !central diff.= c4*(cdf/2)
       Nz%n =2
       iw1=ijk_to_index(i,j,k-1)
       iw2=ijk_to_index(i,j,k+1)
       if (rdis(i,j,k+1)>Rout) then
           iw2=ijk_to_index(i,j,k) !left diff
           Nz%aa(1:2)= (/-1, 1/)*c4
       else
           Nz%aa(1:2)= (/-1, 1/)*c4/2.
       endif
       ix2=(/iw1,iw2/)
       Nz%ja(1:2)= looktable(wdis,iwd,ix2, state)
    ELSE
       Nz%n=5
       iw1=ijk_to_index(i,j,k-2)
       iw2=ijk_to_index(i,j,k-1)
       iw3=ijk_to_index(i,j,k)
       iw4=ijk_to_index(i,j,k+1)
       iw5=ijk_to_index(i,j,k+2)
       if ((rdis(i,j,k-2)<Rin.and.rdis(i,j,k-1)>Rin) .OR.&  !central diff
           (rdis(i,j,k+2)<Rin.and.rdis(i,j,k+1)>Rin) .OR.&
           (rdis(i,j,k-2)>Rout.and.rdis(i,j,k-1)<Rout) .OR.&
           (rdis(i,j,k+2)>Rout.and.rdis(i,j,k+1)<Rout)) then
           Nz%n =2
           iw1=ijk_to_index(i,j,k-1)
           iw2=ijk_to_index(i,j,k+1)
           ix2=(/iw1,iw2/)
           Nz%ja(1:2)= looktable(wdis,iwd,ix2, state)
           if (state(0) > 0) then
              print*,'---abandon this peculiar point with rdis=',rdis(i,j,k)
              cycle
           endif
           Nz%aa(1:2)= (/-1, 1/)*c4/2.
       else if ((rdis(i,j,k-1)<Rin.and.rdis(i,j,k)>Rin) .OR.&  !right diff.
                (rdis(i,j,k-1)>Rout.and.rdis(i,j,k)<Rout)) then
           Nz%n =2
           iw1=ijk_to_index(i,j,k)
           iw2=ijk_to_index(i,j,k+1)
           ix2=(/iw1,iw2/)
           Nz%ja(1:2)= looktable(wdis,iwd,ix2, state)
           Nz%aa(1:2)= (/-1, 1/)*c4
        else if ((rdis(i,j,k+1)<Rin.and.rdis(i,j,k)>Rin) .OR.&  !left diff.
                (rdis(i,j,k+1)>Rout.and.rdis(i,j,k)<Rout)) then
           Nz%n =2
           iw1=ijk_to_index(i,j,k-1)
           iw2=ijk_to_index(i,j,k)
           ix2=(/iw1,iw2/)
           Nz%ja(1:2)= looktable(wdis,iwd,ix2, state)
           Nz%aa(1:2)= (/-1, 1/)*c4
        else
           ix5=(/iw1,iw2,iw3,iw4,iw5/)
           Nz%ja(:)= looktable(wdis,iwd,ix5, state)
           Nz%aa(:)= (/c1,c2,c3,c2,c1/)
           if (state(0) > 0) then
              print*,'---abandon this peculiar point with rdis=',rdis(i,j,k)
              cycle
           endif
        endif
    ENDIF

    IF (nr(1,1) >0) THEN !add the radial weights
       Nx%aa=Nx%aa/ynr(im)
       Ny%aa=Ny%aa/ynr(im)
       Nz%aa=Nz%aa/ynr(im)
    ENDIF

    nja=Nx%n
    call Allocate_Rys(mRys(ii),nja)
    mRys(ii)%ja = Nx%ja(1:nja)
    mRys(ii)%aa = Nx%aa(1:nja)
    IA(ii+1) = IA(ii)+nja
    ii=ii+1

    nja=Ny%n
    call Allocate_Rys(mRys(ii),nja)
    mRys(ii)%ja = Ny%ja(1:nja)
    mRys(ii)%aa = Ny%aa(1:nja)
    IA(ii+1) = IA(ii)+nja
    ii=ii+1

    nja=Nz%n
    call Allocate_Rys(mRys(ii),nja)
    mRys(ii)%ja = Nz%ja(1:nja)
    mRys(ii)%aa = Nz%aa(1:nja)
    IA(ii+1) = IA(ii)+nja
    ii=ii+1
 ENDDO row_loop
print*,'ii=',ii,' initial mm=',mm, ' actual mm=',ii-1

 mm=ii-1
 nnz = IA(mm+1)-1
 call Allocate_SP(R, mm, nn, nnz)
 R%ifst = IA
 DO i=1,mm
    R%icol(IA(i):IA(i+1)-1) = mRys(i)%ja
    R%A(IA(i):IA(i+1)-1) = mRys(i)%aa
    call DeAllocate_Rys(mRys(i))
 ENDDO
 call DeAllocate_Rys(Nx)
 call DeAllocate_Rys(Ny)
 call DeAllocate_Rys(Nz)
RETURN
END FUNCTION Reg2ndsm4cdf_SP

FUNCTION Reg1stsm2cdf_SP (wdis, rdis, nr) RESULT(R)
!+ purpose: First-order smoothing using central differetiate of 2nd-order acccuracy
!dfdx(i)=N(i+1)-N(i-1)
! Input: n - dimension of grid
!       wdis -index for grid points with r>Rin
!       rdis - distance of grids to Sun center for grid [n,n,n]
! return: R -mm x nn sparse matrix for regularization, where
!          nn=ndis and mm=3*nn
! add weight=1/nr for each grid point remove r-dependence effect of
! regularization
USE InternPro
IMPLICIT NONE
INTEGER, INTENT(IN), DIMENSION(:) :: wdis
REAL, INTENT(IN), DIMENSION(:,:,:) :: rdis
REAL (kind=double), INTENT(IN), DIMENSION(:,:) :: nr
TYPE (SparseMat) :: R
TYPE (dataRys) :: Nx, Ny, Nz, mRys(3*size(wdis))
REAL (kind=double) :: Rin, Rout, ds
INTEGER :: nn, mm, nnz, i,j,k, iw1,iw2,iw3, ix2(2)
INTEGER, DIMENSION(size(wdis)) :: iwd, iw,jw,kw
REAL (kind=double), DIMENSION(size(wdis)) :: xnr, ynr
INTEGER, DIMENSION(3*size(wdis)+1) :: IA
INTEGER, ALLOCATABLE :: state(:)
INTEGER :: im, ii, nja

print*,'---Create regularized R_sprs using First-order central differetiate, |df/dx|^2=|(f(i+1)-f(i-1))/2dx|^2'
ds=2.0_double*grid_Rout/dble(grid_sz-1)
Rin=grid_Rin - 2*ds
Rout=grid_Rout + 2.d0*ds
nn=size(wdis)
mm=3*nn
iwd=(/(i,i=1,nn)/)  !column number of R_sparse matrix

call Allocate_Rys(Nx, 2)
call Allocate_Rys(Ny, 2)
call Allocate_Rys(Nz, 2)

IF (nr(1,1)>0) THEN
   print*,'====> Radial weight Nr is applied to R_sparse!'
   call index_to_ijk(wdis,iw,jw,kw)
   DO i=1,nn
      xnr(i)=rdis(iw(i),jw(i),kw(i))
   ENDDO
   call interp(nr(:,1),nr(:,2),xnr,ynr)
ENDIF

  ii=1  ! row number of R_sparse matrix
  IA(1)=1  ! index in R%aa for 1st non-0 element in 1st row of R_sparse
row_loop: DO im=1, nn
   call index_to_ijk(wdis(im),i,j,k)
 !;for df/dx=N(i+1)-N(i-1)
    IF (i==1) THEN    !;right diff
       Nx%n =2
       iw1=ijk_to_index(i,j,k)
       iw2=ijk_to_index(i+1,j,k)
       ix2=(/iw1,iw2/)
       Nx%ja(1:2)= looktable(wdis,iwd,ix2, state)
       Nx%aa(1:2)= (/-1, 1/)*2.
    ELSE IF (i==grid_sz) THEN
       Nx%n =2    !left diff
       iw1=ijk_to_index(i-1,j,k)
       iw2=ijk_to_index(i,j,k)
       ix2=(/iw1,iw2/)
       Nx%ja(1:2)= looktable(wdis,iwd,ix2, state)
       Nx%aa(1:2)= (/-1, 1/)*2.
    ELSE
       Nx%n =2  !central diff
       iw1=ijk_to_index(i-1,j,k)
       iw2=ijk_to_index(i+1,j,k)
       if (rdis(i-1,j,k)<Rin .or. rdis(i-1,j,k)>Rout) then
             !right diff
          iw1=ijk_to_index(i,j,k)
          ix2=(/iw1,iw2/)
          Nx%ja(1:2)= looktable(wdis,iwd,ix2, state)
          Nx%aa(1:2)= (/-1, 1/)*2.
       else if (rdis(i+1,j,k)<Rin .or. rdis(i+1,j,k)>Rout) then
            !left diff
          iw2=ijk_to_index(i,j,k)
          ix2=(/iw1,iw2/)
          Nx%ja(1:2)= looktable(wdis,iwd,ix2, state)
          Nx%aa(1:2)= (/-1, 1/)*2.
       else   !central diff
          ix2=(/iw1,iw2/)
          Nx%ja(:)= looktable(wdis,iwd,ix2, state)
          Nx%aa(:)= (/-1., 1./)
       endif
    ENDIF

   !for dfdy=N(j+1)-N(j-1)
    IF (j==1) THEN    !;right diff
       Ny%n =2
       iw1=ijk_to_index(i,j,k)
       iw2=ijk_to_index(i,j+1,k)
       ix2=(/iw1,iw2/)
       Ny%ja(1:2)= looktable(wdis,iwd,ix2, state)
       Ny%aa(1:2)= (/-1, 1/)*2.
    ELSE IF (j==grid_sz) THEN
       Ny%n =2    !left diff
       iw1=ijk_to_index(i,j-1,k)
       iw2=ijk_to_index(i,j,k)
       ix2=(/iw1,iw2/)
       Ny%ja(1:2)= looktable(wdis,iwd,ix2, state)
       Ny%aa(1:2)= (/-1, 1/)*2.
    ELSE
       Ny%n =2  !central diff
       iw1=ijk_to_index(i,j-1,k)
       iw2=ijk_to_index(i,j+1,k)
       if (rdis(i,j-1,k)<Rin .or. rdis(i,j-1,k)>Rout) then
             !right diff
          iw1=ijk_to_index(i,j,k)
          ix2=(/iw1,iw2/)
          Ny%ja(1:2)= looktable(wdis,iwd,ix2, state)
          Ny%aa(1:2)= (/-1, 1/)*2.
       else if (rdis(i,j+1,k)<Rin .or. rdis(i,j+1,k)>Rout) then
            !left diff
          iw2=ijk_to_index(i,j,k)
          ix2=(/iw1,iw2/)
          Ny%ja(1:2)= looktable(wdis,iwd,ix2, state)
          Ny%aa(1:2)= (/-1, 1/)*2.
       else   !central diff
          ix2=(/iw1,iw2/)
          Ny%ja(:)= looktable(wdis,iwd,ix2, state)
          Ny%aa(:)= (/-1., 1./)
       endif
    ENDIF

  ! for dfdz=N(k+1)-N(k-1)
    IF (k==1) THEN    !;right diff
       Nz%n =2
       iw1=ijk_to_index(i,j,k)
       iw2=ijk_to_index(i,j,k+1)
       ix2=(/iw1,iw2/)
       Nz%ja(1:2)= looktable(wdis,iwd,ix2, state)
       Nz%aa(1:2)= (/-1, 1/)*2.
    ELSE IF (k==grid_sz) THEN
       Nz%n =2    !left diff
       iw1=ijk_to_index(i,j,k-1)
       iw2=ijk_to_index(i,j,k)
       ix2=(/iw1,iw2/)
       Nz%ja(1:2)= looktable(wdis,iwd,ix2, state)
       Nz%aa(1:2)= (/-1, 1/)*2.
    ELSE
       Nz%n =2  !central diff
       iw1=ijk_to_index(i,j,k-1)
       iw2=ijk_to_index(i,j,k+1)
       if (rdis(i,j,k-1)<Rin .or. rdis(i,j,k-1)>Rout) then
             !right diff
          iw1=ijk_to_index(i,j,k)
          ix2=(/iw1,iw2/)
          Nz%ja(1:2)= looktable(wdis,iwd,ix2, state)
          Nz%aa(1:2)= (/-1, 1/)*2.
       else if (rdis(i,j,k+1)<Rin .or. rdis(i,j,k+1)>Rout) then
            !left diff
          iw2=ijk_to_index(i,j,k)
          ix2=(/iw1,iw2/)
          Nz%ja(1:2)= looktable(wdis,iwd,ix2, state)
          Nz%aa(1:2)= (/-1, 1/)*2.
       else   !central diff
          ix2=(/iw1,iw2/)
          Nz%ja(:)= looktable(wdis,iwd,ix2, state)
          Nz%aa(:)= (/-1., 1./)
       endif
    ENDIF

    IF (nr(1,1) >0) THEN !add the radial weights
       Nx%aa=Nx%aa/ynr(im)
       Ny%aa=Ny%aa/ynr(im)
       Nz%aa=Nz%aa/ynr(im)
    ENDIF

    nja=Nx%n
    call Allocate_Rys(mRys(ii),nja)
    mRys(ii)%ja = Nx%ja(1:nja)
    mRys(ii)%aa = Nx%aa(1:nja)
    IA(ii+1) = IA(ii)+nja
    ii=ii+1

    nja=Ny%n
    call Allocate_Rys(mRys(ii),nja)
    mRys(ii)%ja = Ny%ja(1:nja)
    mRys(ii)%aa = Ny%aa(1:nja)
    IA(ii+1) = IA(ii)+nja
    ii=ii+1

    nja=Nz%n
    call Allocate_Rys(mRys(ii),nja)
    mRys(ii)%ja = Nz%ja(1:nja)
    mRys(ii)%aa = Nz%aa(1:nja)
    IA(ii+1) = IA(ii)+nja
    ii=ii+1
ENDDO row_loop

 mm=ii-1
 nnz = IA(mm+1)-1
 call Allocate_SP(R, mm, nn, nnz)
 R%ifst = IA
 DO i=1,mm
    R%icol(IA(i):IA(i+1)-1) = mRys(i)%ja
    R%A(IA(i):IA(i+1)-1) = mRys(i)%aa
    call DeAllocate_Rys(mRys(i))
 ENDDO
 call DeAllocate_Rys(Nx)
 call DeAllocate_Rys(Ny)
 call DeAllocate_Rys(Nz)
RETURN
END FUNCTION Reg1stsm2cdf_SP

FUNCTION Reg1stsm1rdf_SP (wdis, rdis, nr) RESULT(R)
!+ purpose: First-order smoothing using right-side differetiate, e.g.  dfdx(i)=N(i+1)-N(i)
! Input: n - dimension of grid
!       wdis -index for grid points with r>Rin
!       rdis - distance of grids to Sun center for grid [n,n,n]
! return: R -mm x nn sparse matrix for regularization, where
!          nn=ndis and mm=3*nn
! add weight=1/nr for each grid point remove r-dependence effect of
! regularization
USE InternPro
IMPLICIT NONE
INTEGER, INTENT(IN), DIMENSION(:) :: wdis
REAL, INTENT(IN), DIMENSION(:,:,:) :: rdis
REAL (kind=double), INTENT(IN), DIMENSION(:,:) :: nr
TYPE (SparseMat) :: R
TYPE (dataRys) :: Nx, Ny, Nz, mRys(3*size(wdis))
REAL (kind=double) :: Rin, Rout, ds
INTEGER :: nn, mm, nnz, i,j,k, iw1,iw2, ix2(2)
INTEGER, DIMENSION(size(wdis)) :: iwd, iw,jw,kw
REAL (kind=double), DIMENSION(size(wdis)) :: xnr, ynr
INTEGER, DIMENSION(3*size(wdis)+1) :: IA
INTEGER, ALLOCATABLE :: state(:)
INTEGER :: im, ii, nja

print*,'---Create regularized R_sprs using First-order right-side differetiate, |df/dx|=|(f(i+1)-f(i))/dx|'
ds=2.0_double*grid_Rout/dble(grid_sz-1)
Rin=grid_Rin - 2*ds
Rout=grid_Rout + 2.d0*ds
nn=size(wdis)
mm=3*nn
iwd=(/(i,i=1,nn)/)  !column number of R_sparse matrix

call Allocate_Rys(Nx, 2)
call Allocate_Rys(Ny, 2)
call Allocate_Rys(Nz, 2)

IF (nr(1,1)>0) THEN
   print*,'====> Radial weight Nr is applied to R_sparse!'
   call index_to_ijk(wdis,iw,jw,kw)
   DO i=1,nn
      xnr(i)=rdis(iw(i),jw(i),kw(i))
   ENDDO
   call interp(nr(:,1),nr(:,2),xnr,ynr)
ENDIF
  ii=1  ! row number of R_sparse matrix
  IA(1)=1  ! index in R%aa for 1st non-0 element in 1st row of R_sparse
row_loop: DO im=1, nn
   call index_to_ijk(wdis(im),i,j,k)
    !;for dfdx=N(i+1)-N(i)
    IF (i==grid_sz) THEN
       Nx%n =2  !left diff
       iw1=ijk_to_index(i-1,j,k)
       iw2=ijk_to_index(i,j,k)
       ix2=(/iw1,iw2/)
       Nx%ja(1:2)= looktable(wdis,iwd,ix2, state)
       Nx%aa(1:2)= (/-1, 1/)
    ELSE
       Nx%n =2
       iw1=ijk_to_index(i,j,k)
       iw2=ijk_to_index(i+1,j,k)
       if (rdis(i+1,j,k)<Rin .or. rdis(i+1,j,k)>Rout) then
          iw1=ijk_to_index(i-1,j,k)   !left diff.
          iw2=ijk_to_index(i,j,k)
       endif
       ix2=(/iw1,iw2/)
       Nx%ja(1:2)= looktable(wdis,iwd,ix2, state)
       Nx%aa(1:2)= (/-1, 1/)
    ENDIF

 !for dfdy=N(j+1)-N(j)
    IF (j==grid_sz) THEN
       Ny%n =2  !left diff
       iw1=ijk_to_index(i,j-1,k)
       iw2=ijk_to_index(i,j,k)
       ix2=(/iw1,iw2/)
       Ny%ja(1:2)= looktable(wdis,iwd,ix2, state)
       Ny%aa(1:2)= (/-1, 1/)
    ELSE
       Ny%n =2
       iw1=ijk_to_index(i,j,k)
       iw2=ijk_to_index(i,j+1,k)
       if (rdis(i,j+1,k)<Rin .or. rdis(i,j+1,k)>Rout) then
          iw1=ijk_to_index(i,j-1,k)   !left diff.
          iw2=ijk_to_index(i,j,k)
       endif
       ix2=(/iw1,iw2/)
       Ny%ja(1:2)= looktable(wdis,iwd,ix2, state)
       Ny%aa(1:2)= (/-1, 1/)
    ENDIF

   !for dfdz=N(k+1)-N(k)
    IF (k==grid_sz) THEN
       Nz%n =2  !left diff
       iw1=ijk_to_index(i,j,k-1)
       iw2=ijk_to_index(i,j,k)
       ix2=(/iw1,iw2/)
       Nz%ja(1:2)= looktable(wdis,iwd,ix2, state)
       Nz%aa(1:2)= (/-1, 1/)
    ELSE
       Nz%n =2
       iw1=ijk_to_index(i,j,k)
       iw2=ijk_to_index(i,j,k+1)
       if (rdis(i,j,k+1)<Rin .or. rdis(i,j,k+1)>Rout) then
          iw1=ijk_to_index(i,j,k-1)   !left diff.
          iw2=ijk_to_index(i,j,k)
       endif
       ix2=(/iw1,iw2/)
       Nz%ja(1:2)= looktable(wdis,iwd,ix2, state)
       Nz%aa(1:2)= (/-1, 1/)
    ENDIF

    IF (nr(1,1) >0) THEN !add the radial weights
       Nx%aa=Nx%aa/ynr(im)
       Ny%aa=Ny%aa/ynr(im)
       Nz%aa=Nz%aa/ynr(im)
    ENDIF

    nja=Nx%n
    call Allocate_Rys(mRys(ii),nja)
    mRys(ii)%ja = Nx%ja(1:nja)
    mRys(ii)%aa = Nx%aa(1:nja)
    IA(ii+1) = IA(ii)+nja
    ii=ii+1

    nja=Ny%n
    call Allocate_Rys(mRys(ii),nja)
    mRys(ii)%ja = Ny%ja(1:nja)
    mRys(ii)%aa = Ny%aa(1:nja)
    IA(ii+1) = IA(ii)+nja
    ii=ii+1

    nja=Nz%n
    call Allocate_Rys(mRys(ii),nja)
    mRys(ii)%ja = Nz%ja(1:nja)
    mRys(ii)%aa = Nz%aa(1:nja)
    IA(ii+1) = IA(ii)+nja
    ii=ii+1
! print*,'----OK for im=',im
ENDDO row_loop

 mm=ii-1
 nnz = IA(mm+1)-1
 call Allocate_SP(R, mm, nn, nnz)
 R%ifst = IA
 DO i=1,mm
    R%icol(IA(i):IA(i+1)-1) = mRys(i)%ja
    R%A(IA(i):IA(i+1)-1) = mRys(i)%aa
    call DeAllocate_Rys(mRys(i))
 ENDDO
 call DeAllocate_Rys(Nx)
 call DeAllocate_Rys(Ny)
 call DeAllocate_Rys(Nz)
RETURN
END FUNCTION Reg1stsm1rdf_SP

END MODULE TomoPro



!Main ------------------------------------------------
program tomo_cart
 use InternPro
 use TomoPro
 implicit none

! Final version for Cartesian-grid code, modified on 03/18/2025 by Tongjiang Wang
!  the solution of electron density in units of cm^-3
integer, parameter :: nimg=60  
integer :: nt !number of used images
real (kind=8), parameter :: Lse=215.0_double 
real (kind=8), parameter ::  epswt= 1.e-4   !threshold for none-zero elements in matrix A
real (kind=double) :: ds, xcar(grid_sz)
real, dimension(grid_sz,grid_sz,grid_sz) :: xp,yp,zp, rdis, weights, nsolu
integer :: i,j,k, isz, ii, jj, mm, nn, nnz
integer, dimension(:), allocatable:: px1d,tdis,wdis,iw,jw,kw
integer :: ndis,ind, nt0, im,ik, npw0, npw, nxr, ywt=1
real (kind=double),  dimension(:,:), allocatable :: ix,iy, nr
real (kind=double), dimension(:,:), allocatable :: xcrd,ycrd,irp, cobs_xy
real (kind=double) :: Lobs(nimg), cobs(2,nimg), cobs0(2)
type (dataRys) :: arrRys(nimg*grid_sz**2)
integer, dimension(:), allocatable :: IAt, IA
integer :: allocatestatus, Openstatus
type (mapstr) :: map
type (SparseMat) :: Amat, Rmat
real (kind=double), dimension(:), allocatable :: ixw,iyw,izw,irw,irho,ipb,&
                        cx,cy,cz,lon,lat,dse, yd0,ydw, xr,ybk
real, dimension(:), allocatable :: ydat0,ydata
real (kind=double), dimension(:), allocatable:: xe,ye,ze
integer, dimension(:,:), allocatable:: ccell
real, dimension(2) :: xt,yt,zt
integer :: cross, nja, istt,iend 
real (kind=double), dimension(3) :: ni
real, dimension(grid_sz*grid_sz) :: weight0,px0,py0,pz0
integer, dimension(grid_sz*grid_sz) :: wja
real, dimension(:), allocatable:: weightnz,px1,py1,pz1,rp2,sp,Rho,tht,kr
logical :: is_spb=.TRUE.
integer :: write_Abuff, read_Abuff, write_Rbuff, read_Rbuff, wt_nr
real :: mu,tol
real, dimension(:), allocatable :: xsolu
integer :: sav_res,sav_iter,resume
type (dataRys) :: row
integer ::  i0,j0, open_status
real, dimension(:), allocatable :: ydata1
INTEGER :: ia0, ia1, npa, i1,i2
real :: yd00
real :: tstart,tfinish, e,t(2)

real (kind=8) :: amax,amin,amin0=1e10,amax0=1e-10, n_unit=1.0e7
integer, allocatable, dimension(:) :: wmin,wmax
integer :: nmax,nmin

character (len=:),allocatable :: Abuffname, Rbuffname
character (len=:),allocatable :: sdir
character (len=200) :: filename,ybkname,nrbkname,obsname, dir
character (len=200) :: namesolu, pbnames(nimg)
integer :: imgnum
logical :: exist


call cpu_time(tstart)
e=dtime(t)
print*,'elapsed:', e, ', user:', t(1), ', sys:', t(2)

! read parameters for run
  read(*,*) mu
  read(*,*) tol
  read(*,*) resume
  read(*,*) sav_res
  read(*,*) sav_iter
  read(*,*) wt_nr
  read(*,*) write_Abuff
  read(*,*) read_Abuff
  read(*,*) write_Rbuff
  read(*,*) read_Rbuff
  read(*,*) dir
  read(*,*) obsname
  read(*,*) ybkname
  read(*,*) nrbkname
  read(*,*) imgnum
  nt=imgnum
  print*,'the number of used pB images=',nt

  read(*,*) (pbnames(i),i=1,imgnum)

 sdir=trim(adjustL(dir))
 obsname=sdir//trim(adjustL(obsname))
 ybkname=sdir//trim(adjustL(ybkname))
 nrbkname=sdir//trim(adjustL(nrbkname))
 do i=1, imgnum
    pbnames(i)=sdir//trim(adjustL(pbnames(i)))
!   print*,pbnames(i)
 enddo

! set up parameters for running
!wt_nr=1
!write_Abuff=1
!read_Abuff =0
!write_Rbuff=1
!read_Rbuff =0

 Abuffname='Buff_Asparse_sz128.dat'
 if (wt_nr == 1) then
    Rbuffname='Buff_Rsp_sz128_Reg2_nwt.dat'
    namesolu='nsolu_sz128_reg2_nwt.dat'
 else
    Rbuffname='Buff_Rsp_sz128_Reg2.dat'
    namesolu='nsolu_sz128_reg2.dat'
 endif

!  mu=0.3  !best mu for wt_nr=1 from chi_abs CV
! tol=0.01_double !for mu=0.3
! resume=0
! sav_res=1
! sav_iter=1


 ds=2.d0*grid_Rout/dble(grid_sz-1)
 xcar=(/(-grid_Rout+ds*(i-1),i=1,grid_sz)/)
 do 100 j=1,grid_sz
 do 100 i=1,grid_sz
    xp(:,i,j)=xcar
    yp(i,:,j)=xcar
    zp(i,j,:)=xcar
100 continue
  rdis=sqrt(xp**2+yp**2+zp**2)

  allocate(tdis(grid_sz**3), STAT=allocatestatus)
  IF (allocateStatus /= 0) STOP "*** Not enough memory ***"

  ind=1
  do 200 k=1,grid_sz
  do 200 j=1,grid_sz
  do 200 i=1,grid_sz
       if (rdis(i,j,k) >= grid_Rin-2.0*ds .and. rdis(i,j,k) <= grid_Rout+2.0*ds) then
            tdis(ind)= ijk_to_index(i,j,k)
            ind=ind+1
       endif
200 continue
 ndis=ind-1
 ALLOCATE(wdis(ndis))
 wdis=tdis(1:ndis)
 DEALLOCATE(tdis)

 ALLOCATE(iw(ndis),jw(ndis),kw(ndis))
 call index_to_ijk(wdis, iw,jw,kw) 

 IF (read_Abuff ==1) GO TO 400
! [iw,jw,kw] 3D indice for grid-points in the reconstructed domain

!define image coord. for observer
print*,'... Reading pB-images into map-struct----' 
  isz=grid_sz     !sample resolution not more than grid resolution, keeping solu.  underdetermined
  nn=ndis         !set initial value for row number of A_sparse
  mm=nt*isz*isz   !set initial value for column number of A_sparse
  ALLOCATE(px1d(isz),ix(isz,isz),iy(isz,isz),xcrd(isz,isz),ycrd(isz,isz),irp(isz,isz),&
          ixw(isz*isz),iyw(isz*isz),izw(isz*isz),irw(isz*isz),yd0(isz*isz),&
          ydat0(nt*isz*isz), IAt(nt*grid_sz**2))
  px1d=(/(i,i=1,isz)/) !index for image pixel coordinate
! create pix coord. for image
 DO k=1, isz
    ix(:,k)=real(px1d)
    iy(k,:)=real(px1d)
 ENDDO
! Reading (distance, lon, lat) for Observer, 
!    where units are Rsun and degree
 OPEN (UNIT=10, FILE=obsname, ACTION='READ',&
            IOSTAT = Openstatus )
 IF (OpenStatus /= 0) THEN
      STOP "error opening file [cobs_coord_?.txt]"
 ELSE
    print*, '---file ['//obsname//'] is read ....'
 ENDIF

 READ(10,*) nt0
! print*,'nt0=',nt0
 DO k=1,nt0
   READ(10,*) Lobs(k), cobs(:,k)
!   print*,Lobs(k), cobs(:,k)
 ENDDO
 CLOSE(10)

 ! read the radially-dependent background pB obtain from obs
  OPEN (UNIT=20, FILE=ybkname, ACTION='READ', &
        IOSTAT = Openstatus )
  IF (OpenStatus /= 0) THEN
          STOP "error opening file [pbbk_xr_?.txt]"
  ELSE
     print*, '---file ['//ybkname//'] is read ....'
  ENDIF

  READ(20,*) nxr
  ALLOCATE(xr(nxr),ybk(nxr))
  READ(20,*) xr, ybk
  CLOSE(20)

 
 ii=1
 IAt(ii)=1
 npw0=1
 img_loop: DO im=1, nt  !main loop
    filename=pbnames(im)
    print*,'--read img no=',im,'  file=',filename
    map=read_mapstr(filename)
    call getmap_coord(map,ix,iy, xcrd,ycrd)
    xcrd=xcrd/map%rsun !convert img coord from arcsec into Rsun
    ycrd=ycrd/map%rsun
    irp=sqrt(xcrd*xcrd+ycrd*ycrd)
    jj=1
    DO 300 j=1,isz
    DO 300 i=1,isz
       IF (irp(i,j)>=grid_Rin .and. irp(i,j)<=grid_Rout) THEN
           ixw(jj)=xcrd(i,j)
           iyw(jj)=ycrd(i,j)
           irw(jj)=irp(i,j)
           yd0(jj)=map%data(i,j)
           jj=jj+1
       ENDIF
    300 CONTINUE
    npw=jj-1
    izw=Lobs(im)
!Lobs is the distance from Sun center to Observer in unit of Rsun
!   ie., Lobs=Lcor1a or Lcor1b
    
    ALLOCATE(cx(npw),cy(npw),cz(npw),irho(npw),ipb(npw),ydw(npw),&
             lon(npw),lat(npw),dse(npw),cobs_xy(2,npw))
    irho=irw(1:npw)  !distanc from Sun-center to pixel in img-plane
    ydw =yd0(1:npw)  !pB-data for a pixel within given domain
!   Determine the ray direction of each pixel in Carrington coordinates    
    cobs0= cobs(:,im)
    call img2car_coord(cobs0,-ixw(1:npw),-iyw(1:npw),izw(1:npw), cx,cy,cz)     
    call car2sp(cx,cy,cz, lon,lat,dse)  ![lon,lat] in radian
!     dse id the distance from observer to a pixel in img-plane
    cobs_xy(1,:)=mod(lon*180.d0/PI+360.d0,360.d0)
    cobs_xy(2,:)=lat*180.d0/PI
    
!  if it is set to normalize the obsrved and modeled data 
!   by the radial pB background obtained from obs (ywt=1 as default)
!   IF (ywt ==1) THEN
       call interp(xr,ybk, irho,ipb)
       ydw=ydw/ipb  !normalize the data by radial pB background
!   ENDIF
   ydat0(npw0:npw0+npw-1)=ydw
   npw0= npw0+ npw
      pix_loop: DO ik=1,npw  ! second loop
          cobs0=cobs_xy(:,ik)
          call get_cross_points (cobs0, ixw(ik),iyw(ik),0.d0,rdis,is_spb,&
               xe,ye,ze,ccell, cross, ni,xt,yt,zt)     
          cross_if: IF (cross ==1) THEN
! Contribution of the grid to line integration for a ray
             weights=Rsun*ds*ray_weights(xe,ye,ze,ccell) !in units of cm.

! Add the contribution of Thomson scatter function for single electron
             ! find non-zero pixeks in 3d-weights
             k=1
             DO i=1,ndis
                if (weights(iw(i),jw(i),kw(i))>epswt) then
                   weight0(k)=weights(iw(i),jw(i),kw(i))
                   px0(k) = xp(iw(i),jw(i),kw(i)) 
                   py0(k) = yp(iw(i),jw(i),kw(i))   
                   pz0(k) = zp(iw(i),jw(i),kw(i))                                                   
                   wja(k) = i  !indice of non-zero elements in wdis
                   k=k+1   
                endif     
             ENDDO
             nja=k-1  ! number of non-zero elements in weights
!  ;distance from Sun center to point P, rp=|OP(px1,py1,pz1)| =sqrt(rp2)
             ALLOCATE(weightnz(nja),px1(nja),py1(nja),pz1(nja),&
                    rp2(nja),sp(nja),Rho(nja),tht(nja),kr(nja))
             weightnz =weight0(1:nja)
             px1=px0(1:nja)
             py1=py0(1:nja)
             pz1=pz0(1:nja)
             rp2=px1*px1+py1*py1+pz1*pz1

!       Calculate the radial-dependent weight function wi=1/Yi(r)_bk
!  note: apply w=1/ypb here to weightnz to avoid values of "weightnz"
!  too small leading overfloating a range of 10^+/-38
!            IF (ywt ==1) weightnz =weightnz/ipb(ik)
             weightnz =weightnz/ipb(ik)

! sp is the distance from P to POS
!       where POS (plane-of-sky) is a plane that crosses the Sun center and
!        perpendicular to LOS, s=|vec_OP cdot vec_ni|=|rp*cos(alp)|, where alp is
!        the angle between OP and LOS, and vec_ni is the unit vector of LOS.
             sp=abs(px1*ni(1)+py1*ni(2)+pz1*ni(3)) 
             Rho=sqrt(rp2-sp*sp)  ! rho is the distance from the Sun center to LOS
             tht=atan(sp/Rho)     ! tht is the angle between vector OP and POS
             kr =thomsonscatter(Rho,tht,-1.,.true.) !tht in unit of radian
             weightnz = weightnz*(kr*n_unit)    !add effect of Thmoson scatter, I=Sum(Kr*Nr*dl)
                                            !set density in unit of 10^7 cm^-3

         !construct Matrix A in CRS format (CRS: Compressed Row Sparse)
             call Allocate_Rys(arrRys(ii), nja)    
             arrRys(ii)%ja =wja(1:nja)
             arrRys(ii)%aa =weightnz
             IAt(ii+1) = IAt(ii)+nja
             ii=ii+1
          ELSE
             STOP "---error with sub: get_cross_points"
          ENDIF cross_if
         DEALLOCATE(weightnz,px1,py1,pz1,rp2,sp,Rho,tht,kr)
      ENDDO  pix_loop
      DEALLOCATE(cx,cy,cz,irho,ipb,ydw,lon,lat,dse,cobs_xy)
 ENDDO  img_loop
  mm=ii-1  ! total number of used rays or data points
! nn=ndis  ! total number of wdis or grid points in Domain
  ALLOCATE(ydata(mm), IA(mm+1))
  ydata=ydat0(1:mm)
  IA=IAt(1:mm+1)  
  nnz =IA(mm+1)-1   !total number of non-zero elements in A_sparse matrix
  DEALLOCATE(px1d,ix,iy,xcrd,ycrd,irp,ixw,iyw,izw,irw,yd0,ydat0,IAt)

e=dtime(t)
print*,'elapsed:', e, ' (s)  for creating mat A_sparse----'

! save Data for Coeff. Matrix A in a buff file on hard driver
 IF (write_Abuff ==1) THEN
    OPEN(1, FILE=Abuffname, FORM="unformatted",ACTION="write", STATUS="new") 
    WRITE(1) mm,nn,nnz
    WRITE(1) IA
    WRITE(1) (arrRys(i)%ja, i=1,mm)
    WRITE(1) (arrRys(i)%aa, i=1,mm)
    WRITE(1) ydata
    CLOSE(1)
   DeAllocate(IA,ydata)
   DO i=1,mm
     call DeAllocate_Rys (arrRys(i))
   ENDDO
   print*,'Buff file: ['//Abuffname//'] is written! ===>'
 ENDIF


400 IF (write_Abuff==1 .or. read_Abuff==1) THEN
   print*,' ----Reading Buff for matA=====>'
   OPEN(10, FILE=Abuffname,FORM='unformatted',ACTION='read')
   READ(10) mm,nn,nnz
   call Allocate_SP (Amat,mm,nn,nnz)
   Allocate(ydata(mm))
   READ(10) Amat%ifst
   READ(10) Amat%icol
   READ(10) Amat%A
   READ(10) ydata
   CLOSE(10)
   print*,'Buff file: ['//Abuffname//'] is read!==='
ELSE !if both write_Abuff==0 .and. read_Abuff==0
   call Allocate_SP (Amat,mm,nn,nnz)
   Amat%ifst =IA
   DO i=1,mm
      istt=IA(i)
      iend=IA(i+1)-1
      Amat%icol(istt:iend)=arrRys(i)%ja
      Amat%A(istt:iend)   =arrRys(i)%aa
      call DeAllocate_Rys (arrRys(i))
   ENDDO
  DeAllocate(IA)
ENDIF 

 print*,'mm,nn,nnz=', Amat%m, Amat%n, Amat%asiz
 print*,'min max A=',minval(Amat%A),maxval(Amat%A)
 print*,'Amean=',sum(Amat%A)/Amat%asiz
 print*,'min max ydata=',minval(ydata),maxval(ydata) 

 IF (wt_nr ==1) THEN 
     OPEN(20, FILE=nrbkname ,ACTION='read') 
     READ(20,*) nxr
     ALLOCATE(nr(nxr,2))
     READ(20,*) nr(:,1), nr(:,2)
     CLOSE(20)
     print*,'... reading file ['//nrbkname//'] for radial weighting ....'
 ELSE
    print*,'No radial weighting will be applied in Regularization!'
    ALLOCATE(nr(1,2))    
    nr=0.d0
 ENDIF

 IF (read_Rbuff==1) THEN
   print*,'---Reading buff for matR ===>'
   OPEN(10, FILE=Rbuffname,FORM='unformatted',ACTION='read')
   READ(10) mm,nn,nnz
   call Allocate_SP (Rmat,mm,nn,nnz)
   READ(10) Rmat%ifst
   READ(10) Rmat%icol
   READ(10) Rmat%A
   CLOSE(10)
   print*,'Buff file: ['//Rbuffname//'] is read!==='
 ELSE
!  Rmat=Reg1stsm1rdf_SP (wdis, rdis, nr)
!  Rmat=Reg2ndsm4cdf_SP (wdis, rdis, nr)
   Rmat=Reg2ndsm2cdf_SP (wdis, rdis, nr)
   print*,' mat R_sparse is created: mm,nn,nnz=',Rmat%m, Rmat%n, Rmat%asiz
   print*,'min max R=',minval(Rmat%A),maxval(Rmat%A)
! save Data for regularization R_sparse in a buff file 
   IF (write_Rbuff ==1) THEN
      OPEN(1, FILE=Rbuffname, FORM="unformatted",ACTION="write",STATUS="new")
      WRITE(1) Rmat%m, Rmat%n, Rmat%asiz
      WRITE(1) Rmat%ifst
      WRITE(1) Rmat%icol
      WRITE(1) Rmat%A
      CLOSE(1)
     print*,'Buff file: ['//Rbuffname//'] written! ===>'
   ENDIF
 ENDIF
 e=dtime(t)
 print*,'elapsed:', e, ' (s)  for creating mat R_sparse----'

 print *,'----- Find the solution using Conjugate Gradient (CG)'
 print *,' f=min(||Ax-y||^2+mu*||R*x||^2 '

 ALLOCATE(xsolu(Amat%n))
 print*,'.... with mu=',mu, '  tol=',tol

xsolu= ConjugateGrad_SP (Amat,ydata,mu,Rmat, tol,sav_res,sav_iter,resume)


  !Note: indice (iw,jw,kw) for 3D grid in range [1, grid_sz]
if (size(xsolu)/=ndis) STOP "There is an error for wdis and xsolu"
 DO i=1,ndis
    nsolu(iw(i),jw(i),kw(i))=xsolu(i) * n_unit  !density in cm^-3
 ENDDO
print*,'solu min max=',minval(xsolu),maxval(xsolu)

e=dtime(t)
 print*,'elapsed:', e, ' (s)  for obtaining solution using conjugate gradient----'

 INQUIRE(FILE=namesolu, exist=exist)
     if (.not. exist) then
         open(unit=10,file=namesolu, form='unformatted',status="new",action="write",iostat=open_status)
     else
         print*,'The file ['//namesolu//'] already exists!'     
         namesolu='New_'//namesolu
         print*,'The solution will be saved in ['//namesolu//']--------!'
         open(unit=10,file=namesolu, form='unformatted',status="new",action="write",iostat=open_status) 
     endif

if (open_status >0) stop "----Error, File not opened properly. &
                         &it maybe already exist---"
write(10) grid_sz
write(10) xcar
write(10) nsolu  
close(10)
print*, 'Solution file ['//namesolu//'] is created, including electron density in cm^-3!'

call cpu_time(tfinish)
print*,'computing time=',tfinish-tstart

end program tomo_cart 
