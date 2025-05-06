!---------------------------------------------------------------------------------------------
! <Ne3dTomo> provides the tomographic reconstruction of 3D electron density of the corona
! with an option by applying radial weighting factor to the regularization term to provide
! balanced smoothing at different heights and help stablizing the solution. The detailed
! method and applications are referred to the publication "Wang, T., Arge, C.N, & Jones, I.J.
!  Sol.Phys., (2025) 300:46. DOI:10.1007/s11207-025-02454-8 
! Copyright (C) 2025 by Tongjiang Wang (wangt@cua.edu)
!
! Second version released on May 5, 2025: Spheric grid with OpenMP 
! (download available from https://stereo-ssc.nascom.nasa.gov/data/ins_data/secchi/N3D_COR1A)
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

   INTERFACE car2sp
      MODULE PROCEDURE car2sp_1d, car2sp_scl
   END INTERFACE

   INTERFACE sp2car
      MODULE PROCEDURE sp2car_1d, sp2car_scl
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
End Function 


SUBROUTINE locate(xx,x,j)
   IMPLICIT NONE
   REAL (kind=double), INTENT(IN) :: xx(:),x
   INTEGER, INTENT(OUT) :: j
   INTEGER :: n, jl, ju, jm
! Given an array xx(1:n), and a given value x, return a value j such
!that x is between xx(j) and xx(j+1). xx(1:n) must be monotonic, either
!increasing or decreasing. j=0 or j=n is return to indicate that x is out of
!range
! this routine called by interp()
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
use OMP_LIB
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
!$OMP PARALLEL WORKSHARE 
 x=t11*xp + t12*yp +t13*zp
 y=t21*xp + t22*yp +t23*zp
 z=t31*xp + t32*yp +t33*zp
!$OMP END PARALLEL WORKSHARE 
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

SUBROUTINE car2sp_1d (x,y,z,lon,lat,r)
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
USE OMP_LIB
IMPLICIT NONE
REAL (kind=double), INTENT(IN), DIMENSION(:) :: x,y,z
REAL (kind=double), INTENT(OUT), DIMENSION(size(x)) :: lon,lat,r
!$OMP PARALLEL WORKSHARE
 r=sqrt(x*x+y*y+z*z)
 lat=asin(z/r)
 lon=atan2(y,x)
!$OMP END PARALLEL WORKSHARE
END SUBROUTINE car2sp_1d

SUBROUTINE car2sp_scl (x,y,z,lon,lat,r)
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
REAL (kind=double), INTENT(IN) :: x,y,z
REAL (kind=double), INTENT(OUT) :: lon,lat,r

 r=sqrt(x*x+y*y+z*z)
 lat=asin(z/r)
 lon=atan2(y,x)
END SUBROUTINE car2sp_scl


SUBROUTINE  sp2car_1d (lon,lat,r,x,y,z)
!;+ Purpose
!; transform Carrington coord. from Spheric sys. to Cartesian sys.
!;Input:
!; lon - longitude in rad, range in [-pi, pi] or [0, 2*pi], measured from x-axis;
!;       positive in anti-clockwise
!; lat - latitude in rad, range in [-pi/2, pi/2]
!; r   - heliocentric distance in unit same as x,y,or z
!;Output:  x,y,z - cartesian coord. z defined as solar axis toward N-pole
!;      [x,y] in equatial plane
 use OMP_LIB
IMPLICIT NONE
REAL (kind=double), INTENT(IN), DIMENSION(:) :: lon,lat,r
REAL (kind=double), INTENT(OUT), DIMENSION(size(r)) :: x,y,z
!$OMP PARALLEL WORKSHARE
 x=r*cos(lat)*cos(lon)
 y=r*cos(lat)*sin(lon)
 z=r*sin(lat)
!$OMP END PARALLEL WORKSHARE 
END SUBROUTINE sp2car_1d

SUBROUTINE  sp2car_scl (lon,lat,r,x,y,z)
!;+ Purpose
!; transform Carrington coord. from Spheric sys. to Cartesian sys.
!;Input:
!; lon - longitude in rad, range in [-pi, pi] or [0, 2*pi], measured from x-axis;
!;       positive in anti-clockwise
!; lat - latitude in rad, range in [-pi/2, pi/2]
!; r   - heliocentric distance in unit same as x,y,or z
!;Output:  x,y,z - cartesian coord. z defined as solar axis toward N-pole
!;      [x,y] in equatial plane
IMPLICIT NONE
REAL (kind=double), INTENT(IN)  :: lon,lat,r
REAL (kind=double), INTENT(OUT) :: x,y,z

 x=r*cos(lat)*cos(lon)
 y=r*cos(lat)*sin(lon)
 z=r*sin(lat)
END SUBROUTINE sp2car_scl


FUNCTION findmax_D (a,w,n) result(amax)
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
  return
END FUNCTION findmax_D

FUNCTION findmax_R (a,w,n) result(amax)
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
END FUNCTION findmax_R


FUNCTION findmax_I (a,w,n) result(amax)
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
END FUNCTION findmax_I


FUNCTION findmin_D (a,w,n) result(amin)
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
END FUNCTION findmin_D

FUNCTION findmin_R (a,w,n) result(amin)
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
END FUNCTION findmin_R


FUNCTION findmin_I (a,w,n) result(amin)
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
END FUNCTION findmin_I

END MODULE InternPro


MODULE TomoPro
  USE InternPro, only: double
  IMPLICIT NONE

!INTEGER, PARAMETER :: double=SELECTED_REAL_KIND(15)
REAL (kind=double), PARAMETER :: Rsun=6.96d10 !cm  solar radius
REAL (kind=double), PARAMETER :: grid_Rout=4.0_double
REAL (kind=double), PARAMETER :: grid_Rin=1.5_double
INTEGER, PARAMETER :: sz_r=51, sz_th=181, sz_ph=361
!INTEGER, PARAMETER :: sz_r=51, sz_th=91, sz_ph=181
!INTEGER, PARAMETER :: sz_r=6, sz_th=26, sz_ph=51 

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

  REAL(kind=double), INTENT(IN), DIMENSION(:) :: Rin, Tht
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
  REAL(kind=double), INTENT(IN) :: Rin, Tht
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
  USE OMP_LIB        
  IMPLICIT NONE
  TYPE (Mapstr), INTENT(IN) :: map
  REAL (kind=double), INTENT(IN), DIMENSION(:) :: xp,yp
  REAL (kind=double), INTENT(OUT),DIMENSION(:) :: xcoord,ycoord
  REAL :: xdim, ydim
!$OMP PARALLEL WORKSHARE  
  xdim=size(map%data,1)
  ydim=size(map%data,2)
  xcoord=(xp-(xdim+1.d0)/2.d0)*map%dx +map%xc
  ycoord=(yp-(ydim+1.d0)/2.d0)*map%dy +map%yc
!$OMP END PARALLEL WORKSHARE
  RETURN
END SUBROUTINE getmap_coord_1d

SUBROUTINE getmap_coord_2d(map,xp,yp,xcoord,ycoord)
  USE OMP_LIB
  IMPLICIT NONE
  TYPE (Mapstr), INTENT(IN) :: map
  REAL(kind=double), INTENT(IN), DIMENSION(:,:) :: xp,yp
  REAL (kind=double), INTENT(OUT),DIMENSION(:,:) :: xcoord,ycoord
  REAL :: xdim, ydim
!$OMP PARALLEL WORKSHARE  
  xdim=size(map%data,1)
  ydim=size(map%data,2)
  xcoord=(xp-(xdim+1.d0)/2.d0)*map%dx +map%xc
  ycoord=(yp-(ydim+1.d0)/2.d0)*map%dy +map%yc
!$OMP END PARALLEL WORKSHARE  
  RETURN
END SUBROUTINE getmap_coord_2d

SUBROUTINE index_to_ijk_1d(index1d,iw,jw,kw)
 ! convert 1D-index of array(I,J,K) in range=[1:I*J*K] into 3 indice [i,j,k]
! note: where indice variables,index1d,iw,jw,kw, start from 1
! note: index() is an intrinsic func,so replace variable 'index' as 'index1d' 
 IMPLICIT NONE
 INTEGER, INTENT(IN), DIMENSION(:) :: index1d
 INTEGER, INTENT(OUT), DIMENSION(size(index1d)) :: iw,jw,kw 
 INTEGER :: tmp, i
   DO i=1, size(index1d)
      kw(i)=(index1d(i)-1)/(sz_th*sz_ph)+1
      tmp=mod((index1d(i)-1), sz_th*sz_ph)
      jw(i)=tmp/sz_ph+1
      iw(i)=mod(tmp, sz_ph)+1
   ENDDO
END SUBROUTINE index_to_ijk_1d

SUBROUTINE index_to_ijk_scl(index1d,iw,jw,kw)
! convert one index of array(I,J,K) with value in range=[1:I*J*K]
!   into 3 indice [i,j,k] where i=[1,I], j=[1,J], k=[1,K]
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: index1d
 INTEGER, INTENT(OUT) :: iw,jw,kw
 INTEGER :: tmp
      kw=(index1d-1)/(sz_th*sz_ph)+1
      tmp=mod((index1d-1), sz_th*sz_ph)
      jw=tmp/sz_ph+1
      iw=mod(tmp, sz_ph)+1
END SUBROUTINE index_to_ijk_scl

FUNCTION ijk_to_index (i,j,k) RESULT(index1d)
!convert indice of 3D array data[i,j,k] where i=[1,I], j=[1,J], k=[1,K]
! into 1D index data1d(m) with m=0, I*J*K-1
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: i,j,k
  INTEGER :: index1d
  index1d = (i-1) + (j-1)*sz_ph + (k-1)*sz_ph*sz_th
RETURN
END FUNCTION ijk_to_index


SUBROUTINE sph_cross_points (cobs,xi,yi,zi, phe,the,re,ccell, &
              cross, ni, pht,tht,rt)
!;+ for a given point with image coord. [xi,yi,zi] and an observer with Carrington coord
!; cobs[2]=[lon, lat], determine the entry point of LOS (defined from the object to
!; observer) at the boundary of grid box
!; input: cobs[2] -observer's position in degrees
!;        xi,yi,zi - image coord. of a point
!; output:
!;  phe[N],the[N],re[N] - Carrington spherical coords. of N cross-points
!;  ccell[3,N-1] - index for N-1 crossing cells, defined as
!;      the vortex with the smaller ph,th,r
!;  Note: index of iph,ith,or ir start with 0
!; keyword: cross =1, return entry point; cross=0, not crossing.
!;          if cross =1, calculate Carrington coord of entry and exit points
!;      pht[2],tht[2],rt[2] at r=Rout
!;      ni[3] - direction vector of the LOS
!; most updated ver.: add treatment for Af2(i)==0 in Step2
 USE InternPro
 USE m_mrgrnk  !merge sort module
 USE OMP_LIB
 IMPLICIT NONE

 REAL (kind=double), INTENT(IN), DIMENSION(2) :: cobs
 REAL (kind=double), INTENT(IN) :: xi,yi,zi
 REAL (kind=double), INTENT(OUT), DIMENSION(:), allocatable :: phe,the,re
 INTEGER, INTENT(OUT), DIMENSION(:,:), allocatable :: ccell
 INTEGER, INTENT(OUT) :: cross
 REAL (kind=double), INTENT(OUT), DIMENSION(3) :: ni
 REAL (kind=double), INTENT(OUT), DIMENSION(2) :: rt,tht,pht 
 
 INTEGER :: i
 REAL (kind=double) :: dr,dth,dph,x0,y0,z0, ri(sz_r),thi(sz_th),phi(sz_ph)
 REAL (kind=double) :: phc,thc,a,b,c
 REAL (kind=double) :: Af1,Bf1,Cf1(sz_r),dlt1(sz_r),dltw1(sz_r),tn1(sz_r),tp1(sz_r),t1(2*sz_r)
 INTEGER :: nw, w1(sz_r), n_ri, n_thi, n_phi
 REAL (kind=double) :: r_ri(2*sz_r), th_ri(2*sz_r), ph_ri(2*sz_r)
 REAL (kind=double) :: Af2(sz_th),Bf2(sz_th),Cf2(sz_th),ctg2(sz_th),dlt2(sz_th),dltw2(sz_th),&
                 tn2(sz_th),tp2(sz_th),t2(2*sz_th)
 INTEGER :: w21(sz_th),w2(sz_th), thi_cross
 REAL (kind=double) :: r_thi(2*sz_th), th_thi(2*sz_th), ph_thi(2*sz_th)
 INTEGER :: w31(sz_ph), w32(sz_ph),phi_cross, nw1,nw2, pnum, AllocateStatus
 REAL (kind=double) :: ph_min, ph_max, r_phi(sz_ph), th_phi(sz_ph), ph_phi(sz_ph), ti 
 REAL (kind=double), DIMENSION(:), ALLOCATABLE :: r_all, th_all, ph_all, xe,ye,ze
 INTEGER,DIMENSION(:),ALLOCATABLE :: order
 REAL (kind=double) :: xm,ym,zm,phm,thm,rm, slen, xt(2),yt(2),zt(2),tlen
 INTEGER :: iph,ith,ir, thi_singular
 REAL :: t0, r0_thi
 REAL, PARAMETER:: epsn=1e-10

 cross=1

  dr=(grid_Rout-grid_Rin)/dble(sz_r-1)  !grid scale in r-direction
 dth=PI/dble(sz_th-1)                   !grid scale in th-direction
 dph=2*PI/dble(sz_ph-1)                 !grid scale in th-direction

 phi=(/(0.d0 + dph*(i-1),i=1,sz_ph)/)
 thi=(/(-PI/2 + dth*(i-1),i=1,sz_th)/)
 ri=(/(grid_Rin + dr*(i-1),i=1,sz_r)/)

 ! ;transform image coord of the point to Carrington coord.
  call img2car_coord (cobs, xi, yi, zi, x0, y0, z0)

!print*,'[lon, lat] of obs=',cobs(1),cobs(2)
! print*,'----Spherical grid--------------'

! ;calculate the direction vector of the ray
 phc=cobs(1)/radeg
 thc=cobs(2)/radeg
 ni=(/cos(phc)*cos(thc),sin(phc)*cos(thc), sin(thc)/)
!;a,b,c are unit vector of LOS
 a=ni(1)
 b=ni(2)
 c=ni(3)
!print*,'unit vector for Ray (a,b,c)=',a,b,c

!$OMP PARALLEL SECTIONS DEFAULT(PRIVATE) &
!$OMP        SHARED(phi,thi,ri,a,b,c,x0,y0,z0,cross,rt,tht,pht,&
!$OMP        n_ri,n_thi,n_phi, r_ri,th_ri,ph_ri, r_thi,th_thi,ph_thi, r_phi,th_phi,ph_phi)

!$OMP SECTION 
! print*,'Section-1: Thread ID=',omp_get_thread_num(),'Used thread num=',omp_get_num_threads()
! print*,'Max active nested level=',omp_get_max_active_levels()

!; step-1: calculate the cross-points of the ray with spherical shells at r=ri
!; equation of a straight line with parameter-t
!; x=x0+at, y=y0+bt, z=z0+ct
!; (x,y,z)-> (r,th,ph):
!;   x=r*cos(th)*cos(ph)
!;   y=r*cos(th)*sin(ph)
!;   z=r*sin(th)
!; equation of cross-point: (a^2+b^2+c^2)t^2 +2(a*x0+b*y0+c*z0)t + (x0^2+y0^2+z0^2-ri^2)=0
 Af1=a*a+b*b+c*c
 Bf1=2*(a*x0+b*y0+c*z0)
 Cf1=x0*x0+y0*y0+z0*z0-ri*ri  !;Cf is a vector
 dlt1=Bf1*Bf1-4*Af1*Cf1
 nw=0
 DO i=1, sz_r
    IF (dlt1(i) >= 0) THEN 
        nw=nw+1
        w1(nw)=i
        dltw1(nw)=dlt1(i)
    ENDIF
 ENDDO

 IF (nw <= 0) THEN
   STOP "No crossing points-- Error!"
 ELSE 
   tn1=(-Bf1-sqrt(dltw1))/2/Af1
   tp1=(-Bf1+sqrt(dltw1))/2/Af1
 ENDIF

 IF ((nw == 1) .and. (dltw1(1) == 0)) THEN
   print*, 'ri=',ri(w1(1)),' Ray is tangient to the outer boundary!'
   if (ri(w1(1)) /= grid_Rout) then
     STOP "ri/=Rout, Error!"
   endif
   cross=0
   print*, 'Discard this ray by setting cross=0'
!  RETURN
 ELSE IF ((nw > 1) .and. (dltw1(1) ==0)) THEN ! ;cross at least two shells and
 !  r_ri(nw)=ri(w1(1))                         !;tangient to the shell with smallest ri
    print*,'*****tangient to the shell with smallest ri  ********'
    STOP "This special case needs to be treated!"
 ELSE ! ;with at least two cross-points
    n_ri=2*nw
    do i=1,nw
        r_ri(i)   =ri(w1(nw+1-i))
        r_ri(nw+i)=ri(w1(i))
        t1(i)   =tn1(nw+1-i)
        t1(nw+i)=tp1(i)
     enddo
 ENDIF
!     ;Sorting ri[w] is required for Step-3 (to determin entry and exit points for ph_ri
 DO i=1, n_ri
    th_ri(i)=asin((z0+c*t1(i))/r_ri(i))   !;th=[-pi/2, pi/2]
    ph_ri(i)=atan2(y0+b*t1(i), x0+a*t1(i)) !;ph=(-pi, pi]
    ph_ri(i)=mod(ph_ri(i)+2*PI, 2*PI)      !;ph=[0, 2*pi]
 ENDDO

!print*,'----------cross-points for ri:  n_ri=',n_ri
!print*,'r=',r_ri(1:n_ri)
!print*,'th=',th_ri(1:n_ri)*radeg
!print*,'ph=',ph_ri(1:n_ri)*radeg


!$OMP SECTION 
! print*,'Section-2: Thread ID=',omp_get_thread_num(),'Used thread num=',omp_get_num_threads()
! print*,'Max active nested level=',omp_get_max_active_levels()

!;Step-2: calculate the cross-points of the ray with cone surface with
!;  latitude th=thi
!; equation of ray:
!;  (a^2+b^2-c^2*ctg(thi)^2)t^2 +(2a*x0+2b*y0-2c*z0*ctg(thi)^2)t +(x0^2+y0^2-z0^2*ctg(thi)^2=0
 ctg2=1./tan(thi)**2   !; Note: ctg2, Af2,Bf2,Cf2, are all vector
 Af2=a*a+b*b-c*c*ctg2
 Bf2=2*a*x0+2*b*y0-2*c*z0*ctg2
 Cf2=x0*x0+y0*y0-z0*z0*ctg2
 dlt2=Bf2*Bf2-4*Af2*Cf2
 nw=0
 nw1=0
 DO i=1, sz_th !;determine indice of thi-cones crossed by the line
  if (abs(Af2(i)) < epsn) then
       print*,'Ray is parallel to cone surface!'     
       nw1=nw1+1
       w21(nw1)=i
!      STOP "A special treatment is required ****************"
  else !af2 /=0 
    IF (dlt2(i) > 0) THEN
        nw=nw+1
        w2(nw)=i
        dltw2(nw)=dlt2(i)
    ENDIF
    IF (dlt2(i) == 0) THEN
!       nw0=nw0+1
!       w0(nw0)=i
     print*,'Ray is tangient to the cone surface!'
      STOP "A special treatment is required!--------------"
    ENDIF
  endif
 ENDDO

 n_thi =0
 IF (nw > 0) THEN  !for dlt>0 only
   DO i=1,nw
      tn2(i)=(-Bf2(w2(i))-sqrt(dltw2(i)))/2/Af2(w2(i))
      tp2(i)=(-Bf2(w2(i))+sqrt(dltw2(i)))/2/Af2(w2(i))
      th_thi(i)    = thi(w2(i))
      th_thi(nw+i) = thi(w2(i))  
      t2(i)    = tn2(i)
      t2(nw+i) = tp2(i)
   ENDDO  
   thi_cross=1
   n_thi=2*nw
 ENDIF

 IF (nw1 > 0) THEN
    DO i=1,nw1
       th_thi(n_thi+i) =thi(w21(i))
       t2(n_thi+i)     = -Cf2(w21(i))/Bf2(w21(i)) 
    ENDDO
    thi_cross=1
    n_thi=n_thi +nw1
 ENDIF        

 IF (nw+nw1 ==0) THEN 
   print*,'No crossing points with thi-cones!'
   thi_cross=0
 ENDIF

!;excluding the unsatisfied solutions due to ctg^2(th)
!;From sin(thi)=(z0+ct)/r we obtain the restriction for solutions ti:
!; sin(thi)*(z0+c*ti)>=0 related with -/+sqrt(dlt)
! combined with condition for r_thi within Rin<= r <= Rout
  IF (thi_cross == 1) THEN
     r_thi=sqrt((x0+a*t2)**2+(y0+b*t2)**2+(z0+c*t2)**2)
     ph_thi=atan(y0+b*t2, x0+a*t2)
     ph_thi=mod((ph_thi+2*PI), 2*PI)
     nw=0
     DO i=1, n_thi
         if (sin(th_thi(i))*(z0+c*t2(i)) >=0 .and. &
             (r_thi(i) >= grid_Rin .and. r_thi(i) <= grid_Rout)) then
           nw=nw+1
           r_thi(nw) =r_thi(i)
           th_thi(nw)=th_thi(i)
           ph_thi(nw)=ph_thi(i)
         endif
     ENDDO
      if (nw >0) then
          n_thi=nw
      else
!        print*,'no crossing points for thi meeting the following conditions:'
!        print*,'sin(thi)*(z0+c*ti)>=0 and Rin<= ri <= Rout'
         thi_cross=0
         n_thi =0
      endif
   ENDIF

!  ;Check if thi(:) includes thi=0 or very close to 0. (singular point for ctg2()
  DO i=1, sz_th
     if (abs(thi(i)) < epsn) then
         th_thi(n_thi+1)=thi(i)
         t0=-z0/c
         r0_thi=sqrt((x0+a*t0)**2+(y0+b*t0)**2)
         if (r0_thi >= grid_Rin .and. r0_thi <= grid_Rout) then
             r_thi(n_thi+1) = r0_thi
             ph_thi(n_thi+1)=mod((atan(y0+b*t0, x0+a*t0)+2*PI), 2*PI)
             n_thi= n_thi+1
             thi_singular=1
         endif
     endif
  ENDDO


!print*,'----------------cross-points for thi: n_thi=',n_thi
!print*,'r=',r_thi(1:n_thi)
!print*,'th=',th_thi(1:n_thi)*radeg
!print*,'ph=',ph_thi(1:n_thi)*radeg
!$OMP END PARALLEL SECTIONS

! Caution: step-3 cannot be included in SECTIONS because the result from
!  steps-1 and -2 (n_ri and n_thi) is used
!!!!
!;Step-3: Calculate the cross-points of the Ray with longitudinal plane
!; with ph=phi (large circle phi or phi+Pi, where
!; phi should satisfy the condition, abs(ph(in)-ph(out))<=Pi
!print*,'ph for ri-entry=',radeg*ph_ri(1)
!print*,'ph for ri-exit=',radeg*ph_ri(n_ri)
 phi_cross=1
 ph_min=minval((/ph_ri(1),ph_ri(n_ri)/))
 ph_max=maxval((/ph_ri(1),ph_ri(n_ri)/))
 nw=0
 IF (ph_max-ph_min <= PI) THEN 
    DO i=1,sz_ph
       if (phi(i) >= ph_min .and. phi(i) <= ph_max) then
           nw=nw+1
           ph_phi(nw)=phi(i)
       endif
    ENDDO
    n_phi=nw
    if (n_phi ==0) then
        print*,'case 1: not crossing phi ! '
        phi_cross=0
     !; this case did happen when sampled images have large size
    endif
 ELSE !;ph_max-ph_min > !pi
     nw1=0
     nw2=0
     DO i=1, sz_ph
        if (phi(i) >= ph_max .and. phi(i) < 2*PI) then
            nw1=nw1+1
            w31(nw1)=i
        endif
        if (phi(i) >= 0 .and. phi(i) <= ph_min) then
            nw2=nw2+1
            w32(nw2)=i
        endif
     ENDDO
     IF (nw1 >=1 .or. nw2 >=1) THEN
        do i=1, nw1 
            ph_phi(i)= phi(w31(i))
        enddo
        do i=1, nw2
            ph_phi(nw1+i)= phi(w32(i)) + 2*PI
        enddo
        n_phi = nw1+nw2
     ELSE
        STOP "case 2: not crossing phi ! Stop"
     ENDIF

     !;set ph_ri,ph_thi in [0,ph_min] to be in [0,ph_min]+2*!pi
     DO i=1,n_ri
        if (ph_ri(i) <= ph_min .and. ph_ri(i) >=0) then
            ph_ri(i)=ph_ri(i) + 2*PI
        endif
     ENDDO

     DO i=1,n_thi
        if (ph_thi(i) <= ph_min .and. ph_thi(i) >=0) then
            ph_thi(i)=ph_thi(i) + 2*PI
        endif
     ENDDO
 ENDIF

 IF (phi_cross == 1) THEN
     DO i=1, n_phi
        ti      = (y0-x0*tan(ph_phi(i)))/(a*tan(ph_phi(i))-b)
        r_phi(i)=sqrt((x0+a*ti)**2+(y0+b*ti)**2+(z0+c*ti)**2)
        th_phi(i)=asin((z0+c*ti)/r_phi(i))
     ENDDO
 ENDIF 


!;sort cross-points for ri, thi, phi in order of longitude-ph
  pnum=n_ri + n_thi + n_phi
! print*,'Total crossing point number =',pnum

 ALLOCATE(r_all(pnum),th_all(pnum), ph_all(pnum), order(pnum),re(pnum),the(pnum),phe(pnum),&
     ccell(3,pnum-1),xe(pnum),ye(pnum),ze(pnum),STAT = AllocateStatus)
 IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

  IF (pnum >= 1) THEN
      DO i=1, n_ri
         r_all(i)  =r_ri(i)
         th_all(i) =th_ri(i)
         ph_all(i) =ph_ri(i)
      ENDDO
      DO i=1, n_thi
         r_all(n_ri+i)  =r_thi(i)
         th_all(n_ri+i) =th_thi(i)
         ph_all(n_ri+i) =ph_thi(i)
      ENDDO
      DO i=1, n_phi
         r_all(n_ri+n_thi+i)  =r_phi(i)
         th_all(n_ri+n_thi+i) =th_phi(i)
         ph_all(n_ri+n_thi+i) =ph_phi(i)
      ENDDO
  ELSE
    STOP "---No crossing points, errors --Stop!"
  ENDIF

! sort ph_all(:) and return the indice in increasing order
  call mrgrnk(ph_all, order)
  re  = r_all(order)
  the = th_all(order)
  phe = ph_all(order)

! convert phe[:] in rang of [0, 2*PI) after sorting the points
   phe = mod (phe, 2*PI)


!; determine the indice for each crossed cell
 call sp2car(phe,the,re, xe,ye,ze)

IF (1) THEN
! measure the length of ray crossing the spheric box
 slen=0.
 DO i=1,pnum-1
    slen=slen+sqrt((xe(i+1)-xe(i))**2+(ye(i+1)-ye(i))**2+(ze(i+1)-ze(i))**2)
 ENDDO

  rt=(/r_ri(1), r_ri(n_ri)/)
 tht=(/th_ri(1), th_ri(n_ri)/)
 pht=(/ph_ri(1), ph_ri(n_ri)/)

!print*,'-----Summed length from segments=',slen
 !calculate ray length from two crossing points
  call sp2car(pht,tht,rt, xt,yt,zt)
  tlen=sqrt((xt(2)-xt(1))**2+(yt(2)-yt(1))**2+(zt(2)-zt(1))**2)
!print*,'--------Ray length from two ends=',tlen
ENDIF

!print*,'cell index [ir, ith, iph]='
 DO i=1,pnum-1 
    xm=(xe(i)+xe(i+1))/2       !;calculate mid-point position of the segment inside a cell
    ym=(ye(i)+ye(i+1))/2       !; and assume this point should be inside the cell
    zm=(ze(i)+ze(i+1))/2
    call car2sp(xm,ym,zm, phm,thm,rm)
    phm=mod(phm+2*PI, 2*PI)          !;ph in r(0,2*pi)

    iph=int((phm-0.)/dph)   
    ith=int((thm+PI/2)/dth) 
    ir =int((rm-grid_Rin)/dr) 
   ccell(1,i)=mod(iph, (sz_ph-1))     !;keep the cell index in [0,sz_ph-2]
   ccell(2,i)=ith                          !;index in [0, sz_th-1]
   ccell(3,i)=ir                           !;index in [0,sz_r-1]
!   write(*,'(1x, i5, 3(1x,G15.7))') i,ccell(3,i),ccell(2,i),ccell(1,i)
 ENDDO
END SUBROUTINE sph_cross_points


FUNCTION  sph_ray_weights(phe,the,re,ccell,ni) RESULT(wt)
!;+purpose: calculate weights (contribution func to line integral) for a ray
!;passing through a spheric grid box
!; input:
!;  phe,the,re[N] - Spheric coords. of N crossing points for a ray
!; ccell[3,N-1] - indice of N-1 cells which are passed through
!;    Note: index of ccell for iph,ith,ir starts from 0
!; return: weights wt[sz_ph,sz_th,sz_r] for grid contributed by the ray in unit of Rsun
! note: as index() is an intrinsic func, so replace variable 'index' by 'index1d'
USE InternPro, only: PI
IMPLICIT NONE
REAL (kind=double), INTENT(IN), DIMENSION(:) :: phe,the,re
INTEGER, INTENT(IN), DIMENSION(3,size(re)-1) :: ccell
REAL (kind=double), INTENT(IN), DIMENSION(3) :: ni
REAL (kind=double), DIMENSION(1:sz_ph,1:sz_th,1:sz_r) :: wt
INTEGER :: i,j,k,i0,j0,k0, m, index1d(3)
REAL (kind=double), DIMENSION(2,2,2) :: f
REAL (kind=double), DIMENSION(2) :: pht,tht,rt
REAL,PARAMETER :: epsn=1e-10

 wt=0.d0
!note: initial values must be set because wt=wt+f is used!
 DO m=1, size(re)-1
       index1d=ccell(:,m)
    pht=phe(m:m+1)
    tht=the(m:m+1)
    rt=re(m:m+1)
    if (abs(pht(2)) <epsn) pht(2)= 2*PI
     f=sph_cell_weights(index1d,pht,tht,rt,ni)
    DO 100 k=0,1
    DO 100 j=0,1
    DO 100 i=0,1
        i0=index1d(1)+i +1
        j0=index1d(2)+j +1
        k0=index1d(3)+k +1
!  Note: values of index1d(:) are in range [0, sz-1] 
        wt(i0,j0,k0)=wt(i0,j0,k0)+f(i+1,j+1,k+1)
    100 CONTINUE
 ENDDO
!print*,'==Total(wt)=',SUM(wt)
return
END FUNCTION  sph_ray_weights

FUNCTION sph_cell_weights (index1d, pht, tht, rt, ni) RESULT(wt)
!;+purpose: return weights calculated for a spheric cell with given index (i,j,k)
!; input:
!;  index[3] - pixel index of a given cell. Note values of indice are in range of
!;            index(1)-0:sz_ph-1, index(2)-0:sz_th-1, index(3)-0:sz_r-1
!;  pht[2],tht[2],rt[2] - entry[i=1] and exit[i=2] coords. of crossing points of 
!         a ray through a cell
!; output:
!;   return weights f[2,2,2] for 8 vertices with i=0,1;j=0,1;k=0,1
!;-
USE InternPro
IMPLICIT NONE
INTEGER, INTENT(IN), DIMENSION(3) :: index1d
REAL (kind=double), INTENT(IN), DIMENSION(2) :: pht,tht,rt 
REAL (kind=double), INTENT(IN), DIMENSION(3) :: ni
REAL (kind=double), DIMENSION(2,2,2) :: wt
REAL (kind=double), DIMENSION(2) :: xt,yt,zt
REAL (kind=double) :: ds, ds_sp, dr,dth,dph, dxt,dyt,dzt,tlen, nx,ny,nz,&
                      p0,t0,r0, ftt
INTEGER :: i,j,k, nds, tnum, n
REAL (kind=double), DIMENSION(:), ALLOCATABLE :: p,t,r, ids, xi,yi,zi, &
       phi,thi,ri, fp,ft,fr, f
REAL  (kind=double) :: ndt
integer, dimension(:), allocatable:: ifp,ift,ifr
real, parameter :: epsn=1e-10


 dr=(grid_Rout-grid_Rin)/dble(sz_r-1)      !;grid scale in r-direction
 dth=PI/dble(sz_th-1)                     !;grid scale in th-direction
 dph=2*PI/dble(sz_ph-1)                   !;grid scale in th-direction

  ds=minval((/dr, grid_Rin*dth, grid_Rin*dph/))
  nds=100  !50
  ds_sp=ds/nds  !;assume step length for line integral
!;Note: tests show that integral resolution (ie. choice od nds) has only smaller
!; effect on wt(i,j,k) (0.2% - 2%), while has nearly no effect on total(wt)
!print*,'ds_sp=',ds_sp


 call sp2car(pht,tht,rt, xt,yt,zt)
 dxt=xt(2)-xt(1)
 dyt=yt(2)-yt(1)
 dzt=zt(2)-zt(1)
 tlen=sqrt(dxt*dxt+dyt*dyt+dzt*dzt)



! print*,'cell seg-len=',tlen
!;calculate the direct vector [nx,ny,nz]
 nx=dxt/tlen
 ny=dyt/tlen
 nz=dzt/tlen

!;estimate the step number for line integral
 tnum=int(tlen/ds_sp)


!print*,' Tnum= nds=',tnum, nds
 allocate(ids(tnum),xi(tnum),yi(tnum),zi(tnum),phi(tnum),thi(tnum),ri(tnum),&
            p(tnum+2),t(tnum+2),r(tnum+2),fp(tnum+2),ft(tnum+2),fr(tnum+2),&
            f(tnum+2))

 IF (tnum == 0) THEN    !;cell segment < step
    p=pht
    t=tht
    r=rt
 ELSE     !  ;cell segment >= 1 step
    ids=(/(i*ds_sp,i=1,tnum)/)  
    xi=xt(1)+ids*nx
    yi=yt(1)+ids*ny
    zi=zt(1)+ids*nz
   call car2sp(xi,yi,zi, phi,thi,ri)
    phi=mod((phi+2*PI), 2*PI)       !;ph in [0,2*!pi]
    p(1)=pht(1)
    t(1)=tht(1)
    r(1)= rt(1)
    p(2:tnum+1)=phi
    t(2:tnum+1)=thi
    r(2:tnum+1)=ri
    p(tnum+2)=pht(2)
    t(tnum+2)=tht(2)
    r(tnum+2)= rt(2)
 ENDIF

 p0= 0.d0    + index1d(1)*dph
 t0= -PI/2   + index1d(2)*dth
 r0= grid_Rin+ index1d(3)*dr
 fp=(p-p0)/dph
 ft=(t-t0)/dth
 fr=(r-r0)/dr

!; check if 0<=(fp,ft,fr)<=1 and force it to be
!;  accurate strictly at boundary f=0 or 1
 if (abs(fp(1)) < epsn) fp(1)=0.d0
 if (abs(fp(tnum+2)-1)< epsn) fp(tnum+2)=1.d0

 if (abs(ft(1)) < epsn) ft(1)=0.d0
 if (abs(ft(1)-1) < epsn) ft(1)=1.d0
 if (abs(ft(tnum+2))< epsn) ft(tnum+2)=0.d0
 if (abs(ft(tnum+2)-1)< epsn) ft(tnum+2)=1.d0

 if (abs(fr(1)) < epsn) fr(1)=0.d0
 if (abs(fr(1)-1)< epsn) fr(1)=1.d0
 if (abs(fr(tnum+2)) < epsn) fr(tnum+2)=0.d0
 if (abs(fr(tnum+2)-1)< epsn) fr(tnum+2)=1.d0

  wt=0.d0
  DO 100 k=0,1 
  DO 100 j=0,1 
  DO 100 i=0,1 
      f=(-1)**(i+j+k)*(1-i-fp)*(1-j-ft)*(1-k-fr)
      !;calculate the line integral over s=0 to tlen
     if (tnum==0) then
        ftt=tlen*sum(f)/2  !;integral path includes only entry and exit points,
                            !; where value of f is taken as the mean of f(0) and f(1) in integeral
     else if (tnum == 1) then  !; case of integral includes two parts 
        ftt=(tlen-ds_sp)*sum(f(2:3))/2 +&
              ds_sp*sum(f(1:2))/2
     else     !tnum ge 2 for case that integral includes tnum+1 parts
       n=tnum+2
        ftt=(tlen-ds_sp*tnum)*sum(f(n-1:n))/2. +&    !for integral of segment-i, f(i) is taken
             ds_sp*(f(1)/2.+f(n-1)/2.+sum(f(2:n-2)))  ! as f(i)=(f(i)+f(i+1))/2, where i=1, n-1
     endif
     wt(i+1,j+1,k+1)=ftt
  100 CONTINUE
!print*,'total weight=',sum(wt)

 return
END FUNCTION sph_cell_weights


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
!return the index (i,j) of a non-zero element with index (in SP%A) for
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
 use OMP_LIB
IMPLICIT NONE
TYPE (SparseMat), INTENT(IN) :: A
INTEGER, INTENT(IN) :: n
REAL, INTENT(IN), DIMENSION(n) :: x
REAL, DIMENSION(A%m) :: y
INTEGER :: irow,i, istt,iend
INTEGER :: CHUNK 

 CHUNK =100
 Y=0.0 
!$OMP PARALLEL SHARED(A,X,Y,n) PRIVATE(irow,istt,iend,i)
!$OMP DO REDUCTION(+:Y) SCHEDULE(DYNAMIC,CHUNK)
DO irow=1, A%m
   istt=A%ifst(irow)
   iend=A%ifst(irow+1)-1
   DO i=istt, iend
      IF (A%icol(i)<=n) THEN
         y(irow)=y(irow)+A%A(i)*x(A%icol(i))
      ENDIF
   ENDDO
ENDDO
!$OMP END DO 
!$OMP END PARALLEL 
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
 use OMP_LIB
IMPLICIT NONE
INTEGER, INTENT(IN) :: m
REAL, INTENT(IN), DIMENSION(m):: x
TYPE (SparseMat), INTENT(IN) :: A
REAL, DIMENSION(A%n) :: y
INTEGER :: i,j, istt, iend, irow
INTEGER :: CHUNK

CHUNK =100
y=0.0
irow=min(A%m,m)
!$OMP PARALLEL SHARED(A,X,Y,irow) PRIVATE(istt,iend,i,j)
!$OMP DO REDUCTION(+:Y) SCHEDULE(DYNAMIC,CHUNK)
DO i=1, irow
   istt=A%ifst(i)
   iend=A%ifst(i+1)-1
   DO j=istt, iend
      y(A%icol(j))=y(A%icol(j)) + x(i)*A%A(j)
   ENDDO
ENDDO
!$OMP END DO 
!$OMP END PARALLEL
RETURN
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
! keyword: tol -the desired convergence tolerance
!       define tol=|A*x-y|/|y|
!     sav_res - save the middle results for resume. If "resume_iter.sav" exists,
!            save a copy (resume_iter.sav_copy) first.
!     resume - to resume the calculation if set
!     iterr_sav - save data and regularization errors for each step if set
!-
   use OMP_LIB 
IMPLICIT NONE
INTEGER, PARAMETER ::  itmax = 10000
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

!print*,'bvect min max=',minval(bvec),maxval(bvec)
!print*,'btol=',btol

!$OMP PARALLEL 
  !$OMP MASTER
  !$  print*,'Total used OMP thread Num=',omp_get_num_threads()
  !$OMP END MASTER

  !$OMP WORKSHARE
 xk=0. !set initial solution as zero
 rk=bvec !r0=b-A*x0, set initial residual
 pk=rk
  !$OMP END WORKSHARE
!$OMP END PARALLEL 

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
!$OMP PARALLEL WORKSHARE        
       Apk=Apk+mu*Rpk            !Apk=A^T*A*pk + mu*R^T*R*pk
!$OMP END PARALLEL WORKSHARE       
    ENDIF
    alphak=rrk/dot_product(pk,Apk)
!$OMP PARALLEL WORKSHARE    
    xk1=xk +alphak*pk
    rk1=rk -alphak*Apk
!$OMP END PARALLEL WORKSHARE
    rerr=sqrt(sum(rk1**2))/btol  !relative error
  IF (mod(k,100)==0) print*,'k=',k,' res=',rerr, ' (tol=',tol,')'
    IF (rerr <= tol) EXIT
    rrk1=dot_product(rk1,rk1)
    betak=rrk1/rrk
!$OMP PARALLEL WORKSHARE    
    pk1=rk1 +betak*pk
    rk=rk1
    pk=pk1
    xk=xk1
    rrk=rrk1
!$OMP END PARALLEL WORKSHARE     
    k=k + 1
  IF (mod(k,100)==0)  THEN  
     IF (sav_res==1) THEN !save the middle result for resume
        OPEN(1,FILE='resume_iter.dat',FORM='unformatted',ACTION='write',STATUS='replace')
        WRITE(1) k,rrk
        WRITE(1) rk,pk,xk
        CLOSE(1)
     ENDIF
     IF (sav_iter==1) THEN !save data and reg. errors for each iteration step
    !  calculate data and regularization errors
!$OMP PARALLEL WORKSHARE       
       WHERE(xk1 >0) !set xknn including non-negative elements in xk1(:)
          xknn=xk1
       ELSEWHERE
          xknn=0.
       END WHERE
!$OMP END PARALLEL WORKSHARE
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
  ENDIF  
 ENDDO
!$OMP PARALLEL WORKSHARE
  xsolu=xk1
!$OMP END PARALLEL WORKSHARE  
RETURN
END FUNCTION ConjugateGrad_SP


FUNCTION sph_Reg2nd2df_SP (php,thp,rp, nr,pcos2) RESULT(R)
!;+ purpose: Second-order smoothing using central differetiate of 2nd-order accuracy
!; Input:
!;   php,thp,rp[sz1_ph,sz_th,sz_r]-- 3D spheric grid
!; keyword: nr=nr
!;     nr[0],nr[1] gives the dependence of normalized N(r) on radial distance nr[0]
!;     If set, add a radial weight w=1/N(r) to the regularization term Rw^T#Rw
!; return: R -mm x nn sparse matrix for regularization, where
!;          nn=ndis and mm=3*nn
!; add weight=1/nr for each grid point remove r-dependence effect of regularization
USE InternPro
IMPLICIT NONE
REAL (kind=double), INTENT(IN), DIMENSION(:,:,:) :: php,thp,rp
REAL (kind=double), INTENT(IN), DIMENSION(:,:) :: nr
INTEGER, INTENT(IN) :: pcos2
TYPE (SparseMat) :: R
REAL (kind=double) :: dph,dth,dr
INTEGER :: sp,st,sr,mm,nn,nnz, ii, IA(3*size(php)+1), i,j,k, iw1,iw2,iw3,nja
TYPE (dataRys) :: Npi, Nti, Nri, mRys(3*size(php))
REAL (kind=double),DIMENSION(:,:,:), ALLOCATABLE :: nrw

 allocate(nrw, mold=rp)  ! nrw is allocated with the bound of rp

   dr=(grid_Rout-grid_Rin)/dble(sz_r-1)  !grid scale in r-direction
  dth=PI/dble(sz_th-1)                   !grid scale in th-direction
  dph=2*PI/dble(sz_ph-1)                 !grid scale in th-direction
 print*,'Regularization term is chosen as'
 print*,'d^2/dph^2 =1/(r*cos[t]*dph)^2 *(N(i-1)-2*N(i)+N(i+1))'
 print*,'d^2/dth^2 =1/(r*dth)^2*(N(j-1)-2*N(j)+N(j+1))'
 print*,'d^2/dr^2 =1/dr^2 *(N(k-1)-2*N(k)+N(k+1))'

 sp=SIZE(php, DIM=1)
 st=SIZE(php, DIM=2)
 sr=SIZE(php, DIM=3)
 nn =sp*st*sr       !;number of unknown elements in solution vector, X
! mm =3*nn           ! initial row number of regularization matrix, R
print*,'sp,sr,sr=',sp,st,sr

  call Allocate_Rys(Npi, 3)
  call Allocate_Rys(Nti, 3)
  call Allocate_Rys(Nri, 3)

 IF (nr(1,1)>0) THEN
    print*,'--- Regularization with radial weighting w=1/Nr is applied!'
    do k=1,sr
    do j=1,st
       call interp(nr(:,1),nr(:,2), rp(:,j,k),nrw(:,j,k))
    enddo
    enddo
 ENDIF

 IF (pcos2 ==1) THEN
    print*,' a factor 1/cos(th)^2 is added to the term d^2/dph^2!----'
 ENDIF

  ii=1      ! row number of R_sparse matrix
  IA(1)=1  ! index in R%aa for 1st non-0 element in 1st row of R_sparse

  DO k=0, sr-1
  DO j=1, st-2     !;excluding voxels with th=+/-pi/2 (singular points)
  DO i=0, sp-1
                   !;Reg for longitude: d^2/dp^2 = (N(i-1)-2*N(i)+N(i+1))
    if (i == 0) then
       Npi%n=3
       iw1=sp+ j*sp+k*sp*st
       iw2= 1+ j*sp+k*sp*st
       iw3= 2+ j*sp+k*sp*st
       Npi%ja(:)=(/iw2,iw3,iw1/)
       Npi%aa(:)=(/-2, 1, 1/)
    else if (i == sp-1) then
       Npi%n=3
       iw1=(i+0)+j*sp+k*sp*st
       iw2=(i+1)+j*sp+k*sp*st
       iw3= 1   +j*sp+k*sp*st
       Npi%ja(:)=(/iw3,iw1,iw2/)
       Npi%aa(:)=(/1, 1, -2/)
    else
       Npi%n=3
       iw1=(i+0)+j*sp+k*sp*st
       iw2=(i+1)+j*sp+k*sp*st
       iw3=(i+2)+j*sp+k*sp*st
       Npi%ja(:)=(/iw1,iw2,iw3/)
       Npi%aa(:)=(/1, -2, 1/)
    endif
    Npi%aa=Npi%aa/(rp(i+1,j+1,k+1)*dph)**2
    if (pcos2==1) then
        Npi%aa=Npi%aa/cos(thp(i+1,j+1,k+1))**2
    endif

        !;Reg for latitude: for d^2/dt^2 = (N(j-1)-2*N(j)+N(j+1))
     Nti%n=3
     iw1=(i+1)+(j-1)*sp+k*sp*st
     iw2=(i+1)+j*sp+k*sp*st
     iw3=(i+1)+(j+1)*sp+k*sp*st
     Nti%ja(:)=(/iw1,iw2,iw3/)
     Nti%aa(:)=(/1, -2, 1/)
     Nti%aa=Nti%aa/(rp(i+1,j+1,k+1)*dth)**2

!        ;Reg for radial: d^2/dr^2 = (N(k-1)-2*N(k)+N(k+1))
    if (k == 0) then
       Nri%n=2      !;right diff: (N(k+1)-N(k))/dr
       iw1=(i+1)+j*sp+k*sp*st
       iw2=(i+1)+j*sp+(k+1)*sp*st
       Nri%ja(:)=(/iw1,iw2/)
       Nri%aa(:)=(/-1, 1/)/dr
    else if (k == sr-1) then
       Nri%n=2      !;left diff: (N(k)-N(k-1))/dr
       iw1=(i+1)+j*sp+(k-1)*sp*st
       iw2=(i+1)+j*sp+k*sp*st
       Nri%ja(:)=(/iw1,iw2/)
       Nri%aa(:)=(/-1, 1/)/dr
    else
       Nri%n=3
       iw1=(i+1)+j*sp+(k-1)*sp*st
       iw2=(i+1)+j*sp+k*sp*st
       iw3=(i+1)+j*sp+(k+1)*sp*st
       Nri%ja(:)=(/iw1,iw2,iw3/)
       Nri%aa(:)=(/1, -2, 1/)/dr**2
    endif


   if (nr(1,1)>0) then
      Npi%aa=Npi%aa/nrw(i+1,j+1,k+1)
      Nti%aa=Nti%aa/nrw(i+1,j+1,k+1)
      Nri%aa=Nri%aa/nrw(i+1,j+1,k+1)
   endif


      nja=Npi%n
      call Allocate_Rys(mRys(ii),nja)
      mRys(ii)%ja = Npi%ja(1:nja)
      mRys(ii)%aa = Npi%aa(1:nja)
      IA(ii+1) = IA(ii)+nja
      ii=ii+1

      nja=Nti%n
      call Allocate_Rys(mRys(ii),nja)
      mRys(ii)%ja = Nti%ja(1:nja)
      mRys(ii)%aa = Nti%aa(1:nja)
      IA(ii+1) = IA(ii)+nja
      ii=ii+1

      nja=Nri%n
      call Allocate_Rys(mRys(ii),nja)
      mRys(ii)%ja = Nri%ja(1:nja)
      mRys(ii)%aa = Nri%aa(1:nja)
      IA(ii+1) = IA(ii)+nja
      ii=ii+1
  ENDDO
  ENDDO
  ENDDO

 mm=ii-1
 nnz = IA(mm+1)-1

  print*,'ii=',ii,' mm=',mm, '  Should be: ii=mm+1'
  print*, 'nnz=',nnz

 call Allocate_SP(R, mm, nn, nnz)
 R%ifst = IA
 DO i=1,mm
    R%icol(IA(i):IA(i+1)-1) = mRys(i)%ja
    R%A(IA(i):IA(i+1)-1) = mRys(i)%aa
    call DeAllocate_Rys(mRys(i))
 ENDDO
 call DeAllocate_Rys(Npi)
 call DeAllocate_Rys(Nti)
RETURN
END FUNCTION sph_Reg2nd2df_SP


FUNCTION sph_Reg2nd2rtp_SP (php,thp,rp, nr,pcos2) RESULT(R)
!;+ purpose: Second-order smoothing using central differetiate of 2nd-order accuracy
!;     with the form, |RX|^2=Sum( (d^2/dr^2)^2 +(d^2/dth^2)^2 +(d^2/dph^2)^2)
!; Input:
!;   php,thp,rp[sz1_ph,sz_th,sz_r]-- 3D spheric grid
!; keyword: nr=nr
!;     nr[0],nr[1] gives the dependence of normalized N(r) on radial distance nr[0]
!;     If set, add a radial weight w=1/N(r) to the regularization term Rw^T#Rw
!;   pcos2=pcos2
!;      If set, add a coefficient 1/cos^2(phi) to d^2/dph^2)^2, to increase smooth
!;        towards the polar regions
!; return: R -mm x nn sparse matrix for regularization, where
!;          nn=ndis and mm=3*nn
!; add weight=1/nr for each grid point remove r-dependence effect of regularization
USE InternPro
IMPLICIT NONE
REAL (kind=double), INTENT(IN), DIMENSION(:,:,:) :: php,thp,rp
REAL (kind=double), INTENT(IN), DIMENSION(:,:) :: nr
INTEGER, INTENT(IN) :: pcos2
TYPE (SparseMat) :: R
!REAL (kind=double) :: dph,dth,dr
INTEGER :: sp,st,sr,mm,nn,nnz, ii, IA(3*size(php)+1), i,j,k, iw1,iw2,iw3,nja
TYPE (dataRys) :: Npi, Nti, Nri, mRys(3*size(php))
REAL (kind=double),DIMENSION(:,:,:), ALLOCATABLE :: nrw

 allocate(nrw, mold=rp)  ! nrw is allocated with the bound of rp

!  dr=(grid_Rout-grid_Rin)/dble(sz_r-1)  !grid scale in r-direction
! dth=PI/dble(sz_th-1)                   !grid scale in th-direction
! dph=2*PI/dble(sz_ph-1)                 !grid scale in th-direction
 print*,'Regularization term is chosen as'
 print*,'d^2/dph^2 =(N(i-1)-2*N(i)+N(i+1)), d^2/dth^2 =(N(j-1)-2*N(j)+N(j+1)), d^2/dr^2  =(N(k-1)-2*N(k)+N(k+1))'

 sp=SIZE(php, DIM=1)
 st=SIZE(php, DIM=2)
 sr=SIZE(php, DIM=3)
 nn =sp*st*sr       !;number of unknown elements in solution vector, X
! mm =2*nn           ! initial row number of regularization matrix, R
print*,'sp,sr,sr=',sp,st,sr

  call Allocate_Rys(Npi, 3)
  call Allocate_Rys(Nti, 3)
  call Allocate_Rys(Nri, 3)

 IF (nr(1,1)>0) THEN
    print*,'--- Regularization with radial weighting w=1/Nr is applied!'
    do k=1,sr
    do j=1,st
       call interp(nr(:,1),nr(:,2), rp(:,j,k),nrw(:,j,k))
    enddo
    enddo
 ENDIF

 IF (pcos2 ==1) THEN
    print*,' a factor 1/cos(th)^2 is added to the term d^2/dph^2!----'
 ENDIF

  ii=1      ! row number of R_sparse matrix
  IA(1)=1  ! index in R%aa for 1st non-0 element in 1st row of R_sparse

  DO k=0, sr-1
  DO j=1, st-2     !;excluding voxels with th=+/-pi/2 (singular points)
  DO i=0, sp-1
                   !;Reg for longitude: d^2/dp^2 = (N(i-1)-2*N(i)+N(i+1))
    if (i == 0) then
       Npi%n=3
       iw1=sp+ j*sp+k*sp*st
       iw2= 1+ j*sp+k*sp*st
       iw3= 2+ j*sp+k*sp*st
       Npi%ja(:)=(/iw2,iw3,iw1/)
       Npi%aa(:)=(/-2, 1, 1/)
    else if (i == sp-1) then
       Npi%n=3
       iw1=(i+0)+j*sp+k*sp*st
       iw2=(i+1)+j*sp+k*sp*st
       iw3= 1   +j*sp+k*sp*st
       Npi%ja(:)=(/iw3,iw1,iw2/)
       Npi%aa(:)=(/1, 1, -2/)
    else
       Npi%n=3
       iw1=(i+0)+j*sp+k*sp*st
       iw2=(i+1)+j*sp+k*sp*st
       iw3=(i+2)+j*sp+k*sp*st
       Npi%ja(:)=(/iw1,iw2,iw3/)
       Npi%aa(:)=(/1, -2, 1/)
    endif
    if (pcos2==1) then
        Npi%aa=Npi%aa/cos(thp(i+1,j+1,k+1))**2
    endif

        !;Reg for latitude: for d^2/dt^2 = (N(j-1)-2*N(j)+N(j+1))
     Nti%n=3
     iw1=(i+1)+(j-1)*sp+k*sp*st
     iw2=(i+1)+j*sp+k*sp*st
     iw3=(i+1)+(j+1)*sp+k*sp*st
     Nti%ja(:)=(/iw1,iw2,iw3/)
     Nti%aa(:)=(/1, -2, 1/)

!        ;Reg for radial: d^2/dr^2 = (N(k-1)-2*N(k)+N(k+1))
    if (k == 0) then 
       Nri%n=2      !;right diff: (N(k+1)-N(k))/dr
       iw1=(i+1)+j*sp+k*sp*st
       iw2=(i+1)+j*sp+(k+1)*sp*st
       Nri%ja(:)=(/iw1,iw2/)
       Nri%aa(:)=(/-1, 1/)
    else if (k == sr-1) then
       Nri%n=2      !;left diff: (N(k)-N(k-1))/dr
       iw1=(i+1)+j*sp+(k-1)*sp*st
       iw2=(i+1)+j*sp+k*sp*st
       Nri%ja(:)=(/iw1,iw2/)
       Nri%aa(:)=(/-1, 1/)
    else
       Nri%n=3
       iw1=(i+1)+j*sp+(k-1)*sp*st
       iw2=(i+1)+j*sp+k*sp*st
       iw3=(i+1)+j*sp+(k+1)*sp*st
       Nri%ja(:)=(/iw1,iw2,iw3/)
       Nri%aa(:)=(/1, -2, 1/)
    endif


   if (nr(1,1)>0) then
      Npi%aa=Npi%aa/nrw(i+1,j+1,k+1)
      Nti%aa=Nti%aa/nrw(i+1,j+1,k+1)
      Nri%aa=Nri%aa/nrw(i+1,j+1,k+1)
   endif


      nja=Npi%n
      call Allocate_Rys(mRys(ii),nja)
      mRys(ii)%ja = Npi%ja(1:nja)
      mRys(ii)%aa = Npi%aa(1:nja)
      IA(ii+1) = IA(ii)+nja
      ii=ii+1

      nja=Nti%n
      call Allocate_Rys(mRys(ii),nja)
      mRys(ii)%ja = Nti%ja(1:nja)
      mRys(ii)%aa = Nti%aa(1:nja)
      IA(ii+1) = IA(ii)+nja
      ii=ii+1

      nja=Nri%n
      call Allocate_Rys(mRys(ii),nja)
      mRys(ii)%ja = Nri%ja(1:nja)
      mRys(ii)%aa = Nri%aa(1:nja)
      IA(ii+1) = IA(ii)+nja
      ii=ii+1
  ENDDO
  ENDDO
  ENDDO

 mm=ii-1
 nnz = IA(mm+1)-1

  print*,'ii=',ii,' mm=',mm, '  Should be: ii=mm+1'
  print*, 'nnz=',nnz

 call Allocate_SP(R, mm, nn, nnz)
 R%ifst = IA
 DO i=1,mm
    R%icol(IA(i):IA(i+1)-1) = mRys(i)%ja
    R%A(IA(i):IA(i+1)-1) = mRys(i)%aa
    call DeAllocate_Rys(mRys(i))
 ENDDO
 call DeAllocate_Rys(Npi)
 call DeAllocate_Rys(Nti)
RETURN
END FUNCTION sph_Reg2nd2rtp_SP


FUNCTION sph_Reg2nd2smtp_SP (php,thp,rp, nr,pcos2) RESULT(R)
!;+ purpose: Second-order smoothing using central differetiate of 2nd-order accuracy
!;     with the form, |RX|^2=Sum((d^2/dth^2)^2 +(d^2/dph^2)^2)
!;   same as in Frazin etal, 2007, ApJL, 671,L201
!; Input:
!;   php,thp,rp[sz1_ph,sz_th,sz_r]-- 3D spheric grid
!; keyword:
!;   nr=nr
!;     nr[0],nr[1] gives the dependence of normalized N(r) on radial distance nr[0]
!;     If set, add a radial weight w=1/N(r) to the regularization term Rw^T#Rw
!;   pcos2==1
!;      If set, add a smooth factor 1/cos^2(phi) to d^2/dph^2)^2
!; return: R -mm x nn sparse matrix for regularization, where
!;          nn=sz1_ph*sz_th*sz_r and mm=3*nn
!; add weight=1/nr for each grid point remove r-dependence effect of regularization
USE InternPro
IMPLICIT NONE
REAL (kind=double), INTENT(IN), DIMENSION(:,:,:) :: php,thp,rp
REAL (kind=double), INTENT(IN), DIMENSION(:,:) :: nr
INTEGER, INTENT(IN) :: pcos2
TYPE (SparseMat) :: R
!REAL (kind=double) :: dph,dth,dr
INTEGER :: sp,st,sr,mm,nn,nnz, ii, IA(3*size(php)+1), i,j,k, iw1,iw2,iw3,nja
TYPE (dataRys) :: Npi, Nti, mRys(3*size(php)) 
REAL (kind=double),DIMENSION(:,:,:), ALLOCATABLE :: nrw 

 allocate(nrw, mold=rp)  ! nrw is allocated with the bound of rp

!  dr=(grid_Rout-grid_Rin)/dble(sz_r-1)  !grid scale in r-direction
! dth=PI/dble(sz_th-1)                   !grid scale in th-direction
! dph=2*PI/dble(sz_ph-1)                 !grid scale in th-direction
 print*,'Regularization term is chosen as'
 print*,'d^2/dph^2 =(N(i-1)-2*N(i)+N(i+1)), d^2/dth^2 =(N(j-1)-2*N(j)+N(j+1))'

 sp=SIZE(php, DIM=1)
 st=SIZE(php, DIM=2)
 sr=SIZE(php, DIM=3)
 nn =sp*st*sr       !;number of unknown elements in solution vector, X
! mm =2*nn           ! initial row number of regularization matrix, R
print*,'sp,sr,sr=',sp,st,sr

  call Allocate_Rys(Npi, 3)
  call Allocate_Rys(Nti, 3)

 IF (nr(1,1)>0) THEN
    print*,'--- Regularization with radial weighting w=1/Nr is applied!'
    do k=1,sr
    do j=1,st
       call interp(nr(:,1),nr(:,2), rp(:,j,k),nrw(:,j,k)) 
    enddo
    enddo
 ENDIF

 IF (pcos2 ==1) THEN
    print*,' a factor 1/cos(th)^2 is added to the term d^2/dph^2!----'
 ENDIF 

  ii=1      ! row number of R_sparse matrix
  IA(1)=1  ! index in R%aa for 1st non-0 element in 1st row of R_sparse

  DO k=0, sr-1 
  DO j=1, st-2     !;excluding voxels with th=+/-pi/2 (singular points)
  DO i=0, sp-1 
                   !;Reg for longitude: d^2/dp^2 = (N(i-1)-2*N(i)+N(i+1))
    if (i == 0) then
       Npi%n=3
       iw1=sp+ j*sp+k*sp*st
       iw2= 1+ j*sp+k*sp*st
       iw3= 2+ j*sp+k*sp*st
       Npi%ja(:)=(/iw2,iw3,iw1/)
       Npi%aa(:)=(/-2, 1, 1/)
    else if (i == sp-1) then
       Npi%n=3
       iw1=(i+0)+j*sp+k*sp*st
       iw2=(i+1)+j*sp+k*sp*st
       iw3= 1   +j*sp+k*sp*st
       Npi%ja(:)=(/iw3,iw1,iw2/)
       Npi%aa(:)=(/1, 1, -2/)
    else
       Npi%n=3
       iw1=(i+0)+j*sp+k*sp*st
       iw2=(i+1)+j*sp+k*sp*st
       iw3=(i+2)+j*sp+k*sp*st
       Npi%ja(:)=(/iw1,iw2,iw3/)
       Npi%aa(:)=(/1, -2, 1/)
    endif
    if (pcos2==1) then
        Npi%aa=Npi%aa/cos(thp(i+1,j+1,k+1))**2
    endif

        !;Reg for latitude: for d^2/dt^2 = (N(j-1)-2*N(j)+N(j+1))
     Nti%n=3
     iw1=(i+1)+(j-1)*sp+k*sp*st
     iw2=(i+1)+j*sp+k*sp*st
     iw3=(i+1)+(j+1)*sp+k*sp*st
     Nti%ja(:)=(/iw1,iw2,iw3/)
     Nti%aa(:)=(/1, -2, 1/)

   if (nr(1,1)>0) then
      Npi%aa=Npi%aa/nrw(i+1,j+1,k+1)
      Nti%aa=Nti%aa/nrw(i+1,j+1,k+1)
   endif

      nja=Npi%n
      call Allocate_Rys(mRys(ii),nja)
      mRys(ii)%ja = Npi%ja(1:nja)
      mRys(ii)%aa = Npi%aa(1:nja)
      IA(ii+1) = IA(ii)+nja
      ii=ii+1

      nja=Nti%n
      call Allocate_Rys(mRys(ii),nja)
      mRys(ii)%ja = Nti%ja(1:nja)
      mRys(ii)%aa = Nti%aa(1:nja)
      IA(ii+1) = IA(ii)+nja
      ii=ii+1
  ENDDO
  ENDDO
  ENDDO

 mm=ii-1
 nnz = IA(mm+1)-1

   print*,'ii=',ii,' mm=',mm, '  Should be: ii=mm+1'
   print*, 'nnz=',nnz

 call Allocate_SP(R, mm, nn, nnz)
 R%ifst = IA
 DO i=1,mm
    R%icol(IA(i):IA(i+1)-1) = mRys(i)%ja
    R%A(IA(i):IA(i+1)-1) = mRys(i)%aa
    call DeAllocate_Rys(mRys(i))
 ENDDO
 call DeAllocate_Rys(Npi)
 call DeAllocate_Rys(Nti)
RETURN
END FUNCTION sph_Reg2nd2smtp_SP

END MODULE TomoPro


!------------------------------------------------
!  Main program 
! Final version for Spheric grid code with OPenMP, modified on 5/5/2025 by Tongjiang Wang
!  the solution of electron density in units of cm^-3
program tomo_sph_omp
 use InternPro
 use TomoPro
 use OMP_LIB
 implicit none

integer, parameter :: nimg=60, isz=128  !for grid: sz_r=51, sz_th=181, sz_ph=361
! isz: smaple (rebinned) pixel number
! sample resolution not more than grid resolution, keeping solu.  underdetermined
integer :: nt     !number of used images
real (kind=8), parameter :: Lse=215.0_double
real (kind=8), parameter ::  epswt=1.e-4 !;threshold for none-zero elements in matrix A
real (kind=double) :: dr,dth,dph, ph0(sz_ph),th0(sz_th),r0(sz_r), ph1(sz_ph-1)
integer :: i,j,k, sz1_ph, ii,jj,kk,nn,mm   
!note  storage_size(integer)=32, in range [-2^31, 2^31]
real (kind=double), dimension(sz_ph-1,sz_th,sz_r) :: php,thp,rp, weightn
integer, dimension(:), allocatable:: px1d
integer :: im, ik, nt0, npw, npw0, nxr, ywt=1
real (kind=double),  dimension(:,:), allocatable :: ix,iy, nr
type (dataRys) :: arrRys(nimg*isz*isz), Rys(isz*isz)
integer, dimension(:), allocatable :: IAt, IA, Rys_nja
integer :: allocatestatus, Openstatus
type (mapstr) :: map
type (SparseMat) :: Amat, Rmat
real (kind=double), dimension(:), allocatable :: ixw,iyw,izw,irw,irho,ipb,&
                        cx,cy,cz,lon,lat,dse, yd0,ydw, xr,ybk
real (kind=double), dimension(:,:), allocatable :: xcrd,ycrd,irp, cobs_xy
real (kind=double) :: Lobs(nimg), cobs(2,nimg), cobs0(2)
real, dimension(:), allocatable :: ydat0,ydata
real :: ydij
real (kind=double), dimension(:), allocatable:: re,the,phe
integer :: cross, nw, nja, istt,iend
integer, dimension(:,:), allocatable:: ccell
real (kind=double), dimension(3) :: ni
real (kind=double), dimension(2) :: rt,tht,pht
real (kind=double) :: weights (sz_ph,sz_th,sz_r) 
real (kind=double), dimension(sz_ph*sz_th*sz_r) :: weightn0,php0,thp0,rp0
integer :: wja(sz_ph*sz_th*sz_r), nnz
real (kind=double), dimension(:), allocatable:: weightnz,phw,thw,rw,pxw,pyw,pzw,&
              rp2,sp,Rho,theta,kr
integer :: write_Abuff, read_Abuff, write_Rbuff, read_Rbuff, wt_nr, pcos2
real :: mu,tol
real, dimension(:), allocatable :: xsolu
real, dimension(sz_ph, sz_th, sz_r) :: nsolu
integer :: sav_res, sav_iter, resume

real :: tstart,tfinish, e,t(2), e1,e2,t1(2),t2(2)
real :: wstart,wfinish,wrys

INTEGER :: start_clock, end_clock, count_rate
REAL :: elapsed_time

real (kind=8) :: amax,amin,amin0=1e10,amax0=1e-10, n_unit=1.0e7
integer, allocatable, dimension(:) :: wmin,wmax
integer :: nmax,nmin

character (len=:),allocatable :: Abuffname, Rbuffname, gridname
character (len=:),allocatable :: sdir, cmu
character (len=200) :: filename,ybkname,nrbkname,obsname, dir
character (len=200) :: namesolu, pbnames(nimg)
character (len=3) :: snum, sn1,sn2,sn3
integer :: imgnum
logical :: exist

integer :: CHUNK=100, tid

   CALL SYSTEM_CLOCK(count_rate)   ! Get the clock rate
   CALL SYSTEM_CLOCK(start_clock)  ! Get the start time

!$ wstart = OMP_get_wtime()
call cpu_time(tstart)
e=dtime(t1)
e1=e 
print*,'elapsed:', e, ', user:', t1(1), ', sys:', t1(2)

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
 enddo

! set up parameters for running
!wt_nr=1 ! for case with wt_nr=1
!mu=1.0  
!reg0=0  ! for 2nd-order smoothing
!mu=0.001 !for zero-th order smooth
!tol=0.001

!write_Abuff=1
!read_Abuff =0
!write_Rbuff=0
!read_Rbuff =0

! resume=0
! sav_res=1
! sav_iter=1
pcos2=0

write(sn1,'(I3)') sz_ph
write(sn2,'(I3)') sz_th
write(sn3,'(I3)') sz_r
 gridname='p'// trim(adjustL(sn1))// 't' // trim(adjustL(sn2)) &
    //'r'// trim(adjustL(sn3))

 if (wt_nr == 1) then
     Rbuffname='Buff_Rsp_reg2_' // gridname // '_nwt.dat'
     namesolu='nsolu_'// gridname //'_reg2_nwt.dat'
 else
     Rbuffname='Buff_Rsp_reg2_' // gridname // '.dat'
     namesolu='nsolu_'// gridname //'_reg2.dat'
 endif

 Rbuffname=sdir//trim(adjustL(Rbuffname))
 namesolu =sdir//trim(adjustL(namesolu))
 Abuffname=sdir//'Buff_Asp_' // gridname //'.dat'

 print*,'A, R Buff, Solu-name=',Abuffname, ' ', Rbuffname, ' ', namesolu


  dr=(grid_Rout-grid_Rin)/dble(sz_r-1)  !grid scale in r-direction
 dth=PI/dble(sz_th-1)                   !grid scale in th-direction
 dph=2*PI/dble(sz_ph-1)                 !grid scale in th-direction

 ph0=(/(0.d0 + dph*(i-1),i=1,sz_ph)/)
 th0=(/(-PI/2 + dth*(i-1),i=1,sz_th)/)
 r0=(/(grid_Rin + dr*(i-1),i=1,sz_r)/)

 ph1=ph0(1:sz_ph-1)           !correct ph in [0, 2*pi-dph] to
 sz1_ph=sz_ph-1               !remove overlaping of 0 with 2pi in solution
!nn=sz1_ph*sz_th*sz_r         !number of unknown elements in solution vector, X

!$OMP PARALLEL 
!$OMP master
  !$ print*,'Total used thread num=',omp_get_num_threads()
!$OMP end master
!$OMP WORKSHARE
forall(i=1:sz_th, j=1:sz_r)
  php(:,i,j)=ph1
end forall

forall(i=1:sz1_ph, j=1:sz_r)
  thp(i,:,j)=th0
end forall

forall(i=1:sz1_ph, j=1:sz_th)
  rp(i,j,:) =r0
end forall
!$OMP END WORKSHARE
!$OMP END PARALLEL 


 IF (read_Abuff ==1) GO TO 400

!define image coord. for observer
print*,'... Reading pB-images into map-struct----'

  mm=nt*isz*isz   !set initial value for column number of A_sparse
  ALLOCATE(px1d(isz),ix(isz,isz),iy(isz,isz),xcrd(isz,isz),ycrd(isz,isz),irp(isz,isz),&
          ixw(isz*isz),iyw(isz*isz),izw(isz*isz),irw(isz*isz),&
          yd0(isz*isz), ydat0(nt*isz*isz), IAt(nt*isz*isz))
  px1d=(/(i,i=1,isz)/) !index for image pixel coordinate
! create pix coord. for image
 DO k=1, isz
    ix(:,k)=real(px1d)
    iy(k,:)=real(px1d)
 ENDDO
! Reading (distance, lon, lat) for Observer,
!    where units are Rsun and degree
 OPEN (UNIT=10, FILE=obsname,ACTION='READ',&
            IOSTAT = Openstatus )

 IF (OpenStatus /= 0) THEN
      STOP "error opening file [cobs_coord_?.txt]"
 ELSE
    print*, '---file ['//obsname//'] is read ....'
 ENDIF
    

 READ(10,*) nt0
 print*,'nt0=',nt0
 DO k=1,nt0
   READ(10,*) Lobs(k), cobs(:,k)
  !  print*,Lobs(k), cobs(:,k)
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
 img_loop: DO im=1, imgnum  !main loop
    write(snum,'(I3)') im
!   filename=trim(imgname) // trim(adjustL(snum)) // '.dat'
    filename=pbnames(im)
    print*,'--read img no=',im,'  file=',filename
    map=read_mapstr(filename)
    call getmap_coord(map,ix,iy, xcrd,ycrd)
!$OMP PARALLEL WORKSHARE  
    xcrd=xcrd/map%rsun !convert img coord from arcsec into Rsun
    ycrd=ycrd/map%rsun
    irp=sqrt(xcrd*xcrd+ycrd*ycrd)
!$OMP END PARALLEL WORKSHARE
    jj=1
    DO 200 j=1,isz
    DO 200 i=1,isz
       IF (irp(i,j)>=grid_Rin .and. irp(i,j)<=grid_Rout) THEN
          ydij=map%data(i,j)
          if (ydij >= 0) then ! this will exclude data=NaN
           ixw(jj)=xcrd(i,j)
           iyw(jj)=ycrd(i,j)
           irw(jj)=irp(i,j)
           yd0(jj)=ydij
           jj=jj+1
          else if (ydij <0 .or. ydij /=ydij) then
              print*,'data=',ydij, ' at (i=',i,', j=',j,') <0 or is NaN, excluded!---'
          endif
       ENDIF
    200 CONTINUE
    npw=jj-1
    izw=Lobs(im)
    ALLOCATE(cx(npw),cy(npw),cz(npw),irho(npw),ipb(npw),ydw(npw),&
             lon(npw),lat(npw),dse(npw),cobs_xy(2,npw),Rys_nja(npw))
    irho=irw(1:npw)  !distanc from Sun-center to a pixel in img-plane
    ydw =yd0(1:npw)  !pB-data for the pixel within given domain
    call img2car_coord(cobs(:,im),-ixw(1:npw),-iyw(1:npw),izw(1:npw), cx,cy,cz)
!Lobs is the distance from Sun center to Observer in unit of Rsun
!   ie., Lobs=Lcor1a or Lcor1b
    call car2sp(cx,cy,cz, lon,lat,dse)  ![lon,lat] in radian
!     dse id the distance from observer to the pixel in img-plane
!  cobs(1:2,:) give the LOS directions of all pixels in images meeting Rin<r<Rout 
!  from observer
    cobs_xy(1,:)=mod(lon*180.d0/PI+360.d0,360.d0)
    cobs_xy(2,:)=lat*180.d0/PI

!  if it is set to normalize the obsrved and modeled data
!   by rdial distribution of pB background obtained obs
    IF (ywt ==1) THEN
       call interp(xr,ybk, irho,ipb)
       ydw=ydw/ipb  !normalize the data by radial pB background
    ENDIF
   ydat0(npw0:npw0+npw-1)=ydw
   npw0= npw0+ npw

!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(php,thp,rp,cobs_xy,ixw,iyw,npw,Rys,Rys_nja,sz1_ph,ipb,ywt,n_unit) 
!$OMP DO SCHEDULE(DYNAMIC,CHUNK) 
      pix_loop: DO ik=1,npw  ! second loop
          cobs0=cobs_xy(:,ik)
          call sph_cross_points (cobs0, ixw(ik),iyw(ik),0.d0,&
               phe,the,re,ccell, cross,ni,pht,tht,rt)
      cross_if: IF (cross ==1) THEN 

        weights=Rsun*sph_ray_weights(phe,the,re,ccell,ni)  !;in units of cm.

!  ;Contribution of the grid to line integration for a ray
!  ;check its validation: total weights  for a ray should equal to total length
!  ;      between two crossing points within the spheric grid

!       print*,'## Total weights for all voxels in Rsun=',SUM(weights)/Rsun

! ;Important: convert grid from (sz_ph,sz_th,sz_r) to (sz1_ph,sz_th,sz_r)
! ;  to avoid the repetitive voxels php,thp,rp[2*pi, *,*] in solution by
! ; adding the weights for voxels with ph=2*pi to ph=0 as they are same
       weightn=weights(1:sz1_ph, :, :)
       weightn(1,:,:)=weights(1,:,:) + weights(sz_ph,:,:)

! Add the contribution of Thomson scatter function for single electron
             ! find non-zero pixeks in 3d-weights
      nw=0
      DO 300 k=1, sz_r
      DO 300 j=1, sz_th
      DO 300 i=1, sz1_ph
         if (weightn(i,j,k) > epswt) then
            nw=nw + 1
            weightn0(nw) = weightn(i,j,k)
            php0(nw)  = php(i,j,k)
            thp0(nw)  = thp(i,j,k)
            rp0(nw)   = rp(i,j,k)
            wja(nw)   = i + (j-1)*sz1_ph + (k-1)*sz1_ph*sz_th  !1D index must begin with 1
         endif
      300 CONTINUE
      nja=nw  ! number of non-zero elements in weightn
!  ;distance from Sun center to point P, rp=|OP(px1,py1,pz1)| =sqrt(rp2)
      ALLOCATE(weightnz(nja),phw(nja),thw(nja),rw(nja),pxw(nja),pyw(nja),pzw(nja),&
               rp2(nja),sp(nja),Rho(nja),theta(nja),kr(nja))
       weightnz =weightn0(1:nja)

       phw =php0(1:nja)
       thw =thp0(1:nja)
       rw  =rp0(1:nja)
       call sp2car(phw,thw,rw, pxw,pyw,pzw)

!       Calculate the radial-dependent weight function wi=1/Yi(r)_bk
!  note: apply w=1/ypb here to weightnz to avoid values of "weightnz"
!  too small leading overfloating a range of 10^+/-38
       if (ywt ==1) weightnz =weightnz/ipb(ik)

! sp is the distance from P to POS
!       where POS (plane-of-sky) is a plane that crosses the Sun center and
!        perpendicular to LOS, s=|vec_OP cdot vec_ni|=|rp*cos(alp)|, where alp is
!        the angle between OP and LOS, and vec_ni is the unit vector of LOS.
        rp2 =pxw*pxw+pyw*pyw+pzw*pzw
        sp=abs(pxw*ni(1)+pyw*ni(2)+pzw*ni(3))
        Rho=sqrt(rp2-sp*sp)  ! rho is the distance from the Sun center to LOS
        theta=atan(sp/Rho)     ! tht is the angle between vector OP and POS
!        kr =thomsonscatter(Rho,theta,-1.,.true.) !tht in unit of radian
        weightnz = weightnz*(thomsonscatter(Rho,theta,-1.,.true.)*n_unit)    !add effect of Thmoson scatter, I=Sum(Kr*Nr*dl)
                                            !set density in unit of 10.^7 cm^-3

         !construct Matrix A in CRS format (CRS: Compressed Row Sparse)
             call Allocate_Rys(Rys(ik), nja)
             Rys(ik)%ja =wja(1:nja)
             Rys(ik)%aa =weightnz
             Rys_nja(ik) = nja
       ELSE
             STOP "--- cross==0 returned from sub: get_cross_points"// &
                "or the ray crossed grid_Rout, special treatment required!"
       ENDIF cross_if
       DEALLOCATE(weightnz,phw,thw,rw,pxw,pyw,pzw,rp2,sp,Rho,theta,kr)
    ENDDO  pix_loop
!$OMP END DO
!$OMP END PARALLEL
     arrRys(ii:ii+npw-1)=Rys(1:npw)
     DO k=1, npw 
        IAt(ii+k)=IAt(ii+k-1)+Rys_nja(k)
        call DeAllocate_Rys (Rys(k))
     ENDDO
     ii=ii+npw   
    DEALLOCATE(cx,cy,cz,irho,ipb,ydw,lon,lat,dse,cobs_xy,Rys_nja)
    
 ENDDO  img_loop
  mm=ii-1  ! total number of sampled rays or data points
print*,'mm=',mm, 'should be = npw0-1=',npw0-1
  nn=sz1_ph*sz_th*sz_r         !total number of grid points in Domain of solution
                               !or number of unknown elements in solution vector, X
  ALLOCATE(ydata(mm), IA(mm+1))
  ydata=ydat0(1:mm)
  IA=IAt(1:mm+1)
  nnz =IA(mm+1)-1   !total number of non-zero elements in A_sparse matrix
  DEALLOCATE(px1d,ix,iy,xcrd,ycrd,irp,ixw,iyw,izw,irw,yd0,ydat0,IAt)

e=dtime(t)
print*,'elapsed:', e, ' (s)  for creating mat A_sparse----'
!$ wrys=OMP_get_wtime()
!$ print*,'per-thread time=',wrys-wstart



! save Data for Coeff. Matrix A in a buff file on hard driver

 IF (write_Abuff ==1) THEN
    OPEN(1, FILE=Abuffname, FORM="unformatted",ACTION="write", STATUS="replace")
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
   print*,'Buff file:[' // Abuffname // '] is written! ===>'
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
   print*,'Buff file:'//Abuffname// '] is read!==='
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
 print*,'min, max of ydata=',minval(ydata),maxval(ydata)


 IF (wt_nr ==1) THEN
     OPEN(20, FILE=ybkname,ACTION='read')
     READ(20,*) nxr
     ALLOCATE(nr(nxr,2))
     READ(20,*) nr(:,1), nr(:,2)
     CLOSE(20)
    nr(:,2)=nr(:,2)/maxval(nr(:,2))
    print*,'radial weightpb(r): '//ybkname// ' is read===>'
 ELSE
    ALLOCATE(nr(1,2))
    nr=0.d0
   print*,'No radial weight wt=1/Nr is applied---!'
 ENDIF


! Rbuffname='Buff_Rsp_' // trim(adjustL(filename)) // '_reg2nd2smtp.dat'
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
    print*,'---Creating the Regularization Matrix.....'
    Rmat= sph_Reg2nd2rtp_SP(php,thp,rp, nr,pcos2)  !smooth in r, th, and ph-directions
!   Rmat= sph_Reg2nd2smtp_SP(php,thp,rp, nr,pcos2) !smooth only in ph and th-direction
!  save Data from sph_Reg2nd2smtp_SP() for regularization R_sparse in a buff file on hard driver
   IF (write_Rbuff ==1) THEN
      OPEN(1, FILE=Rbuffname, FORM="unformatted",ACTION="write",STATUS="replace")
      WRITE(1) Rmat%m, Rmat%n, Rmat%asiz
      WRITE(1) Rmat%ifst
      WRITE(1) Rmat%icol
      WRITE(1) Rmat%A
      CLOSE(1)
     print*,'Buff file: ['//Rbuffname//'] is written! ===>'
   ENDIF
 ENDIF

   print*,' mat R_sparse is created: mm,nn,nnz=',Rmat%m, Rmat%n, Rmat%asiz
   print*,'min max R=',minval(Rmat%A),maxval(Rmat%A)


 print *,'--------------Test for Conjugate Gradient to minimize'
 print *,' f=min(|Ax-y|^2+mu*|R*x|^2 '

 ALLOCATE(xsolu(Amat%n))
 print*,'.... with mu=',mu, '  tol=',tol

 xsolu= ConjugateGrad_SP (Amat,ydata,mu,Rmat, tol,sav_res,sav_iter,resume)

print*,'solu min max=',minval(xsolu),maxval(xsolu)
 
 ! convert 1D solution into 3D
 nsolu(1:sz1_ph,:,:) = reshape(xsolu, shape(nsolu(1:sz1_ph,:,:)))
 nsolu(sz_ph,:,:)    = nsolu(1,:,:)
!;copy values of voxels with phi=0 to phi=2*pi
! ; return solution with ph=[0,2*pi], th=[-pi/2, pi/2], r=[Rin. Rout]

 INQUIRE(FILE=namesolu, exist=exist)
  if (.not. exist) then
         open(unit=10,file=namesolu, form='unformatted',status="new",action="write",iostat=openstatus)
  else
         print*,'The file ['//namesolu//'] already exists!'
         namesolu='New_'//namesolu
         print*,'The solution will be saved in ['//namesolu//']--------!'
         open(unit=10,file=namesolu, form='unformatted',status="new",action="write",iostat=openstatus)
  endif

  if (Openstatus >0) stop "----Error, File not opened properly. &
                         &it  maybe already exist---"
  
  nsolu= nsolu*n_unit   ! convert unit of electron density into cm^-3             
write(10) sz_ph, sz_th, sz_r
write(10) ph0, th0, r0  !type of double
write(10) nsolu         ! typr of real
close(10)
print*,'Solution file ['//namesolu//'] is created!'

print*,'elapsed:', e, ' (s)  for creating mat A_sparse----'

 e=dtime(t2)
 print*,'end elapsed:', e, ' (s)  for obtaining solution using conjugate gradient----'

call cpu_time(tfinish)
!$ wfinish=OMP_get_wtime()

print*,'computing time=',tfinish-tstart
!$ print*,'per-thread time=',wfinish-wstart


   CALL SYSTEM_CLOCK(end_clock)    ! Get the end time
   elapsed_time = REAL(end_clock - start_clock)  / REAL(count_rate)
   PRINT *, "Elapsed time (seconds): ", elapsed_time

end program tomo_sph_omp

