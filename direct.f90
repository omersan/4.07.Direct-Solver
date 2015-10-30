!-----------------------------------------------------------------------------!
! MAE 5093 DOCUMENTATION | Engineering Numerical Analysis
!-----------------------------------------------------------------------------!
! >>> Direct method to solve boundary value problems 
!     Uses 2nd-order Finite Difference
!-----------------------------------------------------------------------------!
! References: 
! * Fundamentals of Engineering Numerical Analysis by P. Moin (2012) -- Ch.1
! * Numerical Recipes: The Art of Scientific Computing, 2nd Edition (1992) 
!-----------------------------------------------------------------------------!
! Written by Omer San
!            CFDLab, Oklahoma State University, cfdlab.osu@gmail.com
!            www.cfdlab.org
! 
! Last updated: Sep. 29, 2015
!-----------------------------------------------------------------------------!


program shooting
implicit none
real*8 ::h,x0,xL,xe,ye,ya,yb
integer::j,n,np
real*8,allocatable ::y(:)
real*8,allocatable ::a(:),b(:),c(:),r(:),x(:)

!Solve y'' - 4*y = 0, for y(x), x in [0,1]
!y(0) = 0.0d0
!y(1) = 100.0d0 

n=10
allocate(y(0:n))

x0=0.0d0
xL=1.0d0
h=(xL-x0)/dfloat(n)

!Boundary conditions
ya = 0.0d0
yb = 100.0d0

y(0) = ya
y(n) = yb

!Build coefficient matrix:
allocate(a(1:n-1),b(1:n-1),c(1:n-1),r(1:n-1),x(1:n-1))


do j=1,n-1
a(j) = 1.0d0/(h*h)
b(j) =-2.0d0/(h*h) - 4.0d0
c(j) = 1.0d0/(h*h)
r(j) = 0.0d0
end do
!apply boundary conditions
r(1)   = r(1) - a(1)*y(0)     !b.c.
r(n-1) = r(n-1) - c(n-1)*y(n) !b.c.
call tdma(a,b,c,r,x,1,n-1)

!assign solutions to y
do j=1,n-1
y(j)=x(j)
end do

!Make final computation one more time to plot results 
open(12, file="numerical.plt")
write(12,*)'variables ="x","y"'
do j=0,n
write(12,*) dfloat(j)*h,y(j)
end do
close(12)   

    

   
! Writing exact solution using np points
np = 4000 !use 4000 points to plot curve
open(12, file="exact_curve.plt")
write(12,*)'variables ="x","y"'
	do j=0,np
		xe = dfloat(j)*(xL-x0)/(dfloat(np))
        ye = 100.0d0/(dexp(2.0d0)-dexp(-2.0d0))*(dexp(2.0d0*xe)-dexp(-2.0d0*xe))
		write(12,*) xe,ye
	end do
close(12)

end


!------------------------------------------------------------------!
!Tridiagonal matrix algorithm (TDMA)
!Thomas algorithm
!solution tridiagonal systems
!a: lower diagonal
!b: main diagonal
!c: upper diagonal
!r: source vector
!x: solution vector
!   for indices s(start) to e(end)
!   i: s,s+1,s+2, ....,i,....,e 
!
!Note: a(s) and c(e) are dummy coefficients, not used.
!------------------------------------------------------------------!

subroutine tdma(a,b,c,r,x,s,e)
implicit none
integer s,e,i
real*8, dimension(s:e) ::a,b,c,r,x    

! forward elimination phase
do i=s+1,e
b(i) = b(i) - a(i)/b(i-1)*c(i-1)
r(i) = r(i) - a(i)/b(i-1)*r(i-1)
end do
! backward substitution phase 
x(e) = r(e)/b(e)
do i=e-1,s,-1
x(i) = (r(i)-c(i)*x(i+1))/b(i)
end do

return
end







