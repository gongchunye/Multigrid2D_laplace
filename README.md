# 5.06.Multigrid2D
V-cycle Multigrid method for 2D Poisson Equation，代码写得很清楚，是非常好的学习材料


![residual](https://cloud.githubusercontent.com/assets/15114859/10856500/27fbcbd0-7f16-11e5-96af-89f90d1cfa5e.png)
![residual_log](https://cloud.githubusercontent.com/assets/15114859/10856503/299b622a-7f16-11e5-989e-769513abe097.png)
![field](https://cloud.githubusercontent.com/assets/15114859/10856504/2aa99ede-7f16-11e5-840f-220a99dc9485.png)

1)如何编译运行
rm -f a.out

gfortran -O3 poisson2d_MG.f90

./a.out

2)
!----------------------!
!Solver:
!----------------------!
isolver = 1
if (isolver.eq.0) then !Gauss-Seidel scheme
	call GS(nx,ny,dx,dy,f,u,tol)
else
	call MG5(nx,ny,dx,dy,f,u,tol)
end if

这个地方修改isolver = 0 可以转化成Gauss-Seidel迭代

3）

ireo = 2

if (ireo.eq.1) then !simply injection

	do i=1,nxh-1
	do j=1,nyh-1
		f(i,j) = r(2*i,2*j) 							  	
	end do
	end do

else if (ireo.eq.2) then !half-weight

	do i=1,nxh-1
	do j=1,nyh-1
		f(i,j) = 1.0d0/8.0d0*( 4.0d0*r(2*i,2*j) &
			 + 1.0d0*(r(2*i+1,2*j)+r(2*i-1,2*j)+r(2*i,2*j+1)+r(2*i,2*j-1)) )							  	
	end do
	end do

else !full-weight (trapezoidal)

do i=1,nxh-1
do j=1,nyh-1
	f(i,j) = 1.0d0/16.0d0*( 4.0d0*r(2*i,2*j) &
	     + 2.0d0*(r(2*i+1,2*j)+r(2*i-1,2*j)+r(2*i,2*j+1)+r(2*i,2*j-1)) &
	     + 1.0d0*(r(2*i+1,2*j+1)+r(2*i-1,2*j-1)+r(2*i-1,2*j+1)+r(2*i+1,2*j-1)))							  	
end do
end do

修改ireo可以换不现的f2c的插值方式，但是效果都差不多。


