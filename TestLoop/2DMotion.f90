	program bragg

	IMPLICIT NONE
	integer :: i
	common/par/ m,g,alpha
	real m,g,alpha
	common/tole/ tol
	real tol
	common/start/ ystart,ii
	real,dimension(11) :: ystart
	integer ii
	real Pi

	OPEN(UNIT=11,FILE='etest.dat')        
		 Pi=acos(-1.)
         tol=0.01E0
	call parameters
	do ii=1,11
		ystart(ii)=2.*Pi*(ii-1)/10.
		call result
	end do

	close(11)
	STOP
	END

	 

	
	subroutine parameters
	implicit none
	common/par/ m,g,alpha
	real m,g,alpha
	m=1E0
	g=10E0
	alpha=0.1E0
	return
	end



        subroutine derivs(x,y,dydx)
        implicit none
        real*8 x,y(4),dydx(4)
        
        dydx(1)=y(2)
        dydx(2)=0
        dydx(3)=y(4)
        dydx(4)=1.*sin(y(3))

        end


	 
	
	subroutine result
	implicit none

	common/tole/ tol
	real tol

	real endoftim
	integer mxparm,neq
	parameter (mxparm=50,neq=4)
	integer ido,tel,i
	real param(mxparm),t,tend
        real*8 y(neq),yprime(neq),tft
	common/par/ m,g,alpha
	common/start/ ystart,ii
	real,dimension(11) :: ystart
	integer ii
	real m,g,alpha
	

	t=0E0
	do i=1,2
	   y(i)=0E0
	enddo
	y(1)=0E0
	y(2)=10E0
	y(3)=ystart(ii)
	y(4)=0E0       
	param(4)=10000
	param(10)=1.0   
	ido=1
	tel=0   
	tft=0E0
	endoftim=1   
	do i=1,100
	   tft=tft+endoftim/100. 
	   tend=tft
           call derivs(tft,y,yprime)
           call rk4(y,yprime,4,tft,0.01d0,y,derivs)
	enddo
	write(11,111) y(1),y(2),y(3),y(4)
		write(6,111) ystart(ii),y(1),y(2),y(3),y(4)
        ido=3
111     format(5E15.5)   	
        return
	end



! Numerical Recipes subroutine for fourth order Runge-Kutta integration:
!
! Given values for the variables y(1:n) and their derivatives dydx(i:n) 
! known at x, use the fourth-order Runge-Kutta method to advance the 
! solution over an interval h and return the incremented variables as 
! yout(1:n), which need not be a distinct array from y. The user supplies 
! the subroutine derivs(x,y,dydx), which returns derivatives dydx at x.

        subroutine rk4(y,dydx,n,x,h,yout,derivs)
        implicit none
        integer n,nmax,i
        double precision h,x,dydx(n),y(n),yout(n)
        external derivs
        parameter(nmax=50)
        double precision h6,hh,xh,dym(nmax),dyt(nmax),yt(nmax)
        hh=h*0.5d0
        h6=h/6d0
        xh=x+hh
        do i=1,n
           yt(i)=y(i)+hh*dydx(i)
        enddo
        call derivs(xh,yt,dyt)
        do i=1,n
           yt(i)=y(i)+hh*dyt(i)
        enddo
        call derivs(xh,yt,dym)
        do i=1,n
           yt(i)=y(i)+h*dym(i)
           dym(i)=dyt(i)+dym(i)
        enddo
        call derivs(x+h,yt,dyt)
        do i=1,n
           yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2d0*dym(i))
        enddo
        return
        end
