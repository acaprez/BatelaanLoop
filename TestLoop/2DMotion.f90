	program bragg
	use msimsl

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
	common/worksp/rwksp
	real rwksp(43592)
	call iwkin(43592)

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


	subroutine fcn(neq,t,y,yprime)
	implicit none
	common/par/ m,g,alpha
	real m,g,alpha
	integer neq
	real t,y(neq),yprime(neq)
	
	yprime(1)=y(2)   !y(1)=x, y(2)=v
	yprime(2)=0
	yprime(3)=y(4)
	yprime(4)=1.*sin(y(3))
	
	return                       
	end
	


	 
	
	subroutine result
	implicit none

	common/tole/ tol
	real tol

	real endoftim
	integer mxparm,neq
	parameter (mxparm=50,neq=4)
	integer ido,tel,i
	real fcn,param(mxparm),t,tend,tft,y(neq)
	common/par/ m,g,alpha
	common/start/ ystart,ii
	real,dimension(11) :: ystart
	integer ii
	real m,g,alpha
	
	external fcn,ivprk,sset

	t=0E0
	do i=1,2
	   y(i)=0E0
	enddo
	y(1)=0E0
	y(2)=10E0
	y(3)=ystart(ii)
	y(4)=0E0       
!	call sset(mxparm,0.0,param,1)
	param(4)=10000
	param(10)=1.0   
	ido=1
	tel=0   
	tft=0E0
	endoftim=1   
	do i=1,100
	   tft=tft+endoftim/100. 
	   tend=tft
	   call ivprk(ido,neq,fcn,t,tend,tol,param,y)
	enddo
	write(11,111) y(1),y(2),y(3),y(4)
		write(6,111) ystart(ii),y(1),y(2),y(3),y(4)
        ido=3
	call ivprk(ido,neq,fcn,t,tend,tol,param,y)
111     format(5E15.5)   	
        return
	end

