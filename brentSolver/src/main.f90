program main
  use brentSolvMod
  implicit none
  double precision, parameter :: xMin    = 0.d0
  double precision, parameter :: xMax    = 4.d0
  integer         , parameter :: iterMax = 100
  double precision, parameter :: crit    = 1.d-12
  double precision            :: x1, x2, ans
  
  ! ------------------------------------------------------ !
  ! --- [1] solve by bisection Method                  --- !
  ! ------------------------------------------------------ !
  x1  = 0.5d0
  x2  = 1.2d0
  ans = sqrt( 2.d0 )/2.d0
  write(6,*)
  write(6,*) "[Main@brentSolver] solve y=sin(x), where y=sqrt(2)/2"
  write(6,*) "[Main@brentSolver]     =>  ans :: x=pi/4, y=0"
  write(6,*) "[Main@brentSolver] Start Point :: x1, x2"
  write(6,*) x1, x2
  call brentSolver( x1, x2, iterMax, crit )

  write(6,*)
  write(6,*) "[Main@brentSolver] Solver End  :: x1, x2, residual=y=sin(x)-sqrt(2)/2 "
  write(6,*) x1, x2, analyticFunc( x1 )
  write(6,*)
  write(6,*) "[Main@brentSolver] Solved. x1, y=sin(x)"
  write(6,*) "[Main@brentSolver] Answer. xA, y=sin(xA)"
  write(6,*) x1, sin(x1)
  write(6,*) atan(1.d0), sin( atan(1.d0) )
  write(6,*) 
  
  return
end program main
