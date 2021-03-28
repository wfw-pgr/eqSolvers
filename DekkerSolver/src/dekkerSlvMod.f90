module dekkerSlvMod
contains

  ! ====================================================== !
  ! === analyticFunc                                   === !
  ! ====================================================== !
  function analyticFunc( xv )
    implicit none
    double precision, intent(in)  :: xv
    double precision              :: analyticFunc
    
    analyticFunc = sin( xv ) - 0.707d0
    return
  end function analyticFunc
  
  
  ! ====================================================== !
  ! === dekkerSolver                                   === !
  ! ====================================================== !
  subroutine dekkerSolver( x1, x2, iterMax, crit )
    implicit none
    integer         , intent(in)    :: iterMax
    double precision, intent(in)    :: crit
    double precision, intent(inout) :: x1, x2
    integer                         :: iter
    double precision                :: y1, y2
    double precision                :: sval, mval, cval, temp, eps

    ! ------------------------------------------------------ !
    ! --- [1] Main Loop                                  --- !
    ! ------------------------------------------------------ !
    eps  = crit
    y1   = analyticFunc( x1 )
    y2   = analyticFunc( x2 )
    do iter=1, iterMax
       ! -- end of loop judge                                -- !
       write(6,*) x1, y1, x2, y2
       if ( abs( y2 ).lt.crit ) exit
       ! -- calculates solution candidates                   -- !
       sval = ( x1*y2 - x2*y1  ) / ( y2 - y1 )
       mval = 0.5d0 * ( x1 + x2 )
       cval = x2
       ! -- choose closer solution :: mval or sval           -- !
       if ( abs( x2-mval ) < abs( x2-sval ) ) then
          x2 = mval          ! -- bisection method step -- !
       else
          x2 = sval          ! --    secant method step -- !
       endif
       ! -- too small step => at least, proceed eps          -- !
       if ( abs( x2-cval ).lt.eps ) then
          x2 = cval + sign( eps, x2-cval )
       endif
       y1   = analyticFunc( x1 )
       y2   = analyticFunc( x2 )
       ! -- if the new position x2 is in the same side as x1 -- !
       if ( ( y1*y2 ).gt.0.d0 ) then
          x1 = cval
       endif
       ! -- closer position is the main solution :: exchange -- !
       if ( abs(y1).lt.abs(y2) ) then
          temp = x1
          x1   = x2
          x2   = temp
          temp = y1
          y1   = y2
          y2   = temp
       endif
    enddo
    ! ------------------------------------------------------ !
    ! --- [2] Return Answer                              --- !
    ! ------------------------------------------------------ !
    x1 = x2
    x2 = y2
    
    return
  end subroutine dekkerSolver
  
  
end module dekkerSlvMod
