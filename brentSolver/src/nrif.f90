FUNCTION zbrent(func,x1,x2,tol) INTEGER ITMAX
  REAL zbrent,tol,x1,x2,func,EPS EXTERNAL func
  PARAMETER (ITMAX=100,EPS=3.e-8)
  !Using Brent’s method, find the root of a function func known to lie between x1 and x2. The root, returned as zbrent, will be refined until its accuracy is tol.
  !Parameters: Maximum allowed number of iterations, and machine floating-point precision.
  INTEGER iter
  REAL a,b,c,d,e,fa,fb,fc,p,q,r,
  s,tol1,xm
  a=x1
  b=x2
  fa=func(a)
  fb=func(b) if((fa.gt.0..and.fb.gt.0.).or.(fa.lt.0..and.fb.lt.0.))
  pause ’root must be bracketed for zbrent’ fc=fb
  iter=1,ITMAX if((fb.gt.0..and.fc.gt.0.).or.(fb.lt.0..and.fc.lt.0.))then
  c=a Rename a, b, c and adjust bounding interval d. fc=fa
  d=b-a
  e=d
endif if(abs(fc).lt.abs(fb)) then
a=b
b=c
c=a
fa=fb
fb=fc
fc=fa
endif
tol1=2.*EPS*abs(b)+0.5*tol xm=.5*(c-b)
if(abs(xm).le.tol1 .or. fb.eq.0.)then
zbrent=b
return endif
if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
  s=fb/fa Attempt inverse quadratic interpolation. if(a.eq.c) then
  p=2.*xm*s
q=1.-s else
q=fa/fc
r=fb/fc p=s*(2.*xm*q*(q-r)-(b-a)*(r-1.)) q=(q-1.)*(r-1.)*(s-1.)
c=b do 11
Convergence check.
Sample page from NUMERICAL RECIPES IN FORTRAN 77: THE ART OF SCIENTIFIC COMPUTING (ISBN 0-521-43064-X)
Copyright (C) 1986-1992 by Cambridge University Press. Programs Copyright (C) 1986-1992 by Numerical Recipes Software.
Permission is granted for internet users to make one paper copy for their own personal use. Further reproduction, or any copying of machine- readable files (including this one) to any server computer, is strictly prohibited. To order Numerical Recipes books, diskettes, or CDROMs visit website http://www.nr.com or call 1-800-872-7423 (North America only), or send email to trade@cup.cam.ac.uk (outside North America).
9.4 Newton-Raphson Method Using Derivative 355
endif
if(p.gt.0.) q=-q Check whether in bounds. p=abs(p)
if(2.*p .lt. min(3.*xm*q-abs(tol1*q),abs(e*q))) then
 e=d
d=p/q else
d=xm
e=d endif
else d=xm
e=d endif
a=b
fa=fb
if(abs(d) .gt. tol1) then
b=b+d else
b=b+sign(tol1,xm) endif
fb=func(b) enddo 11
Accept interpolation.
Interpolation failed, use bisection.
Bounds decreasing too slowly, use bisection.
Move last best guess to a. Evaluate new trial root.
pause ’zbrent exceeding maximum iterations’ zbrent=b
return
END
