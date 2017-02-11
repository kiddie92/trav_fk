      SUBROUTINE findp0(x,p0)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     solve d(tau)/dp = 0 for p0, tau=p*x+eta*z
c     input:
c	 x --- distance
c	 p0 -- the largest possible p
c     output: p0
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IMPLICIT NONE
      REAL ZERO, dtdp0, x
      COMPLEX dtdp, p0, p1, p2
      ZERO = 1.E-7
      p1 = CMPLX(0.,AIMAG(p0))  !AIMAG复数的虚部
      DO WHILE ( p1.NE.p0 )
	 p2 = p0
	 p0 = 0.5*(p1+p2)
	 dtdp0 = dtdp(x,p0)
	 IF ( ABS(dtdp0).LT.ZERO .OR. p0.EQ.p1 .OR. p0.EQ.p2 ) RETURN
	 IF( dtdp0 .GT. 0. ) THEN
	    p1 = p0
	    p0 = p2
         END IF
      ENDDO
      RETURN
      END
 
      COMPLEX FUNCTION taup(p,x)
c define function tau(p) = p x + eta h
      IMPLICIT NONE
      INCLUDE'aseries.h'
      INTEGER i
      REAL x
      COMPLEX p, pp
      taup = p*x
      pp = p*p
      DO i = topp, bttm
	 taup=taup+SQRT(vps(1,i)-pp)*ray_len(1,i)
     &		  +SQRT(vps(2,i)-pp)*ray_len(2,i)
      ENDDO
      RETURN
      END
 
       COMPLEX FUNCTION dtdp(x,p)
c define d(tau)/dp
      IMPLICIT NONE
      INCLUDE'aseries.h'
      INTEGER j
      REAL x
      COMPLEX p, pp
      pp = p*p
      dtdp = 0.0
      DO j = topp, bttm
	 dtdp=dtdp-ray_len(1,j)/SQRT(vps(1,j)-pp)
     &		  -ray_len(2,j)/SQRT(vps(2,j)-pp)
      ENDDO
      dtdp = x + p*dtdp
      RETURN
      END
 
 