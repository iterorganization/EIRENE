!******************************************************************
!> \brief Approximation of the exponential integral E_1
!>
!> Integral E_1(arg) from x=arg to x=infinity of exp(-x)/x
!> Here, Allen and Hastings approximation of the exponential integral
!> is used.
!> See Handbook of Mathematical Functions: With Formulas, Graphs, and
!> Mathematical Tables, by Milton Abramowitz, Irene A. Stegun, Courier
!> Corporation, 1964, ISBN 0486612724, 9780486612720
!> Formula 5.1.53 (0 <= arg <=1, |error| <= 2 x 10^-7)
!> Formula 5.1.56 (1 <= arg < inf, |error| <= 2 x 10^-8)
      function EIRENE_expint(arg)
      USE EIRMOD_PRECISION
      implicit none
      real(dp)::A0,A1,A2,A3,A4,A5
      real(dp)::   B1,B2,B3,B4
      real(dp)::   C1,C2,C3,C4
      REAL(DP),intent(in)::arg
      REAL(DP)::EIRENE_expint
      data A0/-.57721566_DP/,A1/.99999193_DP/,A2/-.24991055_DP/
      data A3/.05519968_DP/,A4/-.00976004_DP/,A5/.00107857_DP/
      data B1/8.5733287401_DP/,B2/18.0590169730_DP/,B3/8.6347608925_DP/
      data B4/.2677737343_DP/
      data C1/9.5733223454_DP/,C2/25.6329561486_DP/,C3/21.0996530827_DP/
      data C4/3.9584969228_DP/
      if(arg.gt.1.0) then
        EIRENE_expint=exp(-arg)*(B4+arg*(B3+arg*(B2+arg*(B1+arg))))/
     .         (C4+arg*(C3+arg*(C2+arg*(C1+arg))))
        EIRENE_expint=EIRENE_expint/arg
      else
        EIRENE_expint=-log(arg)+A0+
     .                arg*(A1+arg*(A2+arg*(A3+arg*(A4+arg*A5))))
      end if
      return
      end function EIRENE_expint
!****************************** END *******************************


      subroutine EIRENE_expi(x,ei,icon)
      USE EIRMOD_PRECISION
      IMPLICIT NONE
c  exponential integral function EI
c  -int exp(-x)/x, von -x nach unendlich     x<0,  identisch mit
c  +int exp(x)/x,  von -unendl. bis x
c
c   PV +int exp(-x)/x von unendl bis -x      x>0 ,  identisch mit
c   PV +int exp(x)/x von -unendl bis x
c
      INTEGER ICON, IER
      REAL(DP) EIRENE_mmdei,dei,ei,x,xx
      EXTERNAL EIRENE_mmdei
cdr   if (x.gt.0) then
        xx=x
        dei=-EIRENE_mmdei(xx,ier)
        icon=ier
        ei=dei
cdr   else
cdr     write (iunout,*) 'argument in expi lt.0, call exit'
cdr     call eirene_exit_own(1)
cdr   endif
      return
      end subroutine EIRENE_expi

c
      FUNCTION EIRENE_MMDEI (S,IER)
      USE EIRMOD_PRECISION
C  exponential integral function EI
c  Old IMSL routine. There: mmdei(ipot=2,s,ier), also:
c  s gt.0, mmdei is then: integral (s to infinity) exp(-t)/t dt
      IMPLICIT NONE
      REAL(DP) :: EIRENE_MMDEI,S,X,Y,Z,VALU,VAL
      INTEGER IER
      X=S
      Y=ABS(X)
      Z=0.25_DP*Y
      IF (Z.LE.1.0D0) THEN
       VALU=((((((((((((((((((((-.483702D-8*Z+.2685377D-7)*Z-
     1 .11703642D-6)*Z+.585911692D-6)*Z-.2843937873D-5)*Z+
     1 .1284394756D-4)*Z-.547380648948D-4)*Z+.21992775413732D-3)*Z-
     1 .8290046678016D-3)*Z+.0029187885699858)*Z-.0095523774824170)*Z+
     1 .028895941356815)*Z-.080266510032735)*Z+.20317460364863)*Z-
     1 .46439909294939)*Z+.9481481480886)*Z-1.706666666669)*Z+
     1 2.6666666666702)*Z-3.5555555555546)*Z+3.9999999999994)*Z-
     1 3.9999999999996)*Z+.57721566490143+LOG(Y)
       VALU=-VALU
      ELSE
       Z=1.0D0/Z
       VALU=(((((((((((((((((((-.77769007188383D-3*Z+.84295295893998D-2
     1 )*Z-.04272827083185)*Z+.13473041266261)*Z-.29671551148)*Z+
     1 .48618806480858)*Z-.61743468824936)*Z+.62656052544291)*Z-
     1 .52203502518727)*Z+.36771959879483)*Z-.22727998629908)*Z+
     1 .12965447884319)*Z-.72886262704894D-1)*Z+.043418637381012)*Z-
     1 .29244589614262D-1)*Z+.23433749581188D-1)*Z-.023437317333964)*Z+
     1 .03124999448124)*Z-.062499999910765)*Z+.24999999999935)*Z-
     1 .20933859981036D-14
c     if (y.gt.160.D0) then
c       valu=0
c     else
        VALU=VALU*EXP(-Y)
c     endif
      ENDIF
      VAL=VALU
      EIRENE_MMDEI=VAL
      IER=0
      RETURN
      END FUNCTION EIRENE_MMDEI

c
C  next:  FC11A is perhaps same routine as mmdei, but from another external code?

C
      SUBROUTINE EIRENE_EXPGD (X,Y)
      USE EIRMOD_PRECISION
      IMPLICIT NONE
      REAL(DP) :: X, Y, VALUE
      EXTERNAL :: EIRENE_FC11A
C
      CALL EIRENE_FC11A (VALUE,X)
      Y=-VALUE
C
      RETURN
      END SUBROUTINE EIRENE_EXPGD
C
c
      SUBROUTINE EIRENE_FC11A (VAL,S)
      USE EIRMOD_PRECISION
      IMPLICIT NONE
cdr  exponential integral EI, from Harwell Subroutine Library, 1973
C###### 15/05/70 LAST LIBRARY UPDATE
      REAL(DP) :: X,Y,Z,VALU,VAL,S
CDR   X=DBLE(S)
      X=S
      Y=DABS(X)
      Z=0.25D+0*Y
      IF (Z.LE.1.0D0) THEN
        VALU=
     1 ((((((((((((((((((((-.483702D-8*Z+.2685377D-7)*Z-.11703642D-6)
     1 *Z+.585911692D-6)*Z-.2843937873D-5)*Z+.1284394756D-4)*Z-.547380
     1 648948D-4)*Z+.21992775413732D-3)*Z-.8290046678016D-3)*Z+.00291878
     1 85699858D0)*Z-.0095523774824170D0)*Z+.028895941356815D0)*Z-.08026
     1 6510032735D0)*Z+.20317460364863D0)*Z-.46439909294939D0)*Z+.948148
     1 1480886D0)*Z-1.706666666669D0
     1      )*Z+2.6666666666702D0)*Z-3.5555555555546D0)*Z+3.999999999999
     1 4D0)*Z-3.9999999999996D0)*Z+.57721566490143D0+DLOG(Y)
        VALU=-VALU
      ELSE
        Z=1.0D0/Z
        VALU=
     1 (((((((((((((((((((-.77769007188383D-3*Z+.84295295893998D-2)
     1 *Z-.04272827083185D0)*Z+.13473041266261D0)*Z-.29671551148D0)*Z+
     1 .48618806480858D0)*Z-.61743468824936D0)*Z+.62656052544291D0)*Z-
     1 .52203502518727D0
     1 )*Z+.36771959879483D0)*Z-.22727998629908D0)*Z+.12965447884319D0)*
     1 Z-.72886262704894D-1)*Z+.043418637381012D0)*Z-.29244589614262D-1
     1 )*Z+.23433749581188D-1)*Z-.023437317333964D0)*Z+.03124999448124D0
     1 )*Z-.062499999910765D0)*Z+.24999999999935D0)*Z-
     1 .20933859981036D-14
        VALU=VALU*DEXP(-Y)
      END IF
CDR13 VAL=SNGL(VALU)
      VAL=VALU
      RETURN
      END SUBROUTINE EIRENE_FC11A
