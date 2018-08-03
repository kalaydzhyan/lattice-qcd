C**********************************************************************
C*****************************************************************************
      SUBROUTINE ploop 
C              VERSION 12.2006
C TO COMPUTE POLYAKOV LOOP IN ALL DIRECTIONS
C*****************************************************************************
C            THE ENUMERATION OF SITES IS AS FOLLOWS:
C        NS{k}=(k1-1)*NP1+(k2-1)*NP2+(k3-1)*NP3+(k4-1)*NP4
C
C            THE ENUMERATION OF LINKS IS AS FOLLOWS:
C                 NL{k;mue}=NS{k}+4*(mue-1)+1
C*****************************************************************************

      include 'su2lgfl.cmn'
c     include 'param.txt'
c     include 'su2_mc.cmn'
c     PARAMETER (LENGS=20,LENGT=6)
c     PARAMETER (NSITES=LENGS**3*LENGT,NLINKS=4*NSITES,NVAR=4*NLINKS)
      PARAMETER (NVAR=4*NLINKS)
      PARAMETER (NSITE3=LENGS**3)
c     PARAMETER (NP1=16,NP2=NP1*LENGS,NP3=NP2*LENGS,NP4=NP3*LENGS)
c     PARAMETER (LP1=NP1*(1-LENGS),LP2=NP2*(1-LENGS),LP3=NP3*(1-LENGS),
c    *           LP4=NP4*(1-LENGT))
c     COMMON/YNGMIL/ ANEW(NVAR)


c     common /dataploop/ str
      real*8 strt

      DIMENSION NP(4)

      NSTS=3*LENGS**2*LENGT

      STRT=0.

C      NS=-16
C************************************************************
C ploop  IN 4-TH DIRECTION
C************************************************************

      DO 3 I3=1,LENGS
      NP(3)=NP3
      IF(I3.EQ.LENGS) NP(3)=LP3

      DO 2 I2=1,LENGS
      NP(2)=NP2
      IF(I2.EQ.LENGS) NP(2)=LP2

      DO 1 I1=1,LENGS
      NP(1)=NP1
      IF(I1.EQ.LENGS) NP(1)=LP1

      S0=1.
      S1=0.
      S2=0.
      S3=0.

      DO 4 I4=1,LENGT
      NP(4)=NP4
      IF(I4.EQ.LENGT) NP(4)=LP4

C*****************************************************************************
C          HERE IS THE ENUMERATION OF SITES:
C*****************************************************************************

      NS0=(I1-1)*NP1+(I2-1)*NP2+(I3-1)*NP3+(I4-1)*NP4
C      NS=NS+16
C*****************************************************************************
C          HERE IS THE CALCULATION OF NONABELIAN STRING:
C*****************************************************************************

      IL=NS0+13

      T0=ANEW(IL)
      T1=ANEW(IL+1)
      T2=ANEW(IL+2)
      T3=ANEW(IL+3)

      U0=S0*T0-S1*T1-S2*T2-S3*T3
      U1=S0*T1+S1*T0-S2*T3+S3*T2
      U2=S0*T2+S1*T3+S2*T0-S3*T1
      U3=S0*T3-S1*T2+S2*T1+S3*T0

      S0=U0
      S1=U1
      S2=U2
      S3=U3

 4    CONTINUE

      STRT=STRT+S0

 1    CONTINUE
 2    CONTINUE
 3    CONTINUE


      STRT=STRT/NSITE3


      RETURN
      END
C**********************************************************************
C**********************************************************************
