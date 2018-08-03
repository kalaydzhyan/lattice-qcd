      subroutine measure
c
c     * This routine only calls subroutines for measurement.
c
c     * Calculate
c         SU(3) Wilson loops
c
c     * Include file
c         paravp3
c
c       S.Kitahara (96.10.15)
c----------------------------------------------------------------------

      include'paravp3'
      dimension fwt(n0,n0)

      call fwloop(fwt)

      write(3,10) ((fwt(i,j),j=i,n0),i=1,n0)
 10   format(10f12.7)

      return
      end
c*********************************************************************c
