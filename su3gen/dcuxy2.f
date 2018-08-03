c*********************************************************************c
      subroutine dcuxy2
c
c     * for SU(3) Lattice on SX5.
c     --------------------------
c
c     * This routine converts /wrk1/ on an even-odd lattice into 
c       /var3/ on an ordinary lattice for the measurement of 
c       observables.
c
c     * Before this routine is called, you must once call 
c         Dir.
c
c     * Include file
c         parasx5
c
c     * Input arrays 
c         /var2/ : a configuration on an even-odd lattice.
c         /tb1/  : list vectors calculated by Dir.
c
c     * Output arrays 
c         /var3/ : the converted configuration.
c
c     * programmed by S.Kitahara.
c     * modified by S.Kitahara. (98.10.18)
c     * modified by Koma for sx5 (01.5.15)
c----------------------------------------------------------------------
      include 'paravp3'
      common /var2/ u(nvect,3,3,nd),v(nvect,3,3,nd)
      common /var3/ z(nsite,3,3,nd)
      common /tb1/  mioe(nsite), mioo(nsite)
c
      do l=1,nd
      do i=1,3
      do j=1,3
      do m=1,nsite
        if(mioe(m).eq.1) then
          z(m,i,j,l)= u(mioo(m),i,j,l)
        else
          z(m,i,j,l)= v(mioo(m),i,j,l)
        endif
      enddo
      enddo
      enddo
      enddo
c
      return
      end
