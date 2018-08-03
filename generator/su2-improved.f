c*******************************************************************
c                     PROGRAM  "su2-improved.f"
c
c      Version : 01.2007  
c      the code was obtained as modification of the Wilson
c      action code by (V.K.M., M.M.-P., V.B., G.B.)
c---------------------------------------------------------------------    
c     MAIN SUBROUTINES:
c     START         :  defines the starting configuration.
c     UPDATE        :  to update field configuration (heat bath).
c     OVERUPDATE    :  overrelaxation for updating.
c     RANDOM_GT     :  performs random gauge transformation.
c     PLAQUETTE     :  to calculate plaq.
c     SAVE_CONF     :  to save MC configuration.
c     SAVE_RES      :  to save results.
c
c     iran  : random number seed: 4 integers, 1st has to be odd,
c             first 3 have to be < 4096 and 4th < 2048;
c             if init=3: iran(1)<0 means: keep old seed
c
c=======================================================================
      include 'su2lgfl.cmn'

      dimension iran(4)
ctrnc      dimension rcopy_save_sa(ncopy_max)

      character*1  c3
      character*2  c1,c2,c5,crun
      character*3  c4,crunl
      character*17 chmeas_fmt
      character*22 gl_fmt, gh_fmt, chmeas,chres
      character*23 chmeasl

c-----------------------------------------------------------------------
c     Choose parameters :
c-----------------------------------------------------------------------
      open(90,file='param.in',form='formatted',status='old')
      read(90,'(f14.8)') beta
      read(90,'(f14.8)') alfa
      read(90,'(e15.6)') epsilon
      read(90,'(e15.6)') eps_copy
      read(90,'(i14)')   init
      read(90,'(i14)')   nmeas
      read(90,'(i14)')   ntherm
      read(90,'(i14)')   nempty
      read(90,'(i14)')   noverupdate
      read(90,'(i14)')   itmax
      read(90,'(i14)')   nunit
      read(90,'(i14)')   ncopya
      read(90,'(i14)')   ncopyb
      read(90,'(i14)')   iflips
      read(90,'(f14.8)') tmax
      read(90,'(f14.8)') tmc
      read(90,'(f14.8)') tmin
      read(90,'(i14)')   niter
      read(90,'(e15.6)') cg_rest
      read(90,'(i14)')   nkmom
      read(90,'(i14)')   info
      read(90,'(i14)')   nor   
      read(90,'(i14)')   iran(1)
      read(90,'(i14)')   iran(2)
      read(90,'(i14)')   iran(3)
      read(90,'(i14)')   iran(4)
cfVB      read(90,'(i14)')   it_ctf1
      close(90)

c-----------------------------------------------------------------------
      if (nmeas.gt.nmeas_max) nmeas=nmeas_max

c-----------------------------------------------------------------------
c     To initialize and thermalize :
c-----------------------------------------------------------------------

      call neighbour
      call start(init,iran)

      if(nrun.lt.100) then
       write(crun,'(i2.2)') nrun
c     open(80,file='/storage/nfs001/ihep/user/bornvit/runinfo'//crun,
      open(80,file='runinfo'//crun,
     *        form='formatted',status='unknown',
     +     access='append')
      else
       write(crunl,'(i3.3)') nrun
c     open(80,file='/storage/nfs001/ihep/user/bornvit/runinfo'//crunl,
      open(80,file='runinfo'//crunl,
     *        form='formatted',status='unknown',
     +     access='append')
      endif

c-----------------------------------------------------------------------
c     Define auxiliary parameters :
c-----------------------------------------------------------------------
      int_beta=dint(beta)
      diff_beta=beta-dfloat(int_beta)
      nint_beta=(beta-dfloat(int_beta))*1000.0+0.0000001

      write(c1,'(i2.2)') lsize
      write(c2,'(i2.2)') lsizet
      write(c3,'(i1)')   int_beta
      write(c4,'(i3.3)') nint_beta
      write(c5,'(i2.2)') nrun
      write(crun,'(i2.2)') nrun

c-----------------------------------------------------------------------
      if(nrun.lt.100) then
       write(crun,'(i2.2)') nrun
       chmeas='plaq_'//c1//'x'//c2//'_bet'//c3//'p'//c4//'_'//crun
c      open(10,file='/storage/nfs001/ihep/user/bornvit/'//chmeas,
       open(10,file=chmeas,
     *     form='formatted',status='unknown')
       write(10,'(a,a)')'#',chmeas
       chmeas='ploo_'//c1//'x'//c2//'_bet'//c3//'p'//c4//'_'//crun
       open(11,file=chmeas,
     *     form='formatted',status='unknown')
       write(11,'(a,a)')'#',chmeas
      else
       write(crunl,'(i3.3)') nrun
       chmeasl='plaq_'//c1//'x'//c2//'_bet'//c3//'p'//c4//'_'//crunl
c      open(10,file='/storage/nfs001/ihep/user/bornvit/'//chmeasl,
       open(10,file=chmeasl,
     *     form='formatted',status='unknown')
       write(10,'(a,a)')'#',chmeasl
       chmeasl='ploo_'//c1//'x'//c2//'_bet'//c3//'p'//c4//'_'//crunl
       open(11,file=chmeasl,form='formatted',status='unknown')
       write(11,'(a,a)')'#',chmeasl
      endif
 
      write(10,'(a)')       '#'
      write(10,'(a)')       '#     Pure gauge SU(2)'
      write(10,'(a)')       '#     ****************'
      write(10,'(a,i4)')    '#        lsize =',lsize
      write(10,'(a,i4)')    '#       lsizet =',lsizet
      write(10,'(a,f8.4)')  '#         beta =',beta
      write(10,'(a)')       '#'
      write(10,'(a,i7)')    '#       ntherm =',ntherm
      write(10,'(a,i7)')    '#       nempty =',nempty
      write(10,'(a,i7)')    '#  noverupdate =',noverupdate
      write(10,'(a,i7)')    '#        nmeas =',nmeas
      write(10,'(a)')       '#'
      write(10,'(a)')       '#     ****************'

c     if(nrun.lt.100) then
c      write(crun,'(i2.2)') nrun
c      chmeas='SN_'//c1//'x'//c2//'_bet'//c3//'p'//c4//'_'//crun
c      open(17,file='/storage/nfs001/ihep/user/bornvit/'//chmeas,
c      open(17,file=chmeas,
c    *     form='formatted',status='unknown')
c      write(17,'(a,a)')'#',chmeas
c     else
c      write(crunl,'(i3.3)') nrun
c      chmeasl='SN_'//c1//'x'//c2//'_bet'//c3//'p'//c4//'_'//crunl
c      open(17,file='/storage/nfs001/ihep/user/bornvit/'//chmeasl,
c      open(17,file=chmeasl,
c    *     form='formatted',status='unknown')
c      write(17,'(a,a)')'#',chmeasl
c     endif

      
c     write(17,'(a)')       '#'
c     write(17,'(a)')       '#     Pure gauge SU(2)'
c     write(17,'(a)')       '#     ****************'
c     write(17,'(a,i4)')    '#        lsize =',lsize
c     write(17,'(a,i4)')    '#       lsizet =',lsizet
c     write(17,'(a,f8.4)')  '#         beta =',beta
c     write(17,'(a)')       '#'
c     write(17,'(a,i7)')    '#       iflips =',iflips
c     write(17,'(a,i7)')    '#         nrun =',nrun
c     write(17,'(a,i7)')    '#       ntherm =',ntherm
c     write(17,'(a,i7)')    '#       nempty =',nempty
c     write(17,'(a,i7)')    '#  noverupdate =',noverupdate
c     write(17,'(a,i7)')    '#        nmeas =',nmeas
c     write(17,'(a)')       '#'
c     write(17,'(a)')       '#     ****************'
      
cfVB      write(17,'(a)')'imeas icopy  itmx  slg_ct1  slg_fin'
c-----------------------------------------------------------------------
c     The main loop :
c-----------------------------------------------------------------------
      sum_itm=0.d0
      sum_itm2=0.d0
      sum_itm3=0.d0
      sum_itm4=0.d0

      do 900 imeas=1,nmeas
         nrun=nrun+1

         write(80,'(a,i5)') ' imeas =',imeas

         do iempty=1,nempty
            call update_impr
            do ioverupdate=1,noverupdate
               call overupdate_impr
            enddo
          do il=1,nli4
            anew(il)=aold(il)
          enddo
          call plaquette
          call ploop
          write(11,*) imeas,iempty,strt
         enddo
         do il=1,nli4
            anew(il)=aold(il)
         enddo
 
c-----------------------------------------------------------------------
      if(nrun.lt.100) then
       write(crun,'(i2.2)') nrun
c      open(12,file='/storage/nfs001/ihep/user/bornvit/CONFOUT'//crun,
       open(12,file='CON0'//crun,
     *     form='unformatted',status='unknown')
      else
       write(crunl,'(i3.3)') nrun
c      open(12,file='/storage/nfs001/ihep/user/bornvit/CONFOUT'//crunl,
       open(12,file='CON'//crunl,
     *     form='unformatted',status='unknown')
      endif
      write(12) aold
      write(12) nrun
      write(12) lsize,lsizet
      write(12) beta
      write(12) iran
      close(12)
 900  continue

c     itm_av=dint(sum_itm/dfloat(ntot))
c     write(17,'(a,i5)')'#itmax_av for RO w/o flip = ', itm_av      
c     itm_av=dint(sum_itm2/dfloat(ntot))
c     write(17,'(a,i5)')'#itmax_av for SA w/o flip = ', itm_av      
c     itm_av=dint(sum_itm3/dfloat(ntot))
c     write(17,'(a,i5)')'#itmax_av for RO with flip = ', itm_av      
c     itm_av=dint(sum_itm4/dfloat(ntot))
c     write(17,'(a,i5)')'#itmax_av for SA with flip = ', itm_av      
c     
c     close(17)
      close(80)

c-----------------------------------------------------------------------
      call save_conf(init)
c     call save_res(init)
c-----------------------------------------------------------------------
      stop
      end

c***********************************************************************
c***********************************************************************
      subroutine start(init,iran)

c This subroutine initializes the fields.
c Initialization is determined by INIT.
c
c INIT=1 : random gauge configuration
c INIT=2 : ordered gauge configuration
c INIT=3 : read previous configuration from file 'CONFIN'
c INIT=4 : reads some start configuration from file 'CON.LAT'
c=======================================================================
      include 'su2lgfl.cmn'

      dimension iran(4),iranr(4)
      data lran/0/
c-----------------------------------------------------------------------
c     Define certain parameters :
c-----------------------------------------------------------------------
      pi=dacos(-1.d0)
      tpi=2.*pi
      p2i=.5d0/pi

      if(iran(1).lt.0) lran=1
c-----------------------------------------------------------------------
      np1=16
      np2=np1*lsize
      np3=np2*lsize
      np4=np3*lsize

      isp1=1
      isp2=isp1*lsize
      isp3=isp2*lsize
      isp4=isp3*lsize

      lp1=np1*(1-lsize)
      lp2=np2*(1-lsize)
      lp3=np3*(1-lsize)
      lp4=np4*(1-lsizet)

      isl1=isp1*(1-lsize)
      isl2=isp2*(1-lsize)
      isl3=isp3*(1-lsize)
      isl4=isp4*(1-lsizet)
c-----------------------------------------------------------------------
c     Random start :
c-----------------------------------------------------------------------
      if(init.eq.1) then
         call setrn(iran)
         nrun=0

         do 11 i=1,nlinks
            il=4*i-3
            ao0=sranf()-0.5
            ao1=sranf()-0.5
            ao2=sranf()-0.5
            ao3=sranf()-0.5

            aa=dsqrt(ao0*ao0+ao1*ao1+ao2*ao2+ao3*ao3)
            aa1=1.d0/aa

            aold(il)  =ao0*aa1
            aold(il+1)=ao1*aa1
            aold(il+2)=ao2*aa1
            aold(il+3)=ao3*aa1

 11      continue
      endif
c-----------------------------------------------------------------------
c     Ordered start :
c-----------------------------------------------------------------------
      if(init.eq.2) then
         call setrn(iran)
         nrun=0

         do 12 i=1,nlinks
            il=4*i-3
            aold(il)  =1.d0
            aold(il+1)=0.d0
            aold(il+2)=0.d0
            aold(il+3)=0.d0
 12      continue

      endif
c-----------------------------------------------------------------------
c     Read configuration from file :
c-----------------------------------------------------------------------
      if(init.eq.3) then

         open(71,file='CONFIN',form='unformatted',status='old')
         read(71) aold
         read(71) nrun
         read(71) lsizer,lsizetr
         read(71) betar
         read(71) iranr
         close(71)

         if(lsizer.ne.lsize) then
            print 2001,lsizer,lsize
 2001      format(//,2x,'Lattice size incorrect! old=',2i5,', new=',2i5)
            stop 'Lattice size?'
         endif

         if(lsizetr.ne.lsizet) then
            print 2002,lsizetr,lsizet
 2002      format(//,2x,'Lattice size incorrect! old=',2i5,', new=',2i5)
            stop 'Lattice size?'
         endif

         if(betar.ne.beta) then
            print 2003,betar,beta
 2003       format(//,2x,'Beta incorrect! old=',f8.4,', new=',f8.4)
            stop 'Beta ?'
         endif
         
         call setrn(iranr)
c        nrun=nrun+1
         
      endif

      if(init.eq.4) then
         call setrn(iran)
         open(70,status='unknown',file='CON.LAT',form='unformatted')
         read(70) aold
         close(70)
      endif 

c-----------------------------------------------------------------------
c     Thermalization :
c-----------------------------------------------------------------------
      if(ntherm.gt.0.and.init.lt.3) then
         do itherm=1,ntherm
c           call update
c           print *, 'update starts'
            call update_impr
            do il=1,nli4
              anew(il)=aold(il)
            enddo
            call plaquette
c           print *, 'update finished'
            do ioverupdate=1,noverupdate
              call overupdate_impr
            enddo
c           do il=1,nli4
c             anew(il)=aold(il)
c           enddo
c           call plaquette
           if(mod(itherm,10).eq.0) write(80,'(a,i5)') ' itherm =',itherm
         enddo
         call save_conf(init)
      endif
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccc ANTISYMMETRIC TENSOR ccccccccccccccccccccccccccccccc

      nasym(1,2,3) = 1
      nasym(2,3,1) = 1
      nasym(3,1,2) = 1
      nasym(1,3,2) = -1
      nasym(3,2,1) = -1
      nasym(2,1,3) = -1		

ccccccccccccccccccc CONJ GRAD initialization cccccccccccccccccccccccccccccccc

      ncall = 0
      nit_max = 0
      nit_tot = 0	 
c-------------------------------------------------------------------------
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      return
      end
c***********************************************************************
c***********************************************************************
      subroutine update

c     the enumeration of sites is as follows:
c     ns{k}=(k1-1)*np1+(k2-1)*np2+(k3-1)*np3+(k4-1)*np4
c
c     the enumeration of links is as follows:
c     nl{k;mue}=ns{k}+4*(mue-1)+1
c
c     version: 19.03.1989
c=======================================================================
      include 'su2lgfl.cmn'

      dimension np(4),nn(4)
c-----------------------------------------------------------------------
      ns=-16
c-----------------------------------------------------------------------
c     Loops over  all sites are here:
c-----------------------------------------------------------------------
      do 504 i4=1,lsizet
      np(4)=np4
      nn(4)=np4
      if (i4.eq.lsizet) np(4)=lp4
      if (i4.eq.1)      nn(4)=lp4

      do 503 i3=1,lsize
      np(3)=np3
      nn(3)=np3
      if (i3.eq.lsize) np(3)=lp3
      if(i3 .eq. 1)    nn(3)=lp3

      do 502 i2=1,lsize
      np(2)=np2
      nn(2)=np2
      if (i2.eq.lsize) np(2)=lp2
      if (i2.eq.1)     nn(2)=lp2

      do 501 i1=1,lsize
      np(1)=np1
      nn(1)=np1
      if (i1.eq.lsize) np(1)=lp1
      if (i1.eq.1)     nn(1)=lp1
C-----------------------------------------------------------------------
c     here is the enumeration of sites:
c-----------------------------------------------------------------------
      ns=ns+16
c-----------------------------------------------------------------------
c     select link in mue-direction:
c-----------------------------------------------------------------------
      do 400 mue=1,4

      p0=0.0
      p1=0.0
      p2=0.0
      p3=0.0

      il=ns+4*(mue-1)+1
      a0=aold(il)
      a1=aold(il+1)
      a2=aold(il+2)
      a3=aold(il+3)
c-----------------------------------------------------------------------
c   compute the open plaquette in (mue,nue)-plane:
c-----------------------------------------------------------------------
      do 20 nue=1,4
      if (nue.eq.mue) go to 20

      l3=ns+4*(nue-1)+1
      l2=il+np(nue)
      l1=l3+np(mue)
c-----------------------------------------------------------------------
c      multiplication of three links in "positive" direction:
c-----------------------------------------------------------------------
      b0=aold(l1)
      b1=aold(l1+1)
      b2=aold(l1+2)
      b3=aold(l1+3)

      c0=aold(l2)
      c1=-aold(l2+1)
      c2=-aold(l2+2)
      c3=-aold(l2+3)

      bc0=b0*c0-b1*c1-b2*c2-b3*c3
      bc1=b0*c1+b1*c0-b2*c3+b3*c2
      bc2=b0*c2+b1*c3+b2*c0-b3*c1
      bc3=b0*c3-b1*c2+b2*c1+b3*c0

      d0=aold(l3)
      d1=-aold(l3+1)
      d2=-aold(l3+2)
      d3=-aold(l3+3)

      bd0=bc0*d0-bc1*d1-bc2*d2-bc3*d3
      bd1=bc0*d1+bc1*d0-bc2*d3+bc3*d2
      bd2=bc0*d2+bc1*d3+bc2*d0-bc3*d1
      bd3=bc0*d3-bc1*d2+bc2*d1+bc3*d0

      p0=p0+bd0
      p1=p1+bd1
      p2=p2+bd2
      p3=p3+bd3

      l3=ns-nn(nue)+4*(nue-1)+1
      l2=il-nn(nue)
      l1=l3+np(mue)
c-----------------------------------------------------------------------
c     multiplication of three links in "negative" direction:
c-----------------------------------------------------------------------
      b0=aold(l1)
      b1=-aold(l1+1)
      b2=-aold(l1+2)
      b3=-aold(l1+3)

      c0=aold(l2)
      c1=-aold(l2+1)
      c2=-aold(l2+2)
      c3=-aold(l2+3)

      bc0=b0*c0-b1*c1-b2*c2-b3*c3
      bc1=b0*c1+b1*c0-b2*c3+b3*c2
      bc2=b0*c2+b1*c3+b2*c0-b3*c1
      bc3=b0*c3-b1*c2+b2*c1+b3*c0

      d0=aold(l3)
      d1=aold(l3+1)
      d2=aold(l3+2)
      d3=aold(l3+3)

      bd0=bc0*d0-bc1*d1-bc2*d2-bc3*d3
      bd1=bc0*d1+bc1*d0-bc2*d3+bc3*d2
      bd2=bc0*d2+bc1*d3+bc2*d0-bc3*d1
      bd3=bc0*d3-bc1*d2+bc2*d1+bc3*d0

      p0=p0+bd0
      p1=p1+bd1
      p2=p2+bd2
      p3=p3+bd3

20    continue

      s0=beta*p0
      s1=beta*p1
      s2=beta*p2
      s3=beta*p3
c-----------------------------------------------------------------------
c     updating by the heat bath method:
c-----------------------------------------------------------------------
      d=dsqrt(s0*s0+s1*s1+s2*s2+s3*s3)
      dm1=1./d
      ed=dexp(d)
      emd=1./ed
      de=ed-emd
6     x0=sranf()*de+emd
      x0=dm1*dlog(x0)
      r2=1.-x0*x0
      if (r2.le.0) go to 6
      sqr=dsqrt(r2)
      if (sranf().gt.sqr) go to 6
7     x1=sranf()-.5
      x2=sranf()-.5
      v2=x1*x1+x2*x2
      if (v2.gt.0.25) go to 7
      x3=sranf()-.5
      v2=v2+x3*x3
      if (v2.gt.0.25) go to 7
      sn=sqr*dm1/dsqrt(v2)

      x0=x0*dm1
      x1=x1*sn
      x2=x2*sn
      x3=x3*sn
c-----------------------------------------------------------------------
c     it's the output of this subroutine:
c-----------------------------------------------------------------------
      aold(il)  = x0*s0+x1*s1+x2*s2+x3*s3
      aold(il+1)=-x0*s1+x1*s0+x2*s3-x3*s2
      aold(il+2)=-x0*s2-x1*s3+x2*s0+x3*s1
      aold(il+3)=-x0*s3+x1*s2-x2*s1+x3*s0

 350  continue
 400  continue
 501  continue
 502  continue
 503  continue
 504  continue
c-----------------------------------------------------------------------
c ++++++++++++++ Random Number 'Randomizer' +++++++++++++++++++
      jr=int(51.*sranf())+1
      do i=1,jr
      xx=sranf()
      enddo
c ++++++++++++++ Random Number 'Randomizer' +++++++++++++++++++
c-----------------------------------------------------------------------
      return
      end
C***********************************************************************
C***********************************************************************
      subroutine update_impr

c     the enumeration of sites is as follows:
c     ns{k}=(k1-1)*np1+(k2-1)*np2+(k3-1)*np3+(k4-1)*np4
c
c     the enumeration of links is as follows:
c     nl{k;mue}=ns{k}+4*(mue-1)+1
c
c     version: 01.2007
c=======================================================================
      include 'su2lgfl.cmn'

      dimension np(4),nn(4)
c-----------------------------------------------------------------------
      u0=0.89160                        
      gamma=1./(20.*u0*u0)
      ns=-16
      ns0=0
c-----------------------------------------------------------------------
c     Loops over  all sites are here:
c-----------------------------------------------------------------------
      do 504 i4=1,lsizet
      np(4)=np4
      nn(4)=np4
      if (i4.eq.lsizet) np(4)=lp4
      if (i4.eq.1)      nn(4)=lp4

      do 503 i3=1,lsize
      np(3)=np3
      nn(3)=np3
      if (i3.eq.lsize) np(3)=lp3
      if(i3 .eq. 1)    nn(3)=lp3

      do 502 i2=1,lsize
      np(2)=np2
      nn(2)=np2
      if (i2.eq.lsize) np(2)=lp2
      if (i2.eq.1)     nn(2)=lp2

      do 501 i1=1,lsize
      np(1)=np1
      nn(1)=np1
      if (i1.eq.lsize) np(1)=lp1
      if (i1.eq.1)     nn(1)=lp1
C-----------------------------------------------------------------------
c     here is the enumeration of sites:
c-----------------------------------------------------------------------
      ns=ns+16
      ns0=ns0+1
c-----------------------------------------------------------------------
c     select link in mue-direction:
c-----------------------------------------------------------------------
      do 400 mue=1,4

      p0=0.0
      p1=0.0
      p2=0.0
      p3=0.0
      pm0=0.0
      pm1=0.0
      pm2=0.0
      pm3=0.0

      il=ns+4*(mue-1)+1
      a0=aold(il)
      a1=aold(il+1)
      a2=aold(il+2)
      a3=aold(il+3)
c-----------------------------------------------------------------------
c   compute the open plaquette in (mue,nue)-plane:
c-----------------------------------------------------------------------
      do 20 nue=1,4
      if (nue.eq.mue) go to 20

c     print *, mue,nue,ns,ns0

      l3=ns+4*(nue-1)+1
      l2=il+np(nue)
      l1=l3+np(mue)
c-----------------------------------------------------------------------
c      multiplication of three links in "positive" direction:
c-----------------------------------------------------------------------
      b0=aold(l1)
      b1=aold(l1+1)
      b2=aold(l1+2)
      b3=aold(l1+3)

      c0=aold(l2)
      c1=-aold(l2+1)
      c2=-aold(l2+2)
      c3=-aold(l2+3)

      bc0=b0*c0-b1*c1-b2*c2-b3*c3
      bc1=b0*c1+b1*c0-b2*c3+b3*c2
      bc2=b0*c2+b1*c3+b2*c0-b3*c1
      bc3=b0*c3-b1*c2+b2*c1+b3*c0

      d0=aold(l3)
      d1=-aold(l3+1)
      d2=-aold(l3+2)
      d3=-aold(l3+3)

      bd0=bc0*d0-bc1*d1-bc2*d2-bc3*d3
      bd1=bc0*d1+bc1*d0-bc2*d3+bc3*d2
      bd2=bc0*d2+bc1*d3+bc2*d0-bc3*d1
      bd3=bc0*d3-bc1*d2+bc2*d1+bc3*d0

      p0=p0+bd0
      p1=p1+bd1
      p2=p2+bd2
      p3=p3+bd3


c 1st rectangular:
      l6=il-nn(mue)
      l5=l3-nn(mue)
      l4=l6+np(nue)

      d0=aold(l4)
      d1=-aold(l4+1)
      d2=-aold(l4+2)
      d3=-aold(l4+3)

      bd0=bc0*d0-bc1*d1-bc2*d2-bc3*d3
      bd1=bc0*d1+bc1*d0-bc2*d3+bc3*d2
      bd2=bc0*d2+bc1*d3+bc2*d0-bc3*d1
      bd3=bc0*d3-bc1*d2+bc2*d1+bc3*d0

      c0=aold(l5)
      c1=-aold(l5+1)
      c2=-aold(l5+2)
      c3=-aold(l5+3)

      bc0=bd0*c0-bd1*c1-bd2*c2-bd3*c3
      bc1=bd0*c1+bd1*c0-bd2*c3+bd3*c2
      bc2=bd0*c2+bd1*c3+bd2*c0-bd3*c1
      bc3=bd0*c3-bd1*c2+bd2*c1+bd3*c0
      
      d0=aold(l6)
      d1=aold(l6+1)
      d2=aold(l6+2)
      d3=aold(l6+3)

      bd0=bc0*d0-bc1*d1-bc2*d2-bc3*d3
      bd1=bc0*d1+bc1*d0-bc2*d3+bc3*d2
      bd2=bc0*d2+bc1*d3+bc2*d0-bc3*d1
      bd3=bc0*d3-bc1*d2+bc2*d1+bc3*d0

      pm0=pm0+bd0
      pm1=pm1+bd1
      pm2=pm2+bd2
      pm3=pm3+bd3

c 2nd rectangular:
      l7=il+np(mue)
      l9=l7+np(nue)
      ns1=nbp(mue,ns0)
      ns2=nbp(mue,ns1)
      l8=16*(ns2-1)+4*(nue-1)+1

      b0=aold(l7)
      b1=aold(l7+1)
      b2=aold(l7+2)
      b3=aold(l7+3)

      c0=aold(l8)
      c1=aold(l8+1)
      c2=aold(l8+2)
      c3=aold(l8+3)

      bc0=b0*c0-b1*c1-b2*c2-b3*c3
      bc1=b0*c1+b1*c0-b2*c3+b3*c2
      bc2=b0*c2+b1*c3+b2*c0-b3*c1
      bc3=b0*c3-b1*c2+b2*c1+b3*c0

      d0=aold(l9)
      d1=-aold(l9+1)
      d2=-aold(l9+2)
      d3=-aold(l9+3)

      bd0=bc0*d0-bc1*d1-bc2*d2-bc3*d3
      bd1=bc0*d1+bc1*d0-bc2*d3+bc3*d2
      bd2=bc0*d2+bc1*d3+bc2*d0-bc3*d1
      bd3=bc0*d3-bc1*d2+bc2*d1+bc3*d0

      c0=aold(l2)
      c1=-aold(l2+1)
      c2=-aold(l2+2)
      c3=-aold(l2+3)

      bc0=bd0*c0-bd1*c1-bd2*c2-bd3*c3
      bc1=bd0*c1+bd1*c0-bd2*c3+bd3*c2
      bc2=bd0*c2+bd1*c3+bd2*c0-bd3*c1
      bc3=bd0*c3-bd1*c2+bd2*c1+bd3*c0
      
      d0=aold(l3)
      d1=-aold(l3+1)
      d2=-aold(l3+2)
      d3=-aold(l3+3)

      bd0=bc0*d0-bc1*d1-bc2*d2-bc3*d3
      bd1=bc0*d1+bc1*d0-bc2*d3+bc3*d2
      bd2=bc0*d2+bc1*d3+bc2*d0-bc3*d1
      bd3=bc0*d3-bc1*d2+bc2*d1+bc3*d0

      pm0=pm0+bd0
      pm1=pm1+bd1
      pm2=pm2+bd2
      pm3=pm3+bd3

c 3rd rectangular:
      l10=l1+np(nue)
      ns1=nbp(nue,ns0)
      l12=16*(ns1-1)+4*(nue-1)+1
      ns2=nbp(nue,ns1)
      l11=16*(ns2-1)+4*(mue-1)+1

      b0=aold(l1)
      b1=aold(l1+1)
      b2=aold(l1+2)
      b3=aold(l1+3)

      c0=aold(l10)
      c1=aold(l10+1)
      c2=aold(l10+2)
      c3=aold(l10+3)

      bc0=b0*c0-b1*c1-b2*c2-b3*c3
      bc1=b0*c1+b1*c0-b2*c3+b3*c2
      bc2=b0*c2+b1*c3+b2*c0-b3*c1
      bc3=b0*c3-b1*c2+b2*c1+b3*c0

      d0=aold(l11)
      d1=-aold(l11+1)
      d2=-aold(l11+2)
      d3=-aold(l11+3)

      bd0=bc0*d0-bc1*d1-bc2*d2-bc3*d3
      bd1=bc0*d1+bc1*d0-bc2*d3+bc3*d2
      bd2=bc0*d2+bc1*d3+bc2*d0-bc3*d1
      bd3=bc0*d3-bc1*d2+bc2*d1+bc3*d0

      c0=aold(l12)
      c1=-aold(l12+1)
      c2=-aold(l12+2)
      c3=-aold(l12+3)

      bc0=bd0*c0-bd1*c1-bd2*c2-bd3*c3
      bc1=bd0*c1+bd1*c0-bd2*c3+bd3*c2
      bc2=bd0*c2+bd1*c3+bd2*c0-bd3*c1
      bc3=bd0*c3-bd1*c2+bd2*c1+bd3*c0
      
      d0=aold(l3)
      d1=-aold(l3+1)
      d2=-aold(l3+2)
      d3=-aold(l3+3)

      bd0=bc0*d0-bc1*d1-bc2*d2-bc3*d3
      bd1=bc0*d1+bc1*d0-bc2*d3+bc3*d2
      bd2=bc0*d2+bc1*d3+bc2*d0-bc3*d1
      bd3=bc0*d3-bc1*d2+bc2*d1+bc3*d0

      pm0=pm0+bd0
      pm1=pm1+bd1
      pm2=pm2+bd2
      pm3=pm3+bd3

c-----------------------------------------------------------------------
c     multiplication of three links in "negative" direction:
c-----------------------------------------------------------------------

      l3=ns-nn(nue)+4*(nue-1)+1
      l2=il-nn(nue)
      l1=l3+np(mue)

      b0=aold(l1)
      b1=-aold(l1+1)
      b2=-aold(l1+2)
      b3=-aold(l1+3)

      c0=aold(l2)
      c1=-aold(l2+1)
      c2=-aold(l2+2)
      c3=-aold(l2+3)

      bc0=b0*c0-b1*c1-b2*c2-b3*c3
      bc1=b0*c1+b1*c0-b2*c3+b3*c2
      bc2=b0*c2+b1*c3+b2*c0-b3*c1
      bc3=b0*c3-b1*c2+b2*c1+b3*c0

      d0=aold(l3)
      d1=aold(l3+1)
      d2=aold(l3+2)
      d3=aold(l3+3)

      bd0=bc0*d0-bc1*d1-bc2*d2-bc3*d3
      bd1=bc0*d1+bc1*d0-bc2*d3+bc3*d2
      bd2=bc0*d2+bc1*d3+bc2*d0-bc3*d1
      bd3=bc0*d3-bc1*d2+bc2*d1+bc3*d0

      p0=p0+bd0
      p1=p1+bd1
      p2=p2+bd2
      p3=p3+bd3

c 1st rectangular:
c     l6=il-nn(mue)
      l5=l3-nn(mue)
      l4=l6-nn(nue)

      d0=aold(l4)
      d1=-aold(l4+1)
      d2=-aold(l4+2)
      d3=-aold(l4+3)

      bd0=bc0*d0-bc1*d1-bc2*d2-bc3*d3
      bd1=bc0*d1+bc1*d0-bc2*d3+bc3*d2
      bd2=bc0*d2+bc1*d3+bc2*d0-bc3*d1
      bd3=bc0*d3-bc1*d2+bc2*d1+bc3*d0

      c0=aold(l5)
      c1=aold(l5+1)
      c2=aold(l5+2)
      c3=aold(l5+3)

      bc0=bd0*c0-bd1*c1-bd2*c2-bd3*c3
      bc1=bd0*c1+bd1*c0-bd2*c3+bd3*c2
      bc2=bd0*c2+bd1*c3+bd2*c0-bd3*c1
      bc3=bd0*c3-bd1*c2+bd2*c1+bd3*c0
      
      d0=aold(l6)
      d1=aold(l6+1)
      d2=aold(l6+2)
      d3=aold(l6+3)

      bd0=bc0*d0-bc1*d1-bc2*d2-bc3*d3
      bd1=bc0*d1+bc1*d0-bc2*d3+bc3*d2
      bd2=bc0*d2+bc1*d3+bc2*d0-bc3*d1
      bd3=bc0*d3-bc1*d2+bc2*d1+bc3*d0

      pm0=pm0+bd0
      pm1=pm1+bd1
      pm2=pm2+bd2
      pm3=pm3+bd3

c 2nd rectangular:
c     l7=il+np(mue)
      l9=l7-nn(nue)
      ns1=nbp(mue,ns0)
      ns2=nbp(mue,ns1)
      ns3=nbn(nue,ns2)
      l8=16*(ns3-1)+4*(nue-1)+1

      b0=aold(l7)
      b1=aold(l7+1)
      b2=aold(l7+2)
      b3=aold(l7+3)

      c0=aold(l8)
      c1=-aold(l8+1)
      c2=-aold(l8+2)
      c3=-aold(l8+3)

      bc0=b0*c0-b1*c1-b2*c2-b3*c3
      bc1=b0*c1+b1*c0-b2*c3+b3*c2
      bc2=b0*c2+b1*c3+b2*c0-b3*c1
      bc3=b0*c3-b1*c2+b2*c1+b3*c0

      d0=aold(l9)
      d1=-aold(l9+1)
      d2=-aold(l9+2)
      d3=-aold(l9+3)

      bd0=bc0*d0-bc1*d1-bc2*d2-bc3*d3
      bd1=bc0*d1+bc1*d0-bc2*d3+bc3*d2
      bd2=bc0*d2+bc1*d3+bc2*d0-bc3*d1
      bd3=bc0*d3-bc1*d2+bc2*d1+bc3*d0

      c0=aold(l2)
      c1=-aold(l2+1)
      c2=-aold(l2+2)
      c3=-aold(l2+3)

      bc0=bd0*c0-bd1*c1-bd2*c2-bd3*c3
      bc1=bd0*c1+bd1*c0-bd2*c3+bd3*c2
      bc2=bd0*c2+bd1*c3+bd2*c0-bd3*c1
      bc3=bd0*c3-bd1*c2+bd2*c1+bd3*c0
      
      d0=aold(l3)
      d1=aold(l3+1)
      d2=aold(l3+2)
      d3=aold(l3+3)

      bd0=bc0*d0-bc1*d1-bc2*d2-bc3*d3
      bd1=bc0*d1+bc1*d0-bc2*d3+bc3*d2
      bd2=bc0*d2+bc1*d3+bc2*d0-bc3*d1
      bd3=bc0*d3-bc1*d2+bc2*d1+bc3*d0

      pm0=pm0+bd0
      pm1=pm1+bd1
      pm2=pm2+bd2
      pm3=pm3+bd3

c 3rd rectangular:
      ns1=nbn(nue,ns0)
      ns2=nbn(nue,ns1)
      l12=16*(ns2-1)+4*(nue-1)+1
      l11=16*(ns2-1)+4*(mue-1)+1
      ns3=nbp(mue,ns2)
      l10=16*(ns3-1)+4*(nue-1)+1

      b0=aold(l1)
      b1=-aold(l1+1)
      b2=-aold(l1+2)
      b3=-aold(l1+3)

      c0=aold(l10)
      c1=-aold(l10+1)
      c2=-aold(l10+2)
      c3=-aold(l10+3)

      bc0=b0*c0-b1*c1-b2*c2-b3*c3
      bc1=b0*c1+b1*c0-b2*c3+b3*c2
      bc2=b0*c2+b1*c3+b2*c0-b3*c1
      bc3=b0*c3-b1*c2+b2*c1+b3*c0

      d0=aold(l11)
      d1=-aold(l11+1)
      d2=-aold(l11+2)
      d3=-aold(l11+3)

      bd0=bc0*d0-bc1*d1-bc2*d2-bc3*d3
      bd1=bc0*d1+bc1*d0-bc2*d3+bc3*d2
      bd2=bc0*d2+bc1*d3+bc2*d0-bc3*d1
      bd3=bc0*d3-bc1*d2+bc2*d1+bc3*d0

      c0=aold(l12)
      c1=aold(l12+1)
      c2=aold(l12+2)
      c3=aold(l12+3)

      bc0=bd0*c0-bd1*c1-bd2*c2-bd3*c3
      bc1=bd0*c1+bd1*c0-bd2*c3+bd3*c2
      bc2=bd0*c2+bd1*c3+bd2*c0-bd3*c1
      bc3=bd0*c3-bd1*c2+bd2*c1+bd3*c0
      
      d0=aold(l3)
      d1=aold(l3+1)
      d2=aold(l3+2)
      d3=aold(l3+3)

      bd0=bc0*d0-bc1*d1-bc2*d2-bc3*d3
      bd1=bc0*d1+bc1*d0-bc2*d3+bc3*d2
      bd2=bc0*d2+bc1*d3+bc2*d0-bc3*d1
      bd3=bc0*d3-bc1*d2+bc2*d1+bc3*d0

      pm0=pm0+bd0
      pm1=pm1+bd1
      pm2=pm2+bd2
      pm3=pm3+bd3

c-----------------------------------------------------------------------

20    continue

c     s0=beta*p0
c     s1=beta*p1
c     s2=beta*p2
c     s3=beta*p3
      s0=beta*(p0 - gamma*pm0)
      s1=beta*(p1 - gamma*pm1)
      s2=beta*(p2 - gamma*pm2)
      s3=beta*(p3 - gamma*pm3)
c-----------------------------------------------------------------------
c     updating by the heat bath method:
c-----------------------------------------------------------------------
      d=dsqrt(s0*s0+s1*s1+s2*s2+s3*s3)
      dm1=1./d
      ed=dexp(d)
      emd=1./ed
      de=ed-emd
6     x0=sranf()*de+emd
      x0=dm1*dlog(x0)
      r2=1.-x0*x0
      if (r2.le.0) go to 6
      sqr=dsqrt(r2)
      if (sranf().gt.sqr) go to 6
7     x1=sranf()-.5
      x2=sranf()-.5
      v2=x1*x1+x2*x2
      if (v2.gt.0.25) go to 7
      x3=sranf()-.5
      v2=v2+x3*x3
      if (v2.gt.0.25) go to 7
      sn=sqr*dm1/dsqrt(v2)

      x0=x0*dm1
      x1=x1*sn
      x2=x2*sn
      x3=x3*sn
c-----------------------------------------------------------------------
c     it's the output of this subroutine:
c-----------------------------------------------------------------------
      aold(il)  = x0*s0+x1*s1+x2*s2+x3*s3
      aold(il+1)=-x0*s1+x1*s0+x2*s3-x3*s2
      aold(il+2)=-x0*s2-x1*s3+x2*s0+x3*s1
      aold(il+3)=-x0*s3+x1*s2-x2*s1+x3*s0

 350  continue
 400  continue
 501  continue
 502  continue
 503  continue
 504  continue
c-----------------------------------------------------------------------
c ++++++++++++++ Random Number 'Randomizer' +++++++++++++++++++
      jr=int(51.*sranf())+1
      do i=1,jr
      xx=sranf()
      enddo
c ++++++++++++++ Random Number 'Randomizer' +++++++++++++++++++
c-----------------------------------------------------------------------
      return
      end
C***********************************************************************
C***********************************************************************
      subroutine overupdate_impr

c     the enumeration of sites is as follows:
c     ns{k}=(k1-1)*np1+(k2-1)*np2+(k3-1)*np3+(k4-1)*np4
c
c     the enumeration of links is as follows:
c     nl{k;mue}=ns{k}+4*(mue-1)+1
c
c     version: 01.2007
c=======================================================================
      include 'su2lgfl.cmn'

      dimension np(4),nn(4)
c-----------------------------------------------------------------------
      u0=0.90048
      gamma=1./(20.*u0*u0)
      ns=-16
      ns0=0
c-----------------------------------------------------------------------
c     Loops over  all sites are here:
c-----------------------------------------------------------------------
      do 504 i4=1,lsizet
      np(4)=np4
      nn(4)=np4
      if (i4.eq.lsizet) np(4)=lp4
      if (i4.eq.1)      nn(4)=lp4

      do 503 i3=1,lsize
      np(3)=np3
      nn(3)=np3
      if (i3.eq.lsize) np(3)=lp3
      if(i3 .eq. 1)    nn(3)=lp3

      do 502 i2=1,lsize
      np(2)=np2
      nn(2)=np2
      if (i2.eq.lsize) np(2)=lp2
      if (i2.eq.1)     nn(2)=lp2

      do 501 i1=1,lsize
      np(1)=np1
      nn(1)=np1
      if (i1.eq.lsize) np(1)=lp1
      if (i1.eq.1)     nn(1)=lp1
C-----------------------------------------------------------------------
c     here is the enumeration of sites:
c-----------------------------------------------------------------------
      ns=ns+16
      ns0=ns0+1
c-----------------------------------------------------------------------
c     select link in mue-direction:
c-----------------------------------------------------------------------
      do 400 mue=1,4

      p0=0.0
      p1=0.0
      p2=0.0
      p3=0.0
      pm0=0.0
      pm1=0.0
      pm2=0.0
      pm3=0.0

      il=ns+4*(mue-1)+1
      a0=aold(il)
      a1=aold(il+1)
      a2=aold(il+2)
      a3=aold(il+3)
c-----------------------------------------------------------------------
c   compute the open plaquette in (mue,nue)-plane:
c-----------------------------------------------------------------------
      do 20 nue=1,4
      if (nue.eq.mue) go to 20

c     print *, mue,nue,ns,ns0

      l3=ns+4*(nue-1)+1
      l2=il+np(nue)
      l1=l3+np(mue)
c-----------------------------------------------------------------------
c      multiplication of three links in "positive" direction:
c-----------------------------------------------------------------------
      b0=aold(l1)
      b1=aold(l1+1)
      b2=aold(l1+2)
      b3=aold(l1+3)

      c0=aold(l2)
      c1=-aold(l2+1)
      c2=-aold(l2+2)
      c3=-aold(l2+3)

      bc0=b0*c0-b1*c1-b2*c2-b3*c3
      bc1=b0*c1+b1*c0-b2*c3+b3*c2
      bc2=b0*c2+b1*c3+b2*c0-b3*c1
      bc3=b0*c3-b1*c2+b2*c1+b3*c0

      d0=aold(l3)
      d1=-aold(l3+1)
      d2=-aold(l3+2)
      d3=-aold(l3+3)

      bd0=bc0*d0-bc1*d1-bc2*d2-bc3*d3
      bd1=bc0*d1+bc1*d0-bc2*d3+bc3*d2
      bd2=bc0*d2+bc1*d3+bc2*d0-bc3*d1
      bd3=bc0*d3-bc1*d2+bc2*d1+bc3*d0

      p0=p0+bd0
      p1=p1+bd1
      p2=p2+bd2
      p3=p3+bd3


c 1st rectangular:
      l6=il-nn(mue)
      l5=l3-nn(mue)
      l4=l6+np(nue)

      d0=aold(l4)
      d1=-aold(l4+1)
      d2=-aold(l4+2)
      d3=-aold(l4+3)

      bd0=bc0*d0-bc1*d1-bc2*d2-bc3*d3
      bd1=bc0*d1+bc1*d0-bc2*d3+bc3*d2
      bd2=bc0*d2+bc1*d3+bc2*d0-bc3*d1
      bd3=bc0*d3-bc1*d2+bc2*d1+bc3*d0

      c0=aold(l5)
      c1=-aold(l5+1)
      c2=-aold(l5+2)
      c3=-aold(l5+3)

      bc0=bd0*c0-bd1*c1-bd2*c2-bd3*c3
      bc1=bd0*c1+bd1*c0-bd2*c3+bd3*c2
      bc2=bd0*c2+bd1*c3+bd2*c0-bd3*c1
      bc3=bd0*c3-bd1*c2+bd2*c1+bd3*c0
      
      d0=aold(l6)
      d1=aold(l6+1)
      d2=aold(l6+2)
      d3=aold(l6+3)

      bd0=bc0*d0-bc1*d1-bc2*d2-bc3*d3
      bd1=bc0*d1+bc1*d0-bc2*d3+bc3*d2
      bd2=bc0*d2+bc1*d3+bc2*d0-bc3*d1
      bd3=bc0*d3-bc1*d2+bc2*d1+bc3*d0

      pm0=pm0+bd0
      pm1=pm1+bd1
      pm2=pm2+bd2
      pm3=pm3+bd3

c 2nd rectangular:
      l7=il+np(mue)
      l9=l7+np(nue)
      ns1=nbp(mue,ns0)
      ns2=nbp(mue,ns1)
      l8=16*(ns2-1)+4*(nue-1)+1

      b0=aold(l7)
      b1=aold(l7+1)
      b2=aold(l7+2)
      b3=aold(l7+3)

      c0=aold(l8)
      c1=aold(l8+1)
      c2=aold(l8+2)
      c3=aold(l8+3)

      bc0=b0*c0-b1*c1-b2*c2-b3*c3
      bc1=b0*c1+b1*c0-b2*c3+b3*c2
      bc2=b0*c2+b1*c3+b2*c0-b3*c1
      bc3=b0*c3-b1*c2+b2*c1+b3*c0

      d0=aold(l9)
      d1=-aold(l9+1)
      d2=-aold(l9+2)
      d3=-aold(l9+3)

      bd0=bc0*d0-bc1*d1-bc2*d2-bc3*d3
      bd1=bc0*d1+bc1*d0-bc2*d3+bc3*d2
      bd2=bc0*d2+bc1*d3+bc2*d0-bc3*d1
      bd3=bc0*d3-bc1*d2+bc2*d1+bc3*d0

      c0=aold(l2)
      c1=-aold(l2+1)
      c2=-aold(l2+2)
      c3=-aold(l2+3)

      bc0=bd0*c0-bd1*c1-bd2*c2-bd3*c3
      bc1=bd0*c1+bd1*c0-bd2*c3+bd3*c2
      bc2=bd0*c2+bd1*c3+bd2*c0-bd3*c1
      bc3=bd0*c3-bd1*c2+bd2*c1+bd3*c0
      
      d0=aold(l3)
      d1=-aold(l3+1)
      d2=-aold(l3+2)
      d3=-aold(l3+3)

      bd0=bc0*d0-bc1*d1-bc2*d2-bc3*d3
      bd1=bc0*d1+bc1*d0-bc2*d3+bc3*d2
      bd2=bc0*d2+bc1*d3+bc2*d0-bc3*d1
      bd3=bc0*d3-bc1*d2+bc2*d1+bc3*d0

      pm0=pm0+bd0
      pm1=pm1+bd1
      pm2=pm2+bd2
      pm3=pm3+bd3

c 3rd rectangular:
      l10=l1+np(nue)
      ns1=nbp(nue,ns0)
      l12=16*(ns1-1)+4*(nue-1)+1
      ns2=nbp(nue,ns1)
      l11=16*(ns2-1)+4*(mue-1)+1

      b0=aold(l1)
      b1=aold(l1+1)
      b2=aold(l1+2)
      b3=aold(l1+3)

      c0=aold(l10)
      c1=aold(l10+1)
      c2=aold(l10+2)
      c3=aold(l10+3)

      bc0=b0*c0-b1*c1-b2*c2-b3*c3
      bc1=b0*c1+b1*c0-b2*c3+b3*c2
      bc2=b0*c2+b1*c3+b2*c0-b3*c1
      bc3=b0*c3-b1*c2+b2*c1+b3*c0

      d0=aold(l11)
      d1=-aold(l11+1)
      d2=-aold(l11+2)
      d3=-aold(l11+3)

      bd0=bc0*d0-bc1*d1-bc2*d2-bc3*d3
      bd1=bc0*d1+bc1*d0-bc2*d3+bc3*d2
      bd2=bc0*d2+bc1*d3+bc2*d0-bc3*d1
      bd3=bc0*d3-bc1*d2+bc2*d1+bc3*d0

      c0=aold(l12)
      c1=-aold(l12+1)
      c2=-aold(l12+2)
      c3=-aold(l12+3)

      bc0=bd0*c0-bd1*c1-bd2*c2-bd3*c3
      bc1=bd0*c1+bd1*c0-bd2*c3+bd3*c2
      bc2=bd0*c2+bd1*c3+bd2*c0-bd3*c1
      bc3=bd0*c3-bd1*c2+bd2*c1+bd3*c0
      
      d0=aold(l3)
      d1=-aold(l3+1)
      d2=-aold(l3+2)
      d3=-aold(l3+3)

      bd0=bc0*d0-bc1*d1-bc2*d2-bc3*d3
      bd1=bc0*d1+bc1*d0-bc2*d3+bc3*d2
      bd2=bc0*d2+bc1*d3+bc2*d0-bc3*d1
      bd3=bc0*d3-bc1*d2+bc2*d1+bc3*d0

      pm0=pm0+bd0
      pm1=pm1+bd1
      pm2=pm2+bd2
      pm3=pm3+bd3

c-----------------------------------------------------------------------
c     multiplication of three links in "negative" direction:
c-----------------------------------------------------------------------

      l3=ns-nn(nue)+4*(nue-1)+1
      l2=il-nn(nue)
      l1=l3+np(mue)

      b0=aold(l1)
      b1=-aold(l1+1)
      b2=-aold(l1+2)
      b3=-aold(l1+3)

      c0=aold(l2)
      c1=-aold(l2+1)
      c2=-aold(l2+2)
      c3=-aold(l2+3)

      bc0=b0*c0-b1*c1-b2*c2-b3*c3
      bc1=b0*c1+b1*c0-b2*c3+b3*c2
      bc2=b0*c2+b1*c3+b2*c0-b3*c1
      bc3=b0*c3-b1*c2+b2*c1+b3*c0

      d0=aold(l3)
      d1=aold(l3+1)
      d2=aold(l3+2)
      d3=aold(l3+3)

      bd0=bc0*d0-bc1*d1-bc2*d2-bc3*d3
      bd1=bc0*d1+bc1*d0-bc2*d3+bc3*d2
      bd2=bc0*d2+bc1*d3+bc2*d0-bc3*d1
      bd3=bc0*d3-bc1*d2+bc2*d1+bc3*d0

      p0=p0+bd0
      p1=p1+bd1
      p2=p2+bd2
      p3=p3+bd3

c 1st rectangular:
c     l6=il-nn(mue)
      l5=l3-nn(mue)
      l4=l6-nn(nue)

      d0=aold(l4)
      d1=-aold(l4+1)
      d2=-aold(l4+2)
      d3=-aold(l4+3)

      bd0=bc0*d0-bc1*d1-bc2*d2-bc3*d3
      bd1=bc0*d1+bc1*d0-bc2*d3+bc3*d2
      bd2=bc0*d2+bc1*d3+bc2*d0-bc3*d1
      bd3=bc0*d3-bc1*d2+bc2*d1+bc3*d0

      c0=aold(l5)
      c1=aold(l5+1)
      c2=aold(l5+2)
      c3=aold(l5+3)

      bc0=bd0*c0-bd1*c1-bd2*c2-bd3*c3
      bc1=bd0*c1+bd1*c0-bd2*c3+bd3*c2
      bc2=bd0*c2+bd1*c3+bd2*c0-bd3*c1
      bc3=bd0*c3-bd1*c2+bd2*c1+bd3*c0
      
      d0=aold(l6)
      d1=aold(l6+1)
      d2=aold(l6+2)
      d3=aold(l6+3)

      bd0=bc0*d0-bc1*d1-bc2*d2-bc3*d3
      bd1=bc0*d1+bc1*d0-bc2*d3+bc3*d2
      bd2=bc0*d2+bc1*d3+bc2*d0-bc3*d1
      bd3=bc0*d3-bc1*d2+bc2*d1+bc3*d0

      pm0=pm0+bd0
      pm1=pm1+bd1
      pm2=pm2+bd2
      pm3=pm3+bd3

c 2nd rectangular:
c     l7=il+np(mue)
      l9=l7-nn(nue)
      ns1=nbp(mue,ns0)
      ns2=nbp(mue,ns1)
      ns3=nbn(nue,ns2)
      l8=16*(ns3-1)+4*(nue-1)+1

      b0=aold(l7)
      b1=aold(l7+1)
      b2=aold(l7+2)
      b3=aold(l7+3)

      c0=aold(l8)
      c1=-aold(l8+1)
      c2=-aold(l8+2)
      c3=-aold(l8+3)

      bc0=b0*c0-b1*c1-b2*c2-b3*c3
      bc1=b0*c1+b1*c0-b2*c3+b3*c2
      bc2=b0*c2+b1*c3+b2*c0-b3*c1
      bc3=b0*c3-b1*c2+b2*c1+b3*c0

      d0=aold(l9)
      d1=-aold(l9+1)
      d2=-aold(l9+2)
      d3=-aold(l9+3)

      bd0=bc0*d0-bc1*d1-bc2*d2-bc3*d3
      bd1=bc0*d1+bc1*d0-bc2*d3+bc3*d2
      bd2=bc0*d2+bc1*d3+bc2*d0-bc3*d1
      bd3=bc0*d3-bc1*d2+bc2*d1+bc3*d0

      c0=aold(l2)
      c1=-aold(l2+1)
      c2=-aold(l2+2)
      c3=-aold(l2+3)

      bc0=bd0*c0-bd1*c1-bd2*c2-bd3*c3
      bc1=bd0*c1+bd1*c0-bd2*c3+bd3*c2
      bc2=bd0*c2+bd1*c3+bd2*c0-bd3*c1
      bc3=bd0*c3-bd1*c2+bd2*c1+bd3*c0
      
      d0=aold(l3)
      d1=aold(l3+1)
      d2=aold(l3+2)
      d3=aold(l3+3)

      bd0=bc0*d0-bc1*d1-bc2*d2-bc3*d3
      bd1=bc0*d1+bc1*d0-bc2*d3+bc3*d2
      bd2=bc0*d2+bc1*d3+bc2*d0-bc3*d1
      bd3=bc0*d3-bc1*d2+bc2*d1+bc3*d0

      pm0=pm0+bd0
      pm1=pm1+bd1
      pm2=pm2+bd2
      pm3=pm3+bd3

c 3rd rectangular:
      ns1=nbn(nue,ns0)
      ns2=nbn(nue,ns1)
      l12=16*(ns2-1)+4*(nue-1)+1
      l11=16*(ns2-1)+4*(mue-1)+1
      ns3=nbp(mue,ns2)
      l10=16*(ns3-1)+4*(nue-1)+1

      b0=aold(l1)
      b1=-aold(l1+1)
      b2=-aold(l1+2)
      b3=-aold(l1+3)

      c0=aold(l10)
      c1=-aold(l10+1)
      c2=-aold(l10+2)
      c3=-aold(l10+3)

      bc0=b0*c0-b1*c1-b2*c2-b3*c3
      bc1=b0*c1+b1*c0-b2*c3+b3*c2
      bc2=b0*c2+b1*c3+b2*c0-b3*c1
      bc3=b0*c3-b1*c2+b2*c1+b3*c0

      d0=aold(l11)
      d1=-aold(l11+1)
      d2=-aold(l11+2)
      d3=-aold(l11+3)

      bd0=bc0*d0-bc1*d1-bc2*d2-bc3*d3
      bd1=bc0*d1+bc1*d0-bc2*d3+bc3*d2
      bd2=bc0*d2+bc1*d3+bc2*d0-bc3*d1
      bd3=bc0*d3-bc1*d2+bc2*d1+bc3*d0

      c0=aold(l12)
      c1=aold(l12+1)
      c2=aold(l12+2)
      c3=aold(l12+3)

      bc0=bd0*c0-bd1*c1-bd2*c2-bd3*c3
      bc1=bd0*c1+bd1*c0-bd2*c3+bd3*c2
      bc2=bd0*c2+bd1*c3+bd2*c0-bd3*c1
      bc3=bd0*c3-bd1*c2+bd2*c1+bd3*c0
      
      d0=aold(l3)
      d1=aold(l3+1)
      d2=aold(l3+2)
      d3=aold(l3+3)

      bd0=bc0*d0-bc1*d1-bc2*d2-bc3*d3
      bd1=bc0*d1+bc1*d0-bc2*d3+bc3*d2
      bd2=bc0*d2+bc1*d3+bc2*d0-bc3*d1
      bd3=bc0*d3-bc1*d2+bc2*d1+bc3*d0

      pm0=pm0+bd0
      pm1=pm1+bd1
      pm2=pm2+bd2
      pm3=pm3+bd3

c-----------------------------------------------------------------------

20    continue

      p0=p0 - gamma*pm0
      p1=p1 - gamma*pm1
      p2=p2 - gamma*pm2
      p3=p3 - gamma*pm3

      pa0=p0*a0-p1*a1-p2*a2-p3*a3
      pa1=p0*a1+p1*a0-p2*a3+p3*a2
      pa2=p0*a2+p2*a0-p3*a1+p1*a3
      pa3=p0*a3+p3*a0-p1*a2+p2*a1

      pap0=pa0*p0-pa1*p1-pa2*p2-pa3*p3
      pap1=pa0*p1+pa1*p0-pa2*p3+pa3*p2
      pap2=pa0*p2+pa2*p0-pa3*p1+pa1*p3
      pap3=pa0*p3+pa3*p0-pa1*p2+pa2*p1

      pp=dsqrt(pap0*pap0+pap1*pap1+pap2*pap2+pap3*pap3)
      pp=1.d0/pp
c-----------------------------------------------------------------------
c     It's the output of this subroutine:
c-----------------------------------------------------------------------
      aold(il)  = pap0*pp
      aold(il+1)=-pap1*pp
      aold(il+2)=-pap2*pp
      aold(il+3)=-pap3*pp
c-----------------------------------------------------------------------
c     To check :
c-----------------------------------------------------------------------
      x0=a0*p0-a1*p1-a2*p2-a3*p3
      y0=(pap0*p0+pap1*p1+pap2*p2+pap3*p3)*pp
      dxy=dabs(x0-y0)
      if(dxy.gt.0.00001) stop 
c-----------------------------------------------------------------------
 400  continue
 501  continue
 502  continue
 503  continue
 504  continue
c-----------------------------------------------------------------------
      return
      end
C***********************************************************************
C***********************************************************************
      subroutine overupdate

c     the enumeration of sites is as follows:
c     ns{k}=(k1-1)*np1+(k2-1)*np2+(k3-1)*np3+(k4-1)*np4
c
c     the enumeration of links is as follows:
c     nl{k;mue}=ns{k}+4*(mue-1)+1
c
c     version: 18.02.2001
c=======================================================================
      include 'su2lgfl.cmn'

      dimension np(4),nn(4)
c-----------------------------------------------------------------------
      ns=-16
c-----------------------------------------------------------------------
c     Loops over  all sites are here:
c-----------------------------------------------------------------------
      do 504 i4=1,lsizet
      np(4)=np4
      nn(4)=np4
      if (i4.eq.lsizet) np(4)=lp4
      if (i4.eq.1)      nn(4)=lp4

      do 503 i3=1,lsize
      np(3)=np3
      nn(3)=np3
      if (i3.eq.lsize) np(3)=lp3
      if(i3 .eq. 1)    nn(3)=lp3

      do 502 i2=1,lsize
      np(2)=np2
      nn(2)=np2
      if (i2.eq.lsize) np(2)=lp2
      if (i2.eq.1)     nn(2)=lp2

      do 501 i1=1,lsize
      np(1)=np1
      nn(1)=np1
      if (i1.eq.lsize) np(1)=lp1
      if (i1.eq.1)     nn(1)=lp1

      ns=ns+16
c-----------------------------------------------------------------------
c     Select link in mue-direction:
c-----------------------------------------------------------------------
      do 400 mue=1,4

      p0=0.d0
      p1=0.d0
      p2=0.d0
      p3=0.d0

      il=ns+4*(mue-1)+1
      a0=aold(il)
      a1=aold(il+1)
      a2=aold(il+2)
      a3=aold(il+3)
c-----------------------------------------------------------------------
c     Compute the open plaquette in (mue,nue)-plane:
c-----------------------------------------------------------------------
      do 20 nue=1,4
      if (nue.eq.mue) go to 20

      l3=ns+4*(nue-1)+1
      l2=il+np(nue)
      l1=l3+np(mue)
c-----------------------------------------------------------------------
c     Multiplication of three links in "positive" direction:
c-----------------------------------------------------------------------
      b0=aold(l1)
      b1=aold(l1+1)
      b2=aold(l1+2)
      b3=aold(l1+3)

      c0=aold(l2)
      c1=-aold(l2+1)
      c2=-aold(l2+2)
      c3=-aold(l2+3)

      bc0=b0*c0-b1*c1-b2*c2-b3*c3
      bc1=b0*c1+b1*c0-b2*c3+b3*c2
      bc2=b0*c2+b1*c3+b2*c0-b3*c1
      bc3=b0*c3-b1*c2+b2*c1+b3*c0

      d0=aold(l3)
      d1=-aold(l3+1)
      d2=-aold(l3+2)
      d3=-aold(l3+3)

      bd0=bc0*d0-bc1*d1-bc2*d2-bc3*d3
      bd1=bc0*d1+bc1*d0-bc2*d3+bc3*d2
      bd2=bc0*d2+bc1*d3+bc2*d0-bc3*d1
      bd3=bc0*d3-bc1*d2+bc2*d1+bc3*d0

      p0=p0+bd0
      p1=p1+bd1
      p2=p2+bd2
      p3=p3+bd3

      l3=ns-nn(nue)+4*(nue-1)+1
      l2=il-nn(nue)
      l1=l3+np(mue)
c-----------------------------------------------------------------------
c     multiplication of three links in "negative" direction:
c-----------------------------------------------------------------------
      b0=aold(l1)
      b1=-aold(l1+1)
      b2=-aold(l1+2)
      b3=-aold(l1+3)

      c0=aold(l2)
      c1=-aold(l2+1)
      c2=-aold(l2+2)
      c3=-aold(l2+3)

      bc0=b0*c0-b1*c1-b2*c2-b3*c3
      bc1=b0*c1+b1*c0-b2*c3+b3*c2
      bc2=b0*c2+b1*c3+b2*c0-b3*c1
      bc3=b0*c3-b1*c2+b2*c1+b3*c0

      d0=aold(l3)
      d1=aold(l3+1)
      d2=aold(l3+2)
      d3=aold(l3+3)

      bd0=bc0*d0-bc1*d1-bc2*d2-bc3*d3
      bd1=bc0*d1+bc1*d0-bc2*d3+bc3*d2
      bd2=bc0*d2+bc1*d3+bc2*d0-bc3*d1
      bd3=bc0*d3-bc1*d2+bc2*d1+bc3*d0

      p0=p0+bd0
      p1=p1+bd1
      p2=p2+bd2
      p3=p3+bd3

20    continue

      pa0=p0*a0-p1*a1-p2*a2-p3*a3
      pa1=p0*a1+p1*a0-p2*a3+p3*a2
      pa2=p0*a2+p2*a0-p3*a1+p1*a3
      pa3=p0*a3+p3*a0-p1*a2+p2*a1

      pap0=pa0*p0-pa1*p1-pa2*p2-pa3*p3
      pap1=pa0*p1+pa1*p0-pa2*p3+pa3*p2
      pap2=pa0*p2+pa2*p0-pa3*p1+pa1*p3
      pap3=pa0*p3+pa3*p0-pa1*p2+pa2*p1

      pp=dsqrt(pap0*pap0+pap1*pap1+pap2*pap2+pap3*pap3)
      pp=1.d0/pp
c-----------------------------------------------------------------------
c     It's the output of this subroutine:
c-----------------------------------------------------------------------
      aold(il)  = pap0*pp
      aold(il+1)=-pap1*pp
      aold(il+2)=-pap2*pp
      aold(il+3)=-pap3*pp
c-----------------------------------------------------------------------
c     To check :
c-----------------------------------------------------------------------
      x0=a0*p0-a1*p1-a2*p2-a3*p3
      y0=(pap0*p0+pap1*p1+pap2*p2+pap3*p3)*pp
      dxy=dabs(x0-y0)
      if(dxy.gt.0.00001) stop 
c-----------------------------------------------------------------------
400   continue
501   continue
502   continue
503   continue
504   continue
c-----------------------------------------------------------------------
      return
      end
c***********************************************************************
c***********************************************************************
      subroutine random_gt

c     Generates matrices Omega_x and performs random gauge
c     transformations.
c
c     This subroutine does not generate the gauge transformations
c     equally distributed with respect to the Haar measure!
c
c=======================================================================
      include 'su2lgfl.cmn'

      dimension np(4),nn(4)
c-----------------------------------------------------------------------
      nsite1=0
      do 4 i4=1,lsizet
      np(4)=isp4
      nn(4)=isp4
      if (i4.eq.lsizet)  np(4)=isl4
      if (i4.eq.1)       nn(4)=isl4

      do 3 i3=1,lsize
      np(3)=isp3
      nn(3)=isp3
      if (i3.eq.lsize)  np(3)=isl3
      if (i3.eq.1)      nn(3)=isl3

      do 2 i2=1,lsize
      np(2)=isp2
      nn(2)=isp2
      if (i2.eq.lsize)  np(2)=isl2
      if (i2.eq.1)      nn(2)=isl2

      do 1 i1=1,lsize
      np(1)=isp1
      nn(1)=isp1
      if (i1.eq.lsize)  np(1)=isl1
      if (i1.eq.1)      nn(1)=isl1

      nsite1=nsite1+1
      ksl1=(nsite1-1)*16

      a0=sranf()-0.5
      a1=sranf()-0.5
      a2=sranf()-0.5
      a3=sranf()-0.5

      aa=dsqrt(a0*a0+a1*a1+a2*a2+a3*a3)
      aa1=1.d0/aa

      v0=a0*aa1
      v1=a1*aa1
      v2=a2*aa1
      v3=a3*aa1
c-----------------------------------------------------------------------
c     Gauge transformation :
c-----------------------------------------------------------------------
      do 14 j1=1,4
      nsite2=nsite1-nn(j1)

      ksl2=(nsite2-1)*16
      link1=ksl1+j1*4-3

      u0=anew(link1)
      u1=anew(link1+1)
      u2=anew(link1+2)
      u3=anew(link1+3)

      anew(link1)  =v0*u0-v1*u1-v2*u2-v3*u3
      anew(link1+1)=v0*u1+v1*u0-v2*u3+v3*u2
      anew(link1+2)=v0*u2+v1*u3+v2*u0-v3*u1
      anew(link1+3)=v0*u3-v1*u2+v2*u1+v3*u0

      link2=ksl2+j1*4-3

      u0=anew(link2)
      u1=anew(link2+1)
      u2=anew(link2+2)
      u3=anew(link2+3)

      anew(link2)  = u0*v0+u1*v1+u2*v2+u3*v3
      anew(link2+1)=-u0*v1+u1*v0+u2*v3-u3*v2
      anew(link2+2)=-u0*v2-u1*v3+u2*v0+u3*v1
      anew(link2+3)=-u0*v3+u1*v2-u2*v1+u3*v0

14    continue
1     continue
2     continue
3     continue
4     continue
c-----------------------------------------------------------------------
      return
      end
c***********************************************************************
c***********************************************************************
      subroutine plaquette

c     To calculate plaquette.
c=======================================================================
      include 'su2lgfl.cmn'

      dimension np(4)
c-----------------------------------------------------------------------
      sum_t=0.d0
      sum_s=0.d0
      do 504 i4=1,lsizet
         np(4)=np4
         if (i4.eq.lsizet) np(4)=lp4
         do 503 i3=1,lsize
            np(3)=np3
            if (i3.eq.lsize) np(3)=lp3
            do 502 i2=1,lsize
               np(2)=np2
               if (i2.eq.lsize) np(2)=lp2
               do 501 i1=1,lsize
                  np(1)=np1
                  if (i1.eq.lsize) np(1)=lp1
                  ns=(i1-1)*np1+(i2-1)*np2+(i3-1)*np3+(i4-1)*np4
c-----------------------------------------------------------------------
c     compute plaquette in (j1,j2)-plane :
c-----------------------------------------------------------------------
                  do 11 j1=1,3
                     do 12 j2=j1+1,4
                        
                        il1=ns+4*(j1-1)+1
                        il4=ns+4*(j2-1)+1
                        il3=il1+np(j2)
                        il2=il4+np(j1)
                        
                        a0=anew(il1)
                        a1=anew(il1+1)
                        a2=anew(il1+2)
                        a3=anew(il1+3)
                        
                        b0=anew(il2)
                        b1=anew(il2+1)
                        b2=anew(il2+2)
                        b3=anew(il2+3)
                        
                        c0=a0*b0-a1*b1-a2*b2-a3*b3
                        c1=a0*b1+a1*b0-a2*b3+a3*b2
                        c2=a0*b2+a1*b3+a2*b0-a3*b1
                        c3=a0*b3-a1*b2+a2*b1+a3*b0
                        
                        b0=anew(il3)
                        b1=-anew(il3+1)
                        b2=-anew(il3+2)
                        b3=-anew(il3+3)

                        a0=c0*b0-c1*b1-c2*b2-c3*b3
                        a1=c0*b1+c1*b0-c2*b3+c3*b2
                        a2=c0*b2+c1*b3+c2*b0-c3*b1
                        a3=c0*b3-c1*b2+c2*b1+c3*b0

                        b0=anew(il4)
                        b1=-anew(il4+1)
                        b2=-anew(il4+2)
                        b3=-anew(il4+3)

                        c0=a0*b0-a1*b1-a2*b2-a3*b3

                        if(j2.eq.4)then
                           sum_t=sum_t+c0
                        else
                           sum_s=sum_s+c0
                        endif
                        
 12                  continue
 11               continue
c-----------------------------------------------------------------------
 501           continue
 502        continue
 503     continue
 504  continue
c-----------------------------------------------------------------------
      pls=2.d0*sum_s/nplaqs
      plt=2.d0*sum_t/nplaqs
      write(80,*) '  pls=  ',pls,'  plt=  ', plt,'  pl=  ',(pls+plt)*.5
      write(10,*) pls, plt,(pls+plt)*.5
      return
      end
c***********************************************************************
c***********************************************************************
      subroutine save_conf(init)

c=======================================================================
      include 'su2lgfl.cmn'

      character*2  c1,c2,c5,crun
      character*3  crunl
      dimension iran(4),iranr(4)
c-----------------------------------------------------------------------
      call savern(iran)
      if(nrun.lt.100) then
       write(crun,'(i2.2)') nrun
c     write(c5,'(i2.2)') nrun
c      open(12,file='/storage/nfs001/ihep/user/bornvit/CONFOUT'//crun,
       open(12,file='CONFOUT'//crun,
     *     form='unformatted',status='unknown')
      else
       write(crunl,'(i3.3)') nrun
c      open(12,file='/storage/nfs001/ihep/user/bornvit/CONFOUT'//crunl,
       open(12,file='CONFOUT'//crunl,
     *     form='unformatted',status='unknown')
      endif
      write(12) aold
      write(12) nrun
      write(12) lsize,lsizet
      write(12) beta
      write(12) iran
      close(12)
c-----------------------------------------------------------------------
      return
      end
c***********************************************************************
      subroutine save_res(init)

c     To save results of measurements.
c=======================================================================
      include 'su2lgfl.cmn'

      dimension iran(4)

      character*1  c3
      character*2  c1,c2,c5,crun
      character*3  c4,crunl
      character*13 char1(3)
      character*17 chmeas_fmt
      character*22 gl_fmt, gh_fmt, chmeas,chres
      character*22 gl_fmtl, gh_fmtl, chmeasl,chresl

      data char1/' (hot start) ',' (cold start)',' (from disc) '/
c-----------------------------------------------------------------------
c     Define auxiliary parameters :
c-----------------------------------------------------------------------

      int_beta=dint(beta)
      diff_beta=beta-dfloat(int_beta)
      nint_beta=(beta-dfloat(int_beta))*1000.0+0.0000001

      write(c1,'(i2.2)') lsize
      write(c2,'(i2.2)') lsizet
      write(c3,'(i1)')   int_beta
      write(c4,'(i3.3)') nint_beta
      write(c5,'(i2.2)') nrun
c     write(crun,'(i2.2)') nrun

c-----------------------------------------------------------------------
c     Save results of the current run :
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     Create formatted files :
c-----------------------------------------------------------------------
c     Save measurements of GHOST PROPAGATOR :
c-----------------------------------------------------------------------
      

      if(info.gt.1)then

      if(nrun.lt.100) then
       write(crun,'(i2.2)') nrun
       gh_fmt='GH_'//c1//'x'//c2//'_bet'//c3//'p'//c4//'_'//crun
       open(15 ,file='/storage/nfs001/ihep/user/bornvit/'//gh_fmt,
     *        form='formatted',status='unknown')
       write(15,'(a,a)')     '#',   gh_fmt
      else
       write(crunl,'(i3.3)') nrun
       gh_fmtl='GH_'//c1//'x'//c2//'_bet'//c3//'p'//c4//'_'//crunl
       open(15 ,file='/storage/nfs001/ihep/user/bornvit/'//gh_fmtl,
     *        form='formatted',status='unknown')
       write(15,'(a,a)')     '#',   gh_fmtl
      endif

      write(15,'(a)')       '#'
      write(15,'(a)')       '#     Pure gauge SU(2)'
      write(15,'(a)')       '#     ****************'
      write(15,'(a,i4)')    '#        lsize =',lsize
      write(15,'(a,i4)')    '#       lsizet =',lsizet
      write(15,'(a,f8.4)')  '#         beta =',beta
      write(15,'(a)')       '#'
      write(15,'(a,i7)')    '#       iflips =',iflips
      write(15,'(a,i7)')    '#         nrun =',nrun
      write(15,'(a,i7)')    '#       ntherm =',ntherm
      write(15,'(a,i7)')    '#       nempty =',nempty
      write(15,'(a,i7)')    '#  noverupdate =',noverupdate
      write(15,'(a,i7)')    '#        nmeas =',nmeas
      write(15,'(a,i7)')    '#        ncopy =',ncopy
      write(15,'(a,i7)')    '#        ncopya=',ncopya
      write(15,'(a,i7)')    '#        ncopyb=',ncopyb
      write(15,'(a)')       '#'
      write(15,'(a)')       '#     ****************'
      write(15,'(a)')       '#'
      write(15,'(a)')
     *'#imeas  icopy  OR with flip        SA with flip' 

         do 250 imeas=1,nmeas
250      continue

         do 350 imeas=1,nmeas
         do 350 icopy=1,ncopy
c           if(info.gt.3)then
c      write(15,'(i5,2x,i5,2x,e14.6,2x,e14.6
c    *         )') imeas,icopy,
c    *ghpr_all_1(imeas,icopy,1),
c    *ghpr_all_2(imeas,icopy,1)
c           endif
350      continue
         close(15)
      endif
      
c-----------------------------------------------------------------------
c     Save measurements of CORREL :
c-----------------------------------------------------------------------

      if(info.gt.0)then

      if(nrun.lt.100) then
       write(crun,'(i2.2)') nrun
       gl_fmt='GL_'//c1//'x'//c2//'_bet'//c3//'p'//c4//'_'//crun
       open(16,file='/storage/nfs001/ihep/user/bornvit/'//gl_fmt,
     *        form='formatted',status='unknown')
       write(16,'(a,a)')     '#',   gl_fmt
      else
       write(crunl,'(i3.3)') nrun
       gl_fmtl='GL_'//c1//'x'//c2//'_bet'//c3//'p'//c4//'_'//crunl
       open(16,file='/storage/nfs001/ihep/user/bornvit/'//gl_fmtl,
     *        form='formatted',status='unknown')
       write(16,'(a,a)')     '#',   gl_fmtl
      endif

      write(16,'(a)')       '#'
      write(16,'(a)')       '#     Pure gauge SU(2)'
      write(16,'(a)')       '#     ****************'
      write(16,'(a,i4)')    '#        lsize =',lsize
      write(16,'(a,i4)')    '#       lsizet =',lsizet
      write(16,'(a,f8.4)')  '#         beta =',beta
      write(16,'(a)')       '#'
      write(16,'(a,i7)')    '#       iflips =',iflips
      write(16,'(a,i7)')    '#         nrun =',nrun
      write(16,'(a,i7)')    '#       ntherm =',ntherm
      write(16,'(a,i7)')    '#       nempty =',nempty
      write(16,'(a,i7)')    '#  noverupdate =',noverupdate
      write(16,'(a,i7)')    '#        nmeas =',nmeas
      write(16,'(a,i7)')    '#        ncopy =',ncopy
      write(16,'(a,i7)')    '#        ncopya=',ncopya
      write(16,'(a,i7)')    '#        ncopyb=',ncopyb
      write(16,'(a)')       '#'
      write(16,'(a)')       '#     ****************'
      write(16,'(a)')
     *'#imeas  icopy  OR with flip        SA with flip' 

c        do imeas=1,nmeas
c        do icopy=1,ncopy

c      write(16,'(i5,2x,i5,2x,e14.6,2x,e14.6
c    *         )') imeas,icopy,
c    *corr_t_all_1(imeas,icopy,1),
c    *corr_t_all_2(imeas,icopy,1)
c enddo
c        enddo

      write(16,'(a)')       '#     ****************'
      write(16,'(a)')
     *'#imeas  icopy  OR with flip        SA with flip' 

c        do imeas=1,nmeas
c        do icopy=1,ncopy
c      write(16,'(i5,2x,i5,2x,e14.6,2x,e14.6
c    *         )') imeas,icopy,
c    *corr_4_all_1(imeas,icopy,1),
c    *corr_4_all_2(imeas,icopy,1)
c enddo
c enddo

         close(16)

      endif
c-----------------------------------------------------------------------
      return
      end

c***********************************************************************
c***********************************************************************
c VRANF: vector and scalar random number generator
c
c        simplified version by Urs Heller.
c
      function sranf ()

      implicit double precision(a-h,o-z)


*
*  multiple-precision implementation of a linear congruential
*  random number generator. the sequence of random numbers is
*  identical to those produced by the cyber-205 random number
*  generator vranf.
*
*           x(n+1) = a * x(n)  (mod 2**47)
*  where
*           a = 4c65 da2c 866d
*
      integer ran_seed(0:3), ran_mult(0:3)
      common /vrndom/ ran_seed, ran_mult

      real     twom11, twom12
      parameter ( twom11 = 1/2048.0, twom12 = 1/4096.0 )
      integer  i0, i1, i2, i3

      i3 = ran_seed(3)*ran_mult(0) + ran_seed(2)*ran_mult(1)
     &   + ran_seed(1)*ran_mult(2) + ran_seed(0)*ran_mult(3)
      i2 = ran_seed(2)*ran_mult(0) + ran_seed(1)*ran_mult(1)
     &   + ran_seed(0)*ran_mult(2)
      i1 = ran_seed(1)*ran_mult(0) + ran_seed(0)*ran_mult(1)
      i0 = ran_seed(0)*ran_mult(0)

      ran_seed(0) = mod(i0, 4096)
      i1          = i1 + i0/4096
      ran_seed(1) = mod(i1, 4096)
      i2          = i2 + i1/4096
      ran_seed(2) = mod(i2, 4096)
      ran_seed(3) = mod(i3 + i2/4096, 2048)

      sranf = twom11 * ( real(ran_seed(3)) +
     &        twom12 * ( real(ran_seed(2)) +
     &        twom12 * ( real(ran_seed(1)) +
     &        twom12 * ( real(ran_seed(0)) ))))
      return
      end


c**************************************************************************
c**************************************************************************
c
c VRANF: Calculate an array of uniform random numbers.
c
c        On vector machines this could be vectorized.
c
      subroutine vranf ( rand, length)

      implicit double precision(a-h,o-z)


      integer  length

      dimension  rand(length)

      integer ran_seed(0:3), ran_mult(0:3)
      common /vrndom/ ran_seed, ran_mult

      real     twom11, twom12
      parameter ( twom11 = 1/2048.0, twom12 = 1/4096.0 )
      integer  i, i0, i1, i2, i3
c
      do i = 1, length
        i3 = ran_seed(3)*ran_mult(0) + ran_seed(2)*ran_mult(1)
     &     + ran_seed(1)*ran_mult(2) + ran_seed(0)*ran_mult(3)
        i2 = ran_seed(2)*ran_mult(0) + ran_seed(1)*ran_mult(1)
     &     + ran_seed(0)*ran_mult(2)
        i1 = ran_seed(1)*ran_mult(0) + ran_seed(0)*ran_mult(1)
        i0 = ran_seed(0)*ran_mult(0)

        ran_seed(0) = mod(i0, 4096)
        i1          = i1 + i0/4096
        ran_seed(1) = mod(i1, 4096)
        i2          = i2 + i1/4096
        ran_seed(2) = mod(i2, 4096)
        ran_seed(3) = mod(i3 + i2/4096, 2048)

        rand(i) = twom11 * ( real(ran_seed(3)) +
     &            twom12 * ( real(ran_seed(2)) +
     &            twom12 * ( real(ran_seed(1)) +
     &            twom12 * ( real(ran_seed(0)) ))))
      enddo

      return
      end


c**************************************************************************
c SETRN: Initialize the multipliers for the Cyber 205 style random
c        number generator.
c
      subroutine setrn ( lseed )

      implicit double precision(a-h,o-z)

      integer ran_seed(0:3), ran_mult(0:3)
      common /vrndom/ ran_seed, ran_mult

      integer  i, lseed(0:3), a(0:3)
      data  a/ 1645, 712, 1498, 1222 /

      if ( mod(lseed(0),2) .ne. 1 ) then
         stop 'setrn: seed seed(0) must be odd'
      endif

      if ( lseed(0) .ge. 4096 ) then
         stop 'setrn: seed seed(0) must be < 4096'
      endif

      if ( lseed(1) .ge. 4096 ) then
         stop 'setrn: seed seed(1) must be < 4096'
      endif

      if ( lseed(2) .ge. 4096 ) then
         stop 'setrn: seed seed(2) must be < 4096'
      endif

      if ( lseed(3) .ge. 2048 ) then
         stop 'setrn: seed seed(3) must be < 2048'
      endif

c+
c Set the seed and the multiplier info
c-
      do i = 0, 3
         ran_mult(i)     = a(i)
         ran_seed(i)     = lseed(i)
      enddo

      return
      end



c**************************************************************************
c SAVERN: Returns the random number seed that is stored on the front end.
c
      subroutine savern ( lseed )

      implicit double precision(a-h,o-z)

      integer ran_seed(0:3), ran_mult(0:3)
      common /vrndom/ ran_seed, ran_mult

      integer  lseed(0:3)
      integer  i

      do i = 0, 3
         lseed(i) = ran_seed(i)
      enddo

      return
      end
c***********************************************************************
c***********************************************************************
      subroutine avera(nconf,svect,itime,ave,sig,fs,ib,bv)

c      Analyses the data given in SVECT.
c     NCONF :  the number of the data
c     ITIME :  the time separation with corr. < .05
c     AVE   :  gives the average value
c     SIG   :  gives the naive error
c     FS    :  is the integral of the time relaxation
c     BV    :  is the bunched variance vector for up to ib
c              bunches. It is not computed if ib<2 .
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension svect(nconf),bv(ib)
c-----------------------------------------------------------------------
      ave=0.d0
      ave2=0.d0
      fs=0.d0
      itime=0
      do 2 i=1,nconf
      ave=ave+svect(i)
2     ave2=ave2+svect(i)**2
      ave=ave/nconf
      ave2=ave2/nconf
      sig=ave2-ave**2
      if(ave2.eq.0.d0) return
      if(sig/ave2.lt.1.d-10) return
c-----------------------------------------------------------------------
4     itime=itime+1
         if(itime.gt.nconf/3) return
      ax=0.d0
      ay=0.d0
      fi=0.d0
      nc=nconf-itime
      do 5 j=1,nc
      fi=fi+svect(j)*svect(j+itime)
      ax=ax+svect(j)
      ay=ay+svect(j+itime)
5     continue
      fi=(fi-ax*ay/nc)/nc/sig
      fs=fs+fi*nc/nconf
      if(fi.gt..05d0) goto 4
      bv(1)=dsqrt(sig)
      nis=1
      do 23 i=2,ib
      nis=2*nis
      ns=nconf/nis
      nm=nconf-ns*nis
      bav=0.d0
      bav2=0.d0
      do 22 k=1,ns
      sav=0.d0
      do 21 j=1,nis
      nm=nm+1
      sav=sav+svect(nm)
21    continue
      sav=sav/nis
      bav=bav+sav
      bav2=bav2+sav**2
22    continue
      bv(i)=dsqrt((bav2-bav**2/ns)/ns)
23    continue
      sig=dsqrt(sig/(nconf-1))
c-----------------------------------------------------------------------
      return
      end
c***********************************************************************
c***********************************************************************
      subroutine neighbour 
C                                      01.07.96
C**********************************************************************
C     THIS SUBROUTINE creates an array to find neighbour site
C     for any given site.
C**********************************************************************
c     include 'param.txt'
c     include 'proj.cmn'
      include 'su2lgfl.cmn'

      PARAMETER(IP1=1,IP2=IP1*LENGS,IP3=IP2*LENGS,IP4=IP3*LENGS)
      PARAMETER(IL1=IP1*(1-LENGS),IL2=IP2*(1-LENGS),
     *          IL3=IP3*(1-LENGS),IL4=IP4*(1-LENGT))

      dimension NP(4),NN(4)

      NS0=0

      DO 4 I4=1,LENGT
        NN(4)=IP4
        NP(4)=IP4
        IF (I4.EQ.LENGT) NP(4)=IL4
        IF (I4.EQ.1)     NN(4)=IL4
      DO 3 I3=1,LENGS
        NN(3)=IP3
        NP(3)=IP3
        IF (I3.EQ.LENGS) NP(3)=IL3
        IF (I3.EQ.1)     NN(3)=IL3
      DO 2 I2=1,LENGS
        NN(2)=IP2
        NP(2)=IP2
        IF (I2.EQ.LENGS) NP(2)=IL2
        IF (I2.EQ.1)     NN(2)=IL2
      DO 1 I1=1,LENGS
        NN(1)=IP1
        NP(1)=IP1
        IF (I1.EQ.LENGS) NP(1)=IL1
        IF (I1.EQ.1)     NN(1)=IL1

        NS0=NS0+1
        do mu=1,4
          NBP(mu,NS0)=NS0+NP(mu)
          NBN(mu,NS0)=NS0-NN(mu)
          NBP4(mu,NS0)=4*(NBP(mu,NS0)-1)
          NBN4(mu,NS0)=4*(NBN(mu,NS0)-1)
          NBP16(mu,NS0)=16*(NBP(mu,NS0)-1)
          NBN16(mu,NS0)=16*(NBN(mu,NS0)-1)
        enddo

  1    CONTINUE
  2    CONTINUE
  3    CONTINUE
  4    CONTINUE

       return
       end
