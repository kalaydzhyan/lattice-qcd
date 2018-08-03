c*********************************************************************c
      program vpsu3
c---------------------------------------------------------------------c
c     * This is a main routine to perform                             c
c       monte-carlo simulations of SU(3) lattice gauge theory         c
c       (in the maximally abelian gauge).                             c
c---------------------------------------------------------------------c
c
c     * Include file
c         paravp3
c
c     * The following subroutine is only available on VP and VPP.
c         clockv
c         gettod
c
c     * This routine calls
c         Seed, Dirct, Dir, inits, Monte, (Mafix1, Dcuxy,
c         Measure).
c
c     * Input parameters
c         b      : 6/(g*g). g is a coupling constant.
c         init   : 1 = start. 0 = continue.
c         ints   : 0 = hot start. 1 = cold start.
c         iter   : the number of MC sweeps for every measurement.
c         nav    : Iter * Nav is the number of thermalization sweeps.
c         niter  : the number of vacuum configuration.
c         timemax: the time limit (micro sec.).
c         nhit   : see subroutine Monte.
c         ngfmax : see subroutine Mafix1.
c         epsrzz : see subroutine Mafix1.
c         omega  : see subroutine Mafix1.
c       The file, inputvp3, describes a example of these parameters.
c
c       S.Kitahara (98.10.14)
c----------------------------------------------------------------------

      include'paravp3'

      parameter(nbet=30)
      parameter( nkei=256 )

      character*3 num(1000),numc
      common /rndk /k1,k2
      common /rndlp/irand(nkei,250)

      common /var3/ z(nsite,3,3,nd)
      common /var2/ u(nvect,3,3,nd),v(nvect,3,3,nd)
c     dimension am(nsite,3,3,nd),fi(nsite,3)
      dimension nspsw1(nbet)
c     data nspsw1/10,7*4,4*6,8,15,25,50,30,15,12*10/ !325
      data nspsw1/10,7*4,4*6,8,12,20,40,25,12,12*8/ !275

      real*8 time0,time1,time2

      data num /
     *    '000','001','002','003','004','005','006','007','008','009',
     *    '010','011','012','013','014','015','016','017','018','019',
     *    '020','021','022','023','024','025','026','027','028','029',
     *    '030','031','032','033','034','035','036','037','038','039',
     *    '040','041','042','043','044','045','046','047','048','049',
     *    '050','051','052','053','054','055','056','057','058','059',
     *    '060','061','062','063','064','065','066','067','068','069',
     *    '070','071','072','073','074','075','076','077','078','079',
     *    '080','081','082','083','084','085','086','087','088','089',
     *    '090','091','092','093','094','095','096','097','098','099',
     *    '100','101','102','103','104','105','106','107','108','109',
     *    '110','111','112','113','114','115','116','117','118','119',
     *    '120','121','122','123','124','125','126','127','128','129',
     *    '130','131','132','133','134','135','136','137','138','139',
     *    '140','141','142','143','144','145','146','147','148','149',
     *    '150','151','152','153','154','155','156','157','158','159',
     *    '160','161','162','163','164','165','166','167','168','169',
     *    '170','171','172','173','174','175','176','177','178','179',
     *    '180','181','182','183','184','185','186','187','188','189',
     *    '190','191','192','193','194','195','196','197','198','199',
     *    '200','201','202','203','204','205','206','207','208','209',
     *    '210','211','212','213','214','215','216','217','218','219',
     *    '220','221','222','223','224','225','226','227','228','229',
     *    '230','231','232','233','234','235','236','237','238','239',
     *    '240','241','242','243','244','245','246','247','248','249',
     *    '250','251','252','253','254','255','256','257','258','259',
     *    '260','261','262','263','264','265','266','267','268','269',
     *    '270','271','272','273','274','275','276','277','278','279',
     *    '280','281','282','283','284','285','286','287','288','289',
     *    '290','291','292','293','294','295','296','297','298','299',
     *    '300','301','302','303','304','305','306','307','308','309',
     *    '310','311','312','313','314','315','316','317','318','319',
     *    '320','321','322','323','324','325','326','327','328','329',
     *    '330','331','332','333','334','335','336','337','338','339',
     *    '340','341','342','343','344','345','346','347','348','349',
     *    '350','351','352','353','354','355','356','357','358','359',
     *    '360','361','362','363','364','365','366','367','368','369',
     *    '370','371','372','373','374','375','376','377','378','379',
     *    '380','381','382','383','384','385','386','387','388','389',
     *    '390','391','392','393','394','395','396','397','398','399',
     *    '400','401','402','403','404','405','406','407','408','409',
     *    '410','411','412','413','414','415','416','417','418','419',
     *    '420','421','422','423','424','425','426','427','428','429',
     *    '430','431','432','433','434','435','436','437','438','439',
     *    '440','441','442','443','444','445','446','447','448','449',
     *    '450','451','452','453','454','455','456','457','458','459',
     *    '460','461','462','463','464','465','466','467','468','469',
     *    '470','471','472','473','474','475','476','477','478','479',
     *    '480','481','482','483','484','485','486','487','488','489',
     *    '490','491','492','493','494','495','496','497','498','499',
     *    '500','501','502','503','504','505','506','507','508','509',
     *    '510','511','512','513','514','515','516','517','518','519',
     *    '520','521','522','523','524','525','526','527','528','529',
     *    '530','531','532','533','534','535','536','537','538','539',
     *    '540','541','542','543','544','545','546','547','548','549',
     *    '550','551','552','553','554','555','556','557','558','559',
     *    '560','561','562','563','564','565','566','567','568','569',
     *    '570','571','572','573','574','575','576','577','578','579',
     *    '580','581','582','583','584','585','586','587','588','589',
     *    '590','591','592','593','594','595','596','597','598','599',
     *    '600','601','602','603','604','605','606','607','608','609',
     *    '610','611','612','613','614','615','616','617','618','619',
     *    '620','621','622','623','624','625','626','627','628','629',
     *    '630','631','632','633','634','635','636','637','638','639',
     *    '640','641','642','643','644','645','646','647','648','649',
     *    '650','651','652','653','654','655','656','657','658','659',
     *    '660','661','662','663','664','665','666','667','668','669',
     *    '670','671','672','673','674','675','676','677','678','679',
     *    '680','681','682','683','684','685','686','687','688','689',
     *    '690','691','692','693','694','695','696','697','698','699',
     *    '700','701','702','703','704','705','706','707','708','709',
     *    '710','711','712','713','714','715','716','717','718','719',
     *    '720','721','722','723','724','725','726','727','728','729',
     *    '730','731','732','733','734','735','736','737','738','739',
     *    '740','741','742','743','744','745','746','747','748','749',
     *    '750','751','752','753','754','755','756','757','758','759',
     *    '760','761','762','763','764','765','766','767','768','769',
     *    '770','771','772','773','774','775','776','777','778','779',
     *    '780','781','782','783','784','785','786','787','788','789',
     *    '790','791','792','793','794','795','796','797','798','799',
     *    '800','801','802','803','804','805','806','807','808','809',
     *    '810','811','812','813','814','815','816','817','818','819',
     *    '820','821','822','823','824','825','826','827','828','829',
     *    '830','831','832','833','834','835','836','837','838','839',
     *    '840','841','842','843','844','845','846','847','848','849',
     *    '850','851','852','853','854','855','856','857','858','859',
     *    '860','861','862','863','864','865','866','867','868','869',
     *    '870','871','872','873','874','875','876','877','878','879',
     *    '880','881','882','883','884','885','886','887','888','889',
     *    '890','891','892','893','894','895','896','897','898','899',
     *    '900','901','902','903','904','905','906','907','908','909',
     *    '910','911','912','913','914','915','916','917','918','919',
     *    '920','921','922','923','924','925','926','927','928','929',
     *    '930','931','932','933','934','935','936','937','938','939',
     *    '940','941','942','943','944','945','946','947','948','949',
     *    '950','951','952','953','954','955','956','957','958','959',
     *    '960','961','962','963','964','965','966','967','968','969',
     *    '970','971','972','973','974','975','976','977','978','979',
     *    '980','981','982','983','984','985','986','987','988','989',
     *    '990','991','992','993','994','995','996','997','998','999'/

*      call gettod(time0)
      open(1,status='old',file='inputvp3')
      read(1,*) b
      read(1,*) init
      read(1,*) ints
      read(1,*) iter
      read(1,*) nav
      read(1,*) niter
      read(1,*) timemax
      read(1,*) nhit
      read(1,*) ngfmax
      read(1,*) epsrzz
      read(1,*) omega
      read(1,*) nrnd  
      close(1)

      open(2,status='unknown',file='su3_suz.res')
      open(3,status='unknown',file='wloop.dat')
      open(4,status='unknown',file='fmax.dat')
      open(6,status='unknown',file='mess.res')
      open(8,status='unknown',file='monop.dat')
      open(18,status='unknown',file='lclust.dat')
      write(2,300) n1,n2,n3,n4,b,init,ints,nav,iter,niter,nhit
  300 format(1x,'*******************************************'/
     c       1x,' Monte carlo simulation of SU3 LGT '/
     c       1x,'   '/
     c       1x,'*******************************************'/
     c       1x,'lattice size                :',i2,' x',i2,' x',i2
     c                                            ,' x',i2/
     c       1x,'beta                        :',f12.6/
     c       1x,'init (1=start, 0=continue)  :',i8/
     c       1x,'ints (0=hot, 1=cold)        :',i8/
     c       1x,'nav                         :',i8/
     c       1x,'iter                        :',i8/
     c       1x,'niter (# of samples)        :',i8/
     c       1x,'nhit                        :',i8/
     c       1x,'*******************************************')

c     preparation of random numbers, list vectors, constants

      k1=0
      k2=147
      call seed
      call dirct
      call dir
      pi2=8.e0*datan(1.d0)
      pi =4.e0*datan(1.d0)

c     generate a thermalized configuration

      if(init.eq.1)then
        call inits(ints,pi2)
        do iav=1,nav
          call monte(iter,nhit,trace,b,pi2)
          write(2,*)'nsweep=',iav
        enddo
      else
c       input file is written here.
      endif
      if(init.eq.0)then
        OPEN(1,file='CON.LAT',
     *      form='unformatted',status='unknown')
          read(1) z,k1,k2,irand
        CLOSE(1)
       call dcuxy3    ! z --> u,v
      endif

c     main loop start ------------------------------------------

c     do nnn = 1,niter
      do nnn = 1,niter
        numc=num(nnn)
        call monte(iter,nhit,trace,b,pi2)
c       if(iter.ne.0) then
c       OPEN(1,file='CON.LAT',form='unformatted',status='unknown')
c         write(1) u,v,k1,k2,irand
c       CLOSE(1)
c       endif

       call dcuxy2    ! u,v --> z
        OPEN(1,file='CON'//numc//'.LAT',form='unformatted',
     * status='unknown')
          write(1) z,k1,k2,irand
        CLOSE(1)
       call measure
      enddo

c     main loop end --------------------------------------------

  500 continue

   20   format('nnn=',i4,' igf=',i4,' trace=',f12.5)

      end
c*********************************************************************c
      subroutine monte(iter,nhit,trace,b,pi2)
c
c     * for SU(3) Lattice on VP.
c     --------------------------
c
c     * This routine generates a equilibrated configurations(/var2/).
c       using heat bath method.
c
c         action = 1-1/2 Tr\sum_{n,mu,nu} UUUU
c
c     * Before this routine is called, you must once call
c         Dirct.
c
c     * This routine calls
c         Rndpr3.
c
c     * Include file
c         paravp3
c
c     * Input arrays, parameters and constants.
c         /var2/ : input configuration.
c         /tabl/ : the list vectors which are returned by Dirct.
c         iter   : the number of sweeps.
c         nhit   : the number of trial in order to select a component
c                  of SU(2) matrix (usually =1).
c         b      : 6/(g*g). g is the coupling constant.
c         pi2    : 2 times pi (6.28318...).
c
c     * Output arrays and variables
c         /var2/ : updated configuration.
c         trace  : the value of the plaquette action.
c
c     * Work arrays
c         /ran1/ : random numbers generated by Rndpr3.
c
c     * programmed by KEK group.
c     * modified by S.Kitahara. (98.10.14)
c----------------------------------------------------------------------

      include'paravp3'

      common /var2/ u(nvect,3,3,nd),v(nvect,3,3,nd)
      common /ran1/ rm1(nvect),rm2(nvect),rm3(nvect),rm4(nvect)
      common /tabl /mo1(nvect,nd), mo2(nvect,nd),
     &              me1(nvect,nd), me2(nvect,nd),
     &              ldir2(nd-1,nd)

      common /wrk5 /w1(nvect,3,3,nd,nd),  w2(nvect,3,3,nd,nd),
     c              w3(nvect,3,3),w(nvect,3,3),
     c              a0(nvect),dxi2(nvect),delta(nvect)
      dimension     jhit(nvect)

      b3=b/3.e0
      hitp=0.e0

      do 1000 nit=1,iter
        ztr=(0.e0,0.e0)

c   ------- u  -------------------------------------------
      do nu = 2,nd
      do mu = 1,nu-1
      do j  = 1,3
      do i  = 1,3
        do m = 1,nvect
          m2 = mo1(m,mu)
          m4 = mo1(m,nu)
          m6 = mo2(m,nu)
          w1(m,i,j,mu,nu) = v(m2,i,1,nu)*dconjg(v(m4,j,1,mu))
     c                    + v(m2,i,2,nu)*dconjg(v(m4,j,2,mu))
     c                    + v(m2,i,3,nu)*dconjg(v(m4,j,3,mu))
          w2(m,i,j,mu,nu) = dconjg(v(m6,1,i,mu))*v(m6,1,j,nu)
     c                    + dconjg(v(m6,2,i,mu))*v(m6,2,j,nu)
     c                    + dconjg(v(m6,3,i,mu))*v(m6,3,j,nu)
        enddo
      enddo
      enddo
      enddo
      enddo
c
      do nu = 1,nd-1
      do mu = nu+1,nd
      do i  = 1,3
      do j  = 1,3
        do m = 1,nvect
          m5=me2(mo1(m,mu),nu)
          w1(m,i,j,mu,nu) = dconjg(w1(m,j,i,nu,mu))
          w3(m,i,j)       = dconjg(w2(m5,j,i,nu,mu))
        enddo
        do m = 1,nvect
          w2(m,i,j,mu,nu) = w3(m,i,j)
        enddo
      enddo
      enddo
      enddo
      enddo

      do 20 mu = 1,nd

      do j  = 1,3
      do i  = 1,3
        do m = 1,nvect
          w(m,i,j) = 0.e0
        enddo
        do nu0 = 1,nd-1
          nu = ldir2(nu0,mu)
          do m = 1,nvect
            m5=me2(mo1(m,mu),nu)
            w(m,i,j) = w(m,i,j)
     c               + w1(m,i,1,mu,nu)*dconjg(u(m,j,1,nu))
     c               + w1(m,i,2,mu,nu)*dconjg(u(m,j,2,nu))
     c               + w1(m,i,3,mu,nu)*dconjg(u(m,j,3,nu))
     c               + dconjg(u(m5,1,i,nu))*w2(m,1,j,mu,nu)
     c               + dconjg(u(m5,2,i,nu))*w2(m,2,j,mu,nu)
     c               + dconjg(u(m5,3,i,nu))*w2(m,3,j,mu,nu)
          enddo
        enddo
      enddo
      enddo
c - - - - - - - u(1,1)-u(2,2) - - - - - - - - - - - - - -
      do m=1,nvect
        jhit(m)=0
      enddo

      do i=1,2
      do j=1,2
      do m=1,nvect
        w3(m,i,j)=u(m,i,1,mu)*w(m,1,j)
     c           +u(m,i,2,mu)*w(m,2,j)
     c           +u(m,i,3,mu)*w(m,3,j)
      enddo
      enddo
      enddo

      do m=1,nvect
        ww11=w3(m,1,1)+dconjg(w3(m,2,2))
        ww12=w3(m,1,2)-dconjg(w3(m,2,1))
        ww21=w3(m,2,1)-dconjg(w3(m,1,2))
        ww22=w3(m,2,2)+dconjg(w3(m,1,1))
        dxi2(m)=dsqrt(dble(ww11*ww22-ww12*ww21))
      enddo

      do ihit1=1,nhit
        call rndpr3(4)
        do 29 m = 1,nvect
          if(jhit(m).eq.1) goto 29
          delta(m)=-(dlog(rm1(m))+dlog(rm2(m))*dcos(pi2*rm3(m))**2)
     c             /(dxi2(m)*b3)
          if (rm4(m)**2.gt.1.e0-delta(m)/2.e0) goto 29
            a0(m)=1.e0-delta(m)
            jhit(m)=1
  29    continue
      enddo

      call rndpr3(2)

      do m=1,nvect
        if (jhit(m).eq.1)then
          hitp=hitp+1.e0
          c0 =dble (w3(m,1,1)+w3(m,2,2))/dxi2(m)
          c1 =dimag(w3(m,1,2)+w3(m,2,1))/dxi2(m)
          c2 =dble (w3(m,1,2)-w3(m,2,1))/dxi2(m)
          c3 =dimag(w3(m,1,1)-w3(m,2,2))/dxi2(m)
          rad = 1.e0 - a0(m)**2
          a3 = dsqrt(rad)*(2.e0*rm1(m)-1.e0)
          rad2 = dsqrt(dabs(rad-a3**2))
          theta = pi2*rm2(m)
          a1 = rad2*dcos(theta)
          a2 = rad2*dsin(theta)
          d0 = a0(m)*c0+a1*c1+a2*c2+a3*c3
          d1 =-a0(m)*c1+a1*c0+a2*c3-a3*c2
          d2 =-a0(m)*c2+a2*c0+a3*c1-a1*c3
          d3 =-a0(m)*c3+a3*c0+a1*c2-a2*c1
          y11=dcmplx( d0, d3)
          y12=dcmplx( d2, d1)
          y21=dcmplx(-d2, d1)
          y22=dcmplx( d0,-d3)
          u11=u(m,1,1,mu)
          u12=u(m,1,2,mu)
          u13=u(m,1,3,mu)
          u21=u(m,2,1,mu)
          u22=u(m,2,2,mu)
          u23=u(m,2,3,mu)
          u(m,1,1,mu)=y11*u11+y12*u21
          u(m,1,2,mu)=y11*u12+y12*u22
          u(m,1,3,mu)=y11*u13+y12*u23
          u(m,2,1,mu)=y21*u11+y22*u21
          u(m,2,2,mu)=y21*u12+y22*u22
          u(m,2,3,mu)=y21*u13+y22*u23
        endif
      enddo

c - - - - - - - u(2,2)-u(3,3) - - - - - - - - - - - - - -
      do m=1,nvect
        jhit(m)=0
      enddo

      do i=2,3
      do j=2,3
      do m=1,nvect
        w3(m,i,j)=u(m,i,1,mu)*w(m,1,j)
     c           +u(m,i,2,mu)*w(m,2,j)
     c           +u(m,i,3,mu)*w(m,3,j)
      enddo
      enddo
      enddo

      do m=1,nvect
        ww22=w3(m,2,2)+dconjg(w3(m,3,3))
        ww23=w3(m,2,3)-dconjg(w3(m,3,2))
        ww32=w3(m,3,2)-dconjg(w3(m,2,3))
        ww33=w3(m,3,3)+dconjg(w3(m,2,2))
        dxi2(m)=dsqrt(dble(ww22*ww33-ww23*ww32))
      enddo

      do ihit1=1,nhit
        call rndpr3(4)
        do 129 m = 1,nvect
          if(jhit(m).eq.1) goto 129
          delta(m)=-(dlog(rm1(m))+dlog(rm2(m))*dcos(pi2*rm3(m))**2)
     c             /(dxi2(m)*b3)
          if (rm4(m)**2.gt.1.e0-delta(m)/2.e0) goto 129
            a0(m)=1.e0-delta(m)
            jhit(m)=1
 129    continue
      enddo

      call rndpr3(2)

      do m=1,nvect
        if (jhit(m).eq.1)then
          hitp=hitp+1.e0
          c0 =dble (w3(m,2,2)+w3(m,3,3))/dxi2(m)
          c1 =dimag(w3(m,2,3)+w3(m,3,2))/dxi2(m)
          c2 =dble (w3(m,2,3)-w3(m,3,2))/dxi2(m)
          c3 =dimag(w3(m,2,2)-w3(m,3,3))/dxi2(m)
          rad = 1.e0 - a0(m)**2
          a3 = dsqrt(rad)*(2.e0*rm1(m)-1.e0)
          rad2 = dsqrt(dabs(rad-a3**2))
          theta = pi2*rm2(m)
          a1 = rad2*dcos(theta)
          a2 = rad2*dsin(theta)
          d0 = a0(m)*c0+a1*c1+a2*c2+a3*c3
          d1 =-a0(m)*c1+a1*c0+a2*c3-a3*c2
          d2 =-a0(m)*c2+a2*c0+a3*c1-a1*c3
          d3 =-a0(m)*c3+a3*c0+a1*c2-a2*c1
          y22=dcmplx( d0, d3)
          y23=dcmplx( d2, d1)
          y32=dcmplx(-d2, d1)
          y33=dcmplx( d0,-d3)
          u21=u(m,2,1,mu)
          u22=u(m,2,2,mu)
          u23=u(m,2,3,mu)
          u31=u(m,3,1,mu)
          u32=u(m,3,2,mu)
          u33=u(m,3,3,mu)
          u(m,2,1,mu)=y22*u21+y23*u31
          u(m,2,2,mu)=y22*u22+y23*u32
          u(m,2,3,mu)=y22*u23+y23*u33
          u(m,3,1,mu)=y32*u21+y33*u31
          u(m,3,2,mu)=y32*u22+y33*u32
          u(m,3,3,mu)=y32*u23+y33*u33
        endif
      enddo

c - - - - - - - u(1,1)-u(3,3) - - - - - - - - - - - - - -
      do m=1,nvect
        jhit(m)=0
      enddo

      do i=1,3,2
      do j=1,3,2
      do m=1,nvect
        w3(m,i,j)=u(m,i,1,mu)*w(m,1,j)
     c           +u(m,i,2,mu)*w(m,2,j)
     c           +u(m,i,3,mu)*w(m,3,j)
      enddo
      enddo
      enddo

      do m=1,nvect
        ww11=w3(m,1,1)+dconjg(w3(m,3,3))
        ww13=w3(m,1,3)-dconjg(w3(m,3,1))
        ww31=w3(m,3,1)-dconjg(w3(m,1,3))
        ww33=w3(m,3,3)+dconjg(w3(m,1,1))
        dxi2(m)=dsqrt(dble(ww11*ww33-ww13*ww31))
      enddo

      do ihit1=1,nhit
        call rndpr3(4)
        do 229 m = 1,nvect
          if(jhit(m).eq.1) goto 229
          delta(m)=-(dlog(rm1(m))+dlog(rm2(m))*dcos(pi2*rm3(m))**2)
     c             /(dxi2(m)*b3)
          if (rm4(m)**2.gt.1.e0-delta(m)/2.e0) goto 229
            a0(m)=1.e0-delta(m)
            jhit(m)=1
 229    continue
      enddo

      call rndpr3(2)

      do m=1,nvect
        if (jhit(m).eq.1)then
          hitp=hitp+1.e0
          c0 =dble (w3(m,1,1)+w3(m,3,3))/dxi2(m)
          c1 =dimag(w3(m,1,3)+w3(m,3,1))/dxi2(m)
          c2 =dble (w3(m,1,3)-w3(m,3,1))/dxi2(m)
          c3 =dimag(w3(m,1,1)-w3(m,3,3))/dxi2(m)
          rad = 1.e0 - a0(m)**2
          a3 = dsqrt(rad)*(2.e0*rm1(m)-1.e0)
          rad2 = dsqrt(dabs(rad-a3**2))
          theta = pi2*rm2(m)
          a1 = rad2*dcos(theta)
          a2 = rad2*dsin(theta)
          d0 = a0(m)*c0+a1*c1+a2*c2+a3*c3
          d1 =-a0(m)*c1+a1*c0+a2*c3-a3*c2
          d2 =-a0(m)*c2+a2*c0+a3*c1-a1*c3
          d3 =-a0(m)*c3+a3*c0+a1*c2-a2*c1
          y11=dcmplx( d0, d3)
          y13=dcmplx( d2, d1)
          y31=dcmplx(-d2, d1)
          y33=dcmplx( d0,-d3)
          u11=u(m,1,1,mu)
          u12=u(m,1,2,mu)
          u13=u(m,1,3,mu)
          u31=u(m,3,1,mu)
          u32=u(m,3,2,mu)
          u33=u(m,3,3,mu)
          u(m,1,1,mu)=y11*u11+y13*u31
          u(m,1,2,mu)=y11*u12+y13*u32
          u(m,1,3,mu)=y11*u13+y13*u33
          u(m,3,1,mu)=y31*u11+y33*u31
          u(m,3,2,mu)=y31*u12+y33*u32
          u(m,3,3,mu)=y31*u13+y33*u33
        endif
      enddo


      do m=1,nvect
        znrm1=u(m,1,1,mu)*dconjg(u(m,1,1,mu))
     c       +u(m,1,2,mu)*dconjg(u(m,1,2,mu))
     c       +u(m,1,3,mu)*dconjg(u(m,1,3,mu))
*        znrm1=cdsqrt(znrm1)
        znrm1=sqrt(znrm1)
        u(m,1,1,mu)=u(m,1,1,mu)/znrm1
        u(m,1,2,mu)=u(m,1,2,mu)/znrm1
        u(m,1,3,mu)=u(m,1,3,mu)/znrm1
        zipd21=u(m,2,1,mu)*dconjg(u(m,1,1,mu))
     c        +u(m,2,2,mu)*dconjg(u(m,1,2,mu))
     c        +u(m,2,3,mu)*dconjg(u(m,1,3,mu))
        u(m,2,1,mu)=u(m,2,1,mu)-zipd21*u(m,1,1,mu)
        u(m,2,2,mu)=u(m,2,2,mu)-zipd21*u(m,1,2,mu)
        u(m,2,3,mu)=u(m,2,3,mu)-zipd21*u(m,1,3,mu)
        znrm2=u(m,2,1,mu)*dconjg(u(m,2,1,mu))
     c       +u(m,2,2,mu)*dconjg(u(m,2,2,mu))
     c       +u(m,2,3,mu)*dconjg(u(m,2,3,mu))
*        znrm2=cdsqrt(znrm2)
        znrm2=sqrt(znrm2)
        u(m,2,1,mu)=u(m,2,1,mu)/znrm2
        u(m,2,2,mu)=u(m,2,2,mu)/znrm2
        u(m,2,3,mu)=u(m,2,3,mu)/znrm2
        zipd31=u(m,3,1,mu)*dconjg(u(m,1,1,mu))
     c        +u(m,3,2,mu)*dconjg(u(m,1,2,mu))
     c        +u(m,3,3,mu)*dconjg(u(m,1,3,mu))
        zipd32=u(m,3,1,mu)*dconjg(u(m,2,1,mu))
     c        +u(m,3,2,mu)*dconjg(u(m,2,2,mu))
     c        +u(m,3,3,mu)*dconjg(u(m,2,3,mu))
        u(m,3,1,mu)=u(m,3,1,mu)-zipd32*u(m,2,1,mu)-zipd31*u(m,1,1,mu)
        u(m,3,2,mu)=u(m,3,2,mu)-zipd32*u(m,2,2,mu)-zipd31*u(m,1,2,mu)
        u(m,3,3,mu)=u(m,3,3,mu)-zipd32*u(m,2,3,mu)-zipd31*u(m,1,3,mu)
        znrm3=u(m,3,1,mu)*dconjg(u(m,3,1,mu))
     c       +u(m,3,2,mu)*dconjg(u(m,3,2,mu))
     c       +u(m,3,3,mu)*dconjg(u(m,3,3,mu))
*        znrm3=cdsqrt(znrm3)
        znrm3=sqrt(znrm3)
        u(m,3,1,mu)=u(m,3,1,mu)/znrm3
        u(m,3,2,mu)=u(m,3,2,mu)/znrm3
        u(m,3,3,mu)=u(m,3,3,mu)/znrm3
      enddo

      do i=1,3
      do m=1,nvect
        ztr=ztr+w(m,i,1)*u(m,1,i,mu)
     c         +w(m,i,2)*u(m,2,i,mu)
     c         +w(m,i,3)*u(m,3,i,mu)
      enddo
      enddo

  20  continue


c   ------- v  -------------------------------------------
      do nu = 2,nd
      do mu = 1,nu-1
      do j  = 1,3
      do i  = 1,3
        do m = 1,nvect
          m2 = me1(m,mu)
          m4 = me1(m,nu)
          m6 = me2(m,nu)
          w1(m,i,j,mu,nu) = u(m2,i,1,nu)*dconjg(u(m4,j,1,mu))
     c                    + u(m2,i,2,nu)*dconjg(u(m4,j,2,mu))
     c                    + u(m2,i,3,nu)*dconjg(u(m4,j,3,mu))
          w2(m,i,j,mu,nu) = dconjg(u(m6,1,i,mu))*u(m6,1,j,nu)
     c                    + dconjg(u(m6,2,i,mu))*u(m6,2,j,nu)
     c                    + dconjg(u(m6,3,i,mu))*u(m6,3,j,nu)
        enddo
      enddo
      enddo
      enddo
      enddo

      do nu = 1,nd-1
      do mu = nu+1,nd
      do i  = 1,3
      do j  = 1,3
        do m = 1,nvect
          m5=mo2(me1(m,mu),nu)
          w1(m,i,j,mu,nu) = dconjg(w1(m,j,i,nu,mu))
          w3(m,i,j)       = dconjg(w2(m5,j,i,nu,mu))
        enddo
        do m = 1,nvect
          w2(m,i,j,mu,nu) = w3(m,i,j)
        enddo
      enddo
      enddo
      enddo
      enddo

      do 320 mu = 1,nd

      do j  = 1,3
      do i  = 1,3
        do m = 1,nvect
          w(m,i,j) = 0.e0
        enddo
        do nu0 = 1,nd-1
          nu = ldir2(nu0,mu)
          do m = 1,nvect
            m5=mo2(me1(m,mu),nu)
            w(m,i,j) = w(m,i,j)
     c               + w1(m,i,1,mu,nu)*dconjg(v(m,j,1,nu))
     c               + w1(m,i,2,mu,nu)*dconjg(v(m,j,2,nu))
     c               + w1(m,i,3,mu,nu)*dconjg(v(m,j,3,nu))
     c               + dconjg(v(m5,1,i,nu))*w2(m,1,j,mu,nu)
     c               + dconjg(v(m5,2,i,nu))*w2(m,2,j,mu,nu)
     c               + dconjg(v(m5,3,i,nu))*w2(m,3,j,mu,nu)
          enddo
        enddo
      enddo
      enddo
c - - - - - - - v(1,1)-v(2,2) - - - - - - - - - - - - - -
      do m=1,nvect
        jhit(m)=0
      enddo

      do i=1,2
      do j=1,2
      do m=1,nvect
        w3(m,i,j)=v(m,i,1,mu)*w(m,1,j)
     c           +v(m,i,2,mu)*w(m,2,j)
     c           +v(m,i,3,mu)*w(m,3,j)
      enddo
      enddo
      enddo

      do m=1,nvect
        ww11=w3(m,1,1)+dconjg(w3(m,2,2))
        ww12=w3(m,1,2)-dconjg(w3(m,2,1))
        ww21=w3(m,2,1)-dconjg(w3(m,1,2))
        ww22=w3(m,2,2)+dconjg(w3(m,1,1))
        dxi2(m)=dsqrt(dble(ww11*ww22-ww12*ww21))
      enddo

      do ihit1=1,nhit
        call rndpr3(4)
        do 329 m = 1,nvect
          if(jhit(m).eq.1) goto 329
          delta(m)=-(dlog(rm1(m))+dlog(rm2(m))*dcos(pi2*rm3(m))**2)
     c             /(dxi2(m)*b3)
          if (rm4(m)**2.gt.1.e0-delta(m)/2.e0) goto 329
            a0(m)=1.e0-delta(m)
            jhit(m)=1
 329    continue
      enddo
c
      call rndpr3(2)
c
      do m=1,nvect
        if (jhit(m).eq.1)then
          hitp=hitp+1.e0
          c0 =dble (w3(m,1,1)+w3(m,2,2))/dxi2(m)
          c1 =dimag(w3(m,1,2)+w3(m,2,1))/dxi2(m)
          c2 =dble (w3(m,1,2)-w3(m,2,1))/dxi2(m)
          c3 =dimag(w3(m,1,1)-w3(m,2,2))/dxi2(m)
          rad = 1.e0 - a0(m)**2
          a3 = dsqrt(rad)*(2.e0*rm1(m)-1.e0)
          rad2 = dsqrt(dabs(rad-a3**2))
          theta = pi2*rm2(m)
          a1 = rad2*dcos(theta)
          a2 = rad2*dsin(theta)
          d0 = a0(m)*c0+a1*c1+a2*c2+a3*c3
          d1 =-a0(m)*c1+a1*c0+a2*c3-a3*c2
          d2 =-a0(m)*c2+a2*c0+a3*c1-a1*c3
          d3 =-a0(m)*c3+a3*c0+a1*c2-a2*c1
          y11=dcmplx( d0, d3)
          y12=dcmplx( d2, d1)
          y21=dcmplx(-d2, d1)
          y22=dcmplx( d0,-d3)
          v11=v(m,1,1,mu)
          v12=v(m,1,2,mu)
          v13=v(m,1,3,mu)
          v21=v(m,2,1,mu)
          v22=v(m,2,2,mu)
          v23=v(m,2,3,mu)
          v(m,1,1,mu)=y11*v11+y12*v21
          v(m,1,2,mu)=y11*v12+y12*v22
          v(m,1,3,mu)=y11*v13+y12*v23
          v(m,2,1,mu)=y21*v11+y22*v21
          v(m,2,2,mu)=y21*v12+y22*v22
          v(m,2,3,mu)=y21*v13+y22*v23
        endif
      enddo

c - - - - - - - v(2,2)-v(3,3) - - - - - - - - - - - - - -
      do m=1,nvect
        jhit(m)=0
      enddo

      do i=2,3
      do j=2,3
      do m=1,nvect
        w3(m,i,j)=v(m,i,1,mu)*w(m,1,j)
     c           +v(m,i,2,mu)*w(m,2,j)
     c           +v(m,i,3,mu)*w(m,3,j)
      enddo
      enddo
      enddo

      do m=1,nvect
        ww22=w3(m,2,2)+dconjg(w3(m,3,3))
        ww23=w3(m,2,3)-dconjg(w3(m,3,2))
        ww32=w3(m,3,2)-dconjg(w3(m,2,3))
        ww33=w3(m,3,3)+dconjg(w3(m,2,2))
        dxi2(m)=dsqrt(dble(ww22*ww33-ww23*ww32))
      enddo

      do ihit1=1,nhit
        call rndpr3(4)
        do 429 m = 1,nvect
          if(jhit(m).eq.1) goto 429
          delta(m)=-(dlog(rm1(m))+dlog(rm2(m))*dcos(pi2*rm3(m))**2)
     c             /(dxi2(m)*b3)
          if (rm4(m)**2.gt.1.e0-delta(m)/2.e0) goto 429
            a0(m)=1.e0-delta(m)
            jhit(m)=1
 429    continue
      enddo

      call rndpr3(2)

      do m=1,nvect
        if (jhit(m).eq.1)then
        hitp=hitp+1.e0
        c0 =dble (w3(m,2,2)+w3(m,3,3))/dxi2(m)
        c1 =dimag(w3(m,2,3)+w3(m,3,2))/dxi2(m)
        c2 =dble (w3(m,2,3)-w3(m,3,2))/dxi2(m)
        c3 =dimag(w3(m,2,2)-w3(m,3,3))/dxi2(m)
        rad = 1.e0 - a0(m)**2
        a3 = dsqrt(rad)*(2.e0*rm1(m)-1.e0)
        rad2 = dsqrt(dabs(rad-a3**2))
        theta = pi2*rm2(m)
        a1 = rad2*dcos(theta)
        a2 = rad2*dsin(theta)
        d0 = a0(m)*c0+a1*c1+a2*c2+a3*c3
        d1 =-a0(m)*c1+a1*c0+a2*c3-a3*c2
        d2 =-a0(m)*c2+a2*c0+a3*c1-a1*c3
        d3 =-a0(m)*c3+a3*c0+a1*c2-a2*c1
        y22=dcmplx( d0, d3)
        y23=dcmplx( d2, d1)
        y32=dcmplx(-d2, d1)
        y33=dcmplx( d0,-d3)
        v21=v(m,2,1,mu)
        v22=v(m,2,2,mu)
        v23=v(m,2,3,mu)
        v31=v(m,3,1,mu)
        v32=v(m,3,2,mu)
        v33=v(m,3,3,mu)
        v(m,2,1,mu)=y22*v21+y23*v31
        v(m,2,2,mu)=y22*v22+y23*v32
        v(m,2,3,mu)=y22*v23+y23*v33
        v(m,3,1,mu)=y32*v21+y33*v31
        v(m,3,2,mu)=y32*v22+y33*v32
        v(m,3,3,mu)=y32*v23+y33*v33
        endif
      enddo

c - - - - - - - v(1,1)-v(3,3) - - - - - - - - - - - - - -
      do m=1,nvect
        jhit(m)=0
      enddo

      do i=1,3,2
      do j=1,3,2
      do m=1,nvect
        w3(m,i,j)=v(m,i,1,mu)*w(m,1,j)
     c           +v(m,i,2,mu)*w(m,2,j)
     c           +v(m,i,3,mu)*w(m,3,j)
      enddo
      enddo
      enddo

      do m=1,nvect
        ww11=w3(m,1,1)+dconjg(w3(m,3,3))
        ww13=w3(m,1,3)-dconjg(w3(m,3,1))
        ww31=w3(m,3,1)-dconjg(w3(m,1,3))
        ww33=w3(m,3,3)+dconjg(w3(m,1,1))
        dxi2(m)=dsqrt(dble(ww11*ww33-ww13*ww31))
      enddo

      do ihit1=1,nhit
        call rndpr3(4)
        do 529 m = 1,nvect
          if(jhit(m).eq.1) goto 529
          delta(m)=-(dlog(rm1(m))+dlog(rm2(m))*dcos(pi2*rm3(m))**2)
     c             /(dxi2(m)*b3)
          if (rm4(m)**2.gt.1.e0-delta(m)/2.e0) goto 529
            a0(m)=1.e0-delta(m)
            jhit(m)=1
 529    continue
      enddo

      call rndpr3(2)

      do m=1,nvect
        if (jhit(m).eq.1)then
          hitp=hitp+1.e0
          c0 =dble (w3(m,1,1)+w3(m,3,3))/dxi2(m)
          c1 =dimag(w3(m,1,3)+w3(m,3,1))/dxi2(m)
          c2 =dble (w3(m,1,3)-w3(m,3,1))/dxi2(m)
          c3 =dimag(w3(m,1,1)-w3(m,3,3))/dxi2(m)
          rad = 1.e0 - a0(m)**2
          a3 = dsqrt(rad)*(2.e0*rm1(m)-1.e0)
          rad2 = dsqrt(dabs(rad-a3**2))
          theta = pi2*rm2(m)
          a1 = rad2*dcos(theta)
          a2 = rad2*dsin(theta)
          d0 = a0(m)*c0+a1*c1+a2*c2+a3*c3
          d1 =-a0(m)*c1+a1*c0+a2*c3-a3*c2
          d2 =-a0(m)*c2+a2*c0+a3*c1-a1*c3
          d3 =-a0(m)*c3+a3*c0+a1*c2-a2*c1
          y11=dcmplx( d0, d3)
          y13=dcmplx( d2, d1)
          y31=dcmplx(-d2, d1)
          y33=dcmplx( d0,-d3)
          v11=v(m,1,1,mu)
          v12=v(m,1,2,mu)
          v13=v(m,1,3,mu)
          v31=v(m,3,1,mu)
          v32=v(m,3,2,mu)
          v33=v(m,3,3,mu)
          v(m,1,1,mu)=y11*v11+y13*v31
          v(m,1,2,mu)=y11*v12+y13*v32
          v(m,1,3,mu)=y11*v13+y13*v33
          v(m,3,1,mu)=y31*v11+y33*v31
          v(m,3,2,mu)=y31*v12+y33*v32
          v(m,3,3,mu)=y31*v13+y33*v33
        endif
      enddo


      do m=1,nvect
        znrm1=v(m,1,1,mu)*dconjg(v(m,1,1,mu))
     c       +v(m,1,2,mu)*dconjg(v(m,1,2,mu))
     c       +v(m,1,3,mu)*dconjg(v(m,1,3,mu))
*        znrm1=cdsqrt(znrm1)
        znrm1=sqrt(znrm1)
        v(m,1,1,mu)=v(m,1,1,mu)/znrm1
        v(m,1,2,mu)=v(m,1,2,mu)/znrm1
        v(m,1,3,mu)=v(m,1,3,mu)/znrm1
        zipd21=v(m,2,1,mu)*dconjg(v(m,1,1,mu))
     c        +v(m,2,2,mu)*dconjg(v(m,1,2,mu))
     c        +v(m,2,3,mu)*dconjg(v(m,1,3,mu))
        v(m,2,1,mu)=v(m,2,1,mu)-zipd21*v(m,1,1,mu)
        v(m,2,2,mu)=v(m,2,2,mu)-zipd21*v(m,1,2,mu)
        v(m,2,3,mu)=v(m,2,3,mu)-zipd21*v(m,1,3,mu)
        znrm2=v(m,2,1,mu)*dconjg(v(m,2,1,mu))
     c       +v(m,2,2,mu)*dconjg(v(m,2,2,mu))
     c       +v(m,2,3,mu)*dconjg(v(m,2,3,mu))
*        znrm2=cdsqrt(znrm2)
        znrm2=sqrt(znrm2)
        v(m,2,1,mu)=v(m,2,1,mu)/znrm2
        v(m,2,2,mu)=v(m,2,2,mu)/znrm2
        v(m,2,3,mu)=v(m,2,3,mu)/znrm2
        zipd31=v(m,3,1,mu)*dconjg(v(m,1,1,mu))
     c        +v(m,3,2,mu)*dconjg(v(m,1,2,mu))
     c        +v(m,3,3,mu)*dconjg(v(m,1,3,mu))
        zipd32=v(m,3,1,mu)*dconjg(v(m,2,1,mu))
     c        +v(m,3,2,mu)*dconjg(v(m,2,2,mu))
     c        +v(m,3,3,mu)*dconjg(v(m,2,3,mu))
        v(m,3,1,mu)=v(m,3,1,mu)-zipd32*v(m,2,1,mu)-zipd31*v(m,1,1,mu)
        v(m,3,2,mu)=v(m,3,2,mu)-zipd32*v(m,2,2,mu)-zipd31*v(m,1,2,mu)
        v(m,3,3,mu)=v(m,3,3,mu)-zipd32*v(m,2,3,mu)-zipd31*v(m,1,3,mu)
        znrm3=v(m,3,1,mu)*dconjg(v(m,3,1,mu))
     c       +v(m,3,2,mu)*dconjg(v(m,3,2,mu))
     c       +v(m,3,3,mu)*dconjg(v(m,3,3,mu))
*        znrm3=cdsqrt(znrm3)
        znrm3=sqrt(znrm3)
        v(m,3,1,mu)=v(m,3,1,mu)/znrm3
        v(m,3,2,mu)=v(m,3,2,mu)/znrm3
        v(m,3,3,mu)=v(m,3,3,mu)/znrm3
      enddo

      do i=1,3
      do m=1,nvect
        ztr=ztr+w(m,i,1)*v(m,1,i,mu)
     c         +w(m,i,2)*v(m,2,i,mu)
     c         +w(m,i,3)*v(m,3,i,mu)
      enddo
      enddo

 320  continue


 1000 continue

      trace=dreal(ztr)/dble(3*6*nlink)
      hitp =hitp/dble(3*nlink*iter)

      return
      end
c*********************************************************************c
c*********************************************************************c
      subroutine dir
c
c     * for SU(2) and SU(3) Lattice on VP.
c     ------------------------------------
c
c     * Compute list vectors under periodic boundary condition.
c       For the measurement of observables, the list vectors computed
c       here are on the ordinary lattice.
c
c
c     ........       nvect           ........       nsite
c     _     _     _
c     7  7  8  8  9  9               13 14 15 16 17 18
c        _     _     _
c     4  4  5  5  6  6      ===>>    7  8  9  10 11 12
c     _     _     _
c     1  1  2  2  3  3               1  2  3  4  5  6
c
c     even-odd site numbers.         ordinary site numbers.
c     Bar{n} means odd site.
c
c
c     * List vectors Mioe and Mioo are used in order to convert
c       arrays on even-odd lattices into those on ordinary lattices.
c         mioe(n) = 1   (n means odd site)
c                 = 0   (n means even site), n: ordinary site number.
c
c         mioo(n) = even-odd site number   , n: ordinary site number.
c
c     * Im and nm give the neighboring site numbers on ordinary
c       lattice.
c         im(n,mu,i) = # of (n + (i*mu))   , n: ordinary site number
c         nm(n,mu)   = # of (n - mu)       , n: ordinary site number
c
c     * Include file
c         paravp3
c
c     * Output arrays and variables
c         /tb1/
c         /sdb2/
c
c     * Work arrays
c         /sdb1/
c
c     * programmed by KEK group.
c     * modified by S.Kitahara. (96.8.2)
c----------------------------------------------------------------------

      include'paravp3'

      common /tb1  /mioe(nsite), mioo(nsite)
      common /sdb1 /ipower(nd),isize(nd),is(nd),mio(nvect),mie(nvect)
      common /sdb2 /im(nsite,nd,n0), nm(nsite,nd)

      isize(1)=n1
      isize(2)=n2
      isize(3)=n3
      isize(4)=n4
      ipower(1)= 1
      ipower(2)= n1
      ipower(3)= n1*n2
      ipower(4)= n1*n2*n3

      do 1 i4=1,n4
        is(4)= i4
        m4   = (i4-1)*ipower(4)
        mm4  = (i4-1)*n0*n2*n3
      do 1 i3=1,n3
        is(3)= i3
        m3   = m4 +(i3-1)*ipower(3)
        mm3  = mm4+(i3-1)*n0*n2
      do 1 i2=1,n2
        is(2)= i2
        m2   = m3 +(i2-1)*ipower(2)
        mm2  = mm3+(i2-1)*n0
        ioe  = mod(i2+i3+i4,2)
      do 1 i1=1,n0
        is(1)= i1
        m1   = m2 +2*i1
        mm1  = mm2+i1
        mio(mm1) = m1-ioe
        mie(mm1) = m1+ioe-1
    1 continue

      do 2 m=1,nvect
        mioe(mio(m)) = 1
        mioo(mio(m)) = m
        mioe(mie(m)) = 0
        mioo(mie(m)) = m
    2 continue

      do 3 i4=1,n4
        is(4)= i4
        m4   = (i4-1)*ipower(4)
      do 3 i3=1,n3
        is(3)= i3
        m3   = m4 +(i3-1)*ipower(3)
      do 3 i2=1,n2
        is(2)= i2
        m2   = m3 +(i2-1)*ipower(2)
      do 3 i1=1,n1
        is(1)= i1
        m1   = m2 +i1
      do 3 l=1,nd
        nm(m1,l)=m1+(mod(is(l)-2+isize(l),isize(l))+1-is(l))*ipower(l)
      do 3 k=1,n0
        im(m1,l,k) = m1+(mod(is(l)+k-1,isize(l))+1-is(l))*ipower(l)
    3 continue
      return
      end
c*********************************************************************c
      subroutine dirct
c
c     * for SU(2) and SU(3) Lattice on VP.
c     ------------------------------------
c
c     * Compute list vectors under periodic boundary condition.
c       For vector processing, the entire lattice is divided into
c       two parts: even lattice sites and odd lattice sites.
c
c     * The list vectors are mainly used by Monte and Mafix1.
c
c
c     ........       nvect
c     _     _     _
c     7  7  8  8  9  9
c        _     _     _
c     4  4  5  5  6  6
c     _     _     _
c     1  1  2  2  3  3
c
c     even-odd site numbers.
c     Bar{n} means odd site.
c
c
c     * mo1,mo2,me1 and me2 give the neighboring site numbers.
c         mo1(n,mu) = # of (n + mu)   (n: odd,  n+mu: even)
c         mo2(n,mu) = # of (n - mu)   (n: odd,  n-mu: even)
c         me1(n,mu) = # of (n + mu)   (n: even, n+mu: odd )
c         me2(n,mu) = # of (n - mu)   (n: even, n-mu: odd )
c
c     * ldir2(i,mu) (i=1,2 and 3) gives a direction perpendicular to mu.
c       ldir2(i,1),ldir2(i,2),ldir2(i,3) and ldir2(i,4) are
c       perpendicular to each other.
c
c     * Include file
c         paravp3
c
c     * Output arrays and variables
c         /tabl/
c
c     * Work arrays
c         /vard/
c
c     * programmed by KEK group.
c     * modified by S.Kitahara. (96.8.2)
c----------------------------------------------------------------------

      include'paravp3'

      common /tabl /mo1(nvect,nd), mo2(nvect,nd),
     &              me1(nvect,nd), me2(nvect,nd),
     &              ldir2(nd-1,nd)
      common /vard /mcycle(0:n0+1), mup(50,2:nd), mdn(50,2:nd),
     &              ipower(2:nd),   is(2:nd),     isize(nd)

      do i =1,3
      do mu=1,nd
        ldir2(i,mu) = mod(i+mu-1,nd) + 1
      enddo
      enddo

      isize(1)   =  n1
      isize(2)   =  n2
      isize(3)   =  n3
      isize(4)   =  n4

      ipower(2)  =  n0
      ipower(3)  =  n0*n2
      ipower(4)  =  n0*n2*n3

      do 1 m  = 0,n0+1
    1    mcycle(m) = mod(m-1+n0,n0)+1

      do 2 i = 2,nd
      do 2 m = 1,50
         mup(m,i)  = mod(m,isize(i))+1
    2    mdn(m,i)  = mod(m-2+isize(i),isize(i))+1
      do 3 i4  = 1,n4
         is(4) = i4
         m4    = (i4-1)*ipower(4)
      do 3 i3  = 1,n3
         is(3) = i3
         m3    = m4+(i3-1)*ipower(3)
      do 3 i2  = 1,n2
         is(2) = i2
         m2    = m3+(i2-1)*ipower(2)
         ioe   = mod(i2+i3+i4,2)
      do 3 i1  = 1,n0
         m1    = m2+i1
         mo1(m1,1) = m1 + mcycle(i1+1-ioe)-i1
         mo2(m1,1) = m1 + mcycle(i1  -ioe)-i1
         me1(m1,1) = m1 + mcycle(i1  +ioe)-i1
         me2(m1,1) = m1 + mcycle(i1-1+ioe)-i1
       do 11 i = 2,nd
           mo1(m1,i) = m1 + (mup(is(i),i)-is(i))*ipower(i)
           me1(m1,i) = mo1(m1,i)
           mo2(m1,i) = m1 + (mdn(is(i),i)-is(i))*ipower(i)
           me2(m1,i) = mo2(m1,i)
   11  continue
    3  continue
       return
       end
c*********************************************************************c
c*********************************************************************c
      subroutine inits(ints,pi2)
c
c     * for SU(3) Lattice on VP.
c     --------------------------
c
c     * This routine returns an initial configuration.
c
c     * This routine calls
c         Rndpr3.
c
c     * Include file
c         paravp3
c
c     * Input parameters and constants
c         ints   : 1 = cold start. 0 = hot start.
c         pi2    : 2 times pi (6.28318...).
c
c     * Output arrays and variables
c         /var2/ : generated initial configuration.
c
c     * Work arrays
c         /ran1/ : random numbers generated by Rndpr3.
c
c     * programmed by S.Kitahara.
c     * modified by S.Kitahara. (98.10.14)
c----------------------------------------------------------------------

      include'paravp3'

      common /var2/ u(nvect,3,3,nd),v(nvect,3,3,nd)
      common /ran1/ rm1(nvect),rm2(nvect),rm3(nvect),rm4(nvect)

      if(ints.eq.0)then
        do i = 1,3
        do j = 1,3
        do mu = 1,nd
        do m = 1,nvect
          u(m,i,j,mu) = (0.e0,0.e0)
          v(m,i,j,mu) = (0.e0,0.e0)
        enddo
        enddo
        enddo
        enddo

        do i = 1,3
        do mu = 1,nd
        do m = 1,nvect
          u(m,i,i,mu) = (1.e0,0.e0)
          v(m,i,i,mu) = (1.e0,0.e0)
        enddo
        enddo
        enddo

      elseif(ints.eq.1)then

        do mu=1,nd
          call rndpr3(4)
          do m=1,nvect
            a0  =rm1(m)
            rad = 1.e0 - a0**2
            a3  = dsqrt(rad)*(2.e0*rm2(m)-1.e0)
            rad2 = dsqrt(dabs(rad-a3**2))
            theta= pi2*rm3(m)
            a1 = rad2*dcos(theta)
            a2 = rad2*dsin(theta)
            y11=dcmplx( a0, a3)
            y12=dcmplx( a2, a1)
            y21=dcmplx(-a2, a1)
            y22=dcmplx( a0,-a3)
            u(m,1,1,mu)=y11
            u(m,1,2,mu)=y12
            u(m,2,1,mu)=y21
            u(m,2,2,mu)=y22
          enddo
        enddo

        do mu=1,nd
          call rndpr3(4)
          do m=1,nvect
            a0  =rm1(m)
            rad = 1.e0 - a0**2
            a3  = dsqrt(rad)*(2.e0*rm2(m)-1.e0)
            rad2 = dsqrt(dabs(rad-a3**2))
            theta= pi2*rm3(m)
            a1 = rad2*dcos(theta)
            a2 = rad2*dsin(theta)
            y22=dcmplx( a0, a3)
            y23=dcmplx( a2, a1)
            y32=dcmplx(-a2, a1)
            y33=dcmplx( a0,-a3)
            u12=u(m,1,2,mu)
            u22=u(m,2,2,mu)
            u(m,1,2,mu)=u12*y22
            u(m,1,3,mu)=u12*y23
            u(m,2,2,mu)=u22*y22
            u(m,2,3,mu)=u22*y23
            u(m,3,2,mu)=y32
            u(m,3,3,mu)=y33
          enddo
        enddo

        do mu=1,nd
          call rndpr3(4)
          do m=1,nvect
            a0  =rm1(m)
            rad = 1.e0 - a0**2
            a3  = dsqrt(rad)*(2.e0*rm2(m)-1.e0)
            rad2 = dsqrt(dabs(rad-a3**2))
            theta= pi2*rm3(m)
            a1 = rad2*dcos(theta)
            a2 = rad2*dsin(theta)
            y11=dcmplx( a0, a3)
            y13=dcmplx( a2, a1)
            y31=dcmplx(-a2, a1)
            y33=dcmplx( a0,-a3)
            u11=u(m,1,1,mu)
            u13=u(m,1,3,mu)
            u21=u(m,2,1,mu)
            u23=u(m,2,3,mu)
            u33=u(m,3,3,mu)
            u(m,1,1,mu)=u11*y11+u13*y31
            u(m,1,3,mu)=u11*y13+u13*y33
            u(m,2,1,mu)=u21*y11+u23*y31
            u(m,2,3,mu)=u21*y13+u23*y33
            u(m,3,1,mu)=u33*y31
            u(m,3,3,mu)=u33*y33
          enddo
        enddo

        do mu=1,nd
          call rndpr3(4)
          do m=1,nvect
            a0  =rm1(m)
            rad = 1.e0 - a0**2
            a3  = dsqrt(rad)*(2.e0*rm2(m)-1.e0)
            rad2 = dsqrt(dabs(rad-a3**2))
            theta= pi2*rm3(m)
            a1 = rad2*dcos(theta)
            a2 = rad2*dsin(theta)
            y11=dcmplx( a0, a3)
            y12=dcmplx( a2, a1)
            y21=dcmplx(-a2, a1)
            y22=dcmplx( a0,-a3)
            v(m,1,1,mu)=y11
            v(m,1,2,mu)=y12
            v(m,2,1,mu)=y21
            v(m,2,2,mu)=y22
          enddo
        enddo

        do mu=1,nd
          call rndpr3(4)
          do m=1,nvect
            a0  =rm1(m)
            rad = 1.e0 - a0**2
            a3  = dsqrt(rad)*(2.e0*rm2(m)-1.e0)
            rad2 = dsqrt(dabs(rad-a3**2))
            theta= pi2*rm3(m)
            a1 = rad2*dcos(theta)
            a2 = rad2*dsin(theta)
            y22=dcmplx( a0, a3)
            y23=dcmplx( a2, a1)
            y32=dcmplx(-a2, a1)
            y33=dcmplx( a0,-a3)
            v12=v(m,1,2,mu)
            v22=v(m,2,2,mu)
            v(m,1,2,mu)=v12*y22
            v(m,1,3,mu)=v12*y23
            v(m,2,2,mu)=v22*y22
            v(m,2,3,mu)=v22*y23
            v(m,3,2,mu)=y32
            v(m,3,3,mu)=y33
          enddo
        enddo

        do mu=1,nd
          call rndpr3(4)
          do m=1,nvect
            a0  =rm1(m)
            rad = 1.e0 - a0**2
            a3  = dsqrt(rad)*(2.e0*rm2(m)-1.e0)
            rad2 = dsqrt(dabs(rad-a3**2))
            theta= pi2*rm3(m)
            a1 = rad2*dcos(theta)
            a2 = rad2*dsin(theta)
            y11=dcmplx( a0, a3)
            y13=dcmplx( a2, a1)
            y31=dcmplx(-a2, a1)
            y33=dcmplx( a0,-a3)
            v11=v(m,1,1,mu)
            v13=v(m,1,3,mu)
            v21=v(m,2,1,mu)
            v23=v(m,2,3,mu)
            v33=v(m,3,3,mu)
            v(m,1,1,mu)=v11*y11+v13*y31
            v(m,1,3,mu)=v11*y13+v13*y33
            v(m,2,1,mu)=v21*y11+v23*y31
            v(m,2,3,mu)=v21*y13+v23*y33
            v(m,3,1,mu)=v33*y31
            v(m,3,3,mu)=v33*y33
          enddo
        enddo
      endif

      return
      end
c*********************************************************************c
      subroutine rndpr3(nrn)
c
c     * for SU(2) and SU(3) Lattice on VP.
c     ------------------------------------
c
c     * random number generator ( m-series,p=250,q=103 )
c
c     * Before this routine is called, you must once call
c         Seed.
c
c     * Include file
c         paravp3
c
c     * Input arrays and parameters
c         nrn    : must be 2 or 4.
c                  Rm1 and Rm2 are changed(nrn = 2).
c                  Rm1, Rm2, Rm3 and Rm4 are changed(nrn = 4).
c         /rndk/ : Initial value of K1 and K2 are 0 and 147
c                  respectively which you must give elsewhere.
c         /rndlp/: The array, Irand is computed by Seed.
c
c     * Output arrays
c         /ran1/
c
c     * Work arrays
c         /rand/
c
c     * programmed by J.Makino
c     * modified   by S.Hioki
c     * modified by S.Kitahara. (96.8.2)
c----------------------------------------------------------------------

      include'paravp3'

      parameter( nkei=256 )

      common /rndlp/irand(nkei,250)
      common /rand /rn(nsite*2+nkei)
      common /rndk /k1,k2
      common /ran1 /rm1(nvect),rm2(nvect),rm3(nvect),rm4(nvect)

      ic = 2*(2 ** 30 - 1) + 1
      rc = 1.0/ic

c --  main production

      mkei=(nvect*nrn-1)/nkei
      do 10  j0 = 0 , mkei
        j1=nkei*j0
        k1=mod(k1,250)+1
        k2=mod(k2,250)+1
      do 20 j = 1 , nkei
        irr           = ieor( irand(j,k1),irand(j,k2) )
        rn(j+j1)      = dble( irr )*rc
        irand(j,k1)= irr
   20 continue
   10 continue

      if(nrn.eq.2) then
        do m=1,nvect
          rm1(m)=rn(m)
          rm2(m)=rn(m+nvect)
        enddo
      else
        do m=1,nvect
          rm1(m)=rn(m)
          rm2(m)=rn(m+nvect)
          rm3(m)=rn(m+nvect2)
          rm4(m)=rn(m+nvect3)
        enddo
      endif

      return
      end
c*********************************************************************c
      subroutine seed
c
c     * for SU(2) and SU(3) Lattice on VP.
c     ------------------------------------
c
c     * random number seed preparation.
c
c     * Output arrays
c         /rndlp/
c
c     * programmed by J.Makino
c----------------------------------------------------------------------

      parameter( nkei=256 )
      common /rndlp/irand(nkei,250)
      dimension irand0(500),l(0:499)
      data ix/1774315169/
      do 110 i1=1,250
      do 110 i2=1,nkei
        irand(i2,i1)=0
  110 continue
      do 10 i=1,250
        ix=iand(ix*48828125,2147483647)
        irand0(i)=ix
   10 continue
      do 20 i=251,500
        irand0(i)=ieor(irand0(i-250),irand0(i-103))
   20 continue
      do 900 k=1,nkei
        do 100 i=0,499
          l(i)=0
  100   continue
        if(k.le.249) then
          l(k)=1
        else
          l(k-250)=1
          l(k-103)=1
        endif
        do 500 j=242,1,-1
          do 200 i=249,0,-1
            l(2*i+1)=0
            l(2*i)=l(i)
  200     continue
          do 300 i=498,250,-1
            l(i-250)=ieor(l(i-250),l(i))
            l(i-103)=ieor(l(i-103),l(i))
  300     continue
  500   continue
        do 700 j=1,250
          do 600 i=0,249
            irand(k,j)=ieor(irand(k,j),l(i)*irand0(i+j))
  600     continue
  700   continue
  900 continue
      return
      end
