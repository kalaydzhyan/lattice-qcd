      subroutine start_log_arnoldi()
      include 'debug.h'



      logical              :: exist


      ndigit = -3
      logfil = 5
      mcaupd = 1
      mnaupd = 1
      mcaup2 = 2
      mnaup2 = 2
      mcaitr = 0
      mnaitr = 0
      mceupd = 1
      mneupd = 1


      inquire(file='arnoldi.log', exist=exist)

      if (exist) then
        open(unit=5, access='append', err=10, file='arnoldi.log', 
     c        form='formatted', status='old')

      else
         open(unit=5, access='sequential', err=10, file='arnoldi.log',
     c        form='formatted', status='new')
      end if


      return

 10   print *,'error opening file'

      end

      subroutine end_log_arnoldi()


      close(5)

      return



      end
