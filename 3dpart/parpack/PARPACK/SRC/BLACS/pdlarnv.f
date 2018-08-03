c\BeginDoc
c
c\Name: pdlarnv
c
c Message Passing Layer: BLACS
c
c\Description:
c
c  Parallel Version of ARPACK utility routine dlarnv
c
c  PSLARNV returns a vector of n (nloc) random real numbers from a uniform or
c  normal distribution. It is assumed that X is distributed across a 1-D array 
c  of processors ( nprocs < 1000 )
c
c\Arguments
c  COMM    BLACS Communicator for the processor grid
c
c  IDIST   (input) INTEGER
c          Specifies the distribution of the random numbers:
c          = 1:  uniform (0,1)
c          = 2:  uniform (-1,1)
c          = 3:  normal (0,1)
c
c  ISEED   (input/output) INTEGER array, dimension (4)
c          On entry, the seed of the random number generator; the array
c          elements must be between 0 and 4095, and ISEED(4) must be
c          odd.
c          On exit, the seed is updated.
c
c  N       (input) INTEGER
c          The number of random numbers to be generated.
c
c  X       (output) Double precision array, dimension (N)
c          The generated random numbers.
c
c\Author: Kristi Maschhoff
c
c\Details
c
c  Simple parallel version of LAPACK auxiliary routine dlarnv 
c  for X distributed across a 1-D array of processors.
c  This routine calls the auxiliary routine SLARNV to generate random
c  real numbers from a uniform (0,1) distribution. Output is consistent
c  with serial version. 
c
c\SCCS Information: 
c FILE: larnv.F   SID: 1.1   DATE OF SID: 1/23/96   
c
c-----------------------------------------------------------------------
c
      subroutine pdlarnv( comm, idist, iseed, n, x )
c
c     .. BLACS VARIABLES AND FUNCTIONS ..
      integer    comm, nprow, npcol, myprow, mypcol
c
c     .. External Functions ..
      external   BLACS_GRIDINFO
c     ..
c     .. Scalar Arguments ..
      integer			i, n
c     ..
c     .. Array Arguments ..
      integer			iseed( 4 )
      Double precision			
     &                  x( * )
c     ..
c     .. Local Array Arguments ..
      integer			dist(1000)
c     ..
c     .. External Subroutines ..
      external			dlarnv
c     ..
c     .. Executable Statements ..
c
      call BLACS_GRIDINFO( comm, nprow, npcol, myprow, mypcol )
c
      do 10 i=1,nprow
         dist(i) = 0
 10   continue
c
      dist( myprow + 1 ) = n
      call IGSUM2D( comm, 'ALL', ' ', nprow, 1, dist, nprow, -1, -1 )
c
      do 15 i = 1, myprow + 1
         call dlarnv ( idist, iseed, dist(i), x )
 15   continue
c
      return
      end
