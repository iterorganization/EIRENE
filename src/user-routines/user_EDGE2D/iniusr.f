      SUBROUTINE EIRENE_INIUSR
      USE EIRMOD_CINIT
      USE EIRMOD_PARMMOD, ONLY: IFOFF, NLIM, NSTS
      USE EIRMOD_COMUSR,  ONLY: NPLSI, NPLS_FIX
      USE EIRMOD_USRDATA
      IMPLICIT NONE

      character(256) :: filename, line
      integer :: ll,nlines,ind,itec1,itec2,j,ier
      integer, parameter :: fp=9998

      integer :: idum1,idum2,idum3,idum4
      logical :: lex

c     set default namelist values
      eirene_ftransfer = 'eirene.transfer'
      eirene_njetto=0
cdmh
      eirene_fstoreneutflux = 'eirene.chemFluxDep'
      eirene_wallFluxModel = 1
cdmh
c     get namelist config
      open(unit=fp,file='eirene_user.namelist')
      read(fp,eirene_user)
      close(fp)
      ll=len_trim(casename)

! Read triangle map file
      filename=casename(1:ll) // '.trimap'
      open (unit=fp,file=filename)
      read(fp,'(3i8)',iostat=ier) nlines, nrpla
      close(fp)

! Read surfaces file to identify targets
      filename=casename(1:ll) // '.surfaces'
      open (unit=fp,file=filename)
      ! count lines
      nlines = 0
      read(fp,'(a)',iostat=ier) line ! skip headers
      do while (ier==0)
        read(fp,'(a)',iostat=ier) line
        nlines = nlines + 1
      enddo
      call eirene_alloc_usrdata(max(NLIM+NSTS,nlines))

      ! read surfaces
      rewind(fp)
      read(fp,'(a)',iostat=ier) line ! skip headers
      j = 0
      do while (ier==0)
        read(fp,'(3i7,a)',iostat=ier) ind, itec1, itec2, line
        if(itec1<=PUMP.or.itec1>WALL) cycle
        j = j+1
        surftype(j) = itec1
      enddo
      close(fp)

!     set npls_fix to handle NN collisions
      ll=len_trim(casename)
      filename=casename(1:ll) // '.species'
      inquire(file=filename, exist=lex)
      if(lex) then
        open (unit=fp,file=filename,access='sequential',
     &        form='formatted')
        read(fp,*) line
        read(fp,*) idum1,idum2,idum3,idum4,nbkgusr
        close(fp)
      endif

      npls_fix = nplsi
      if(nbkgusr>0) npls_fix = nplsi - nbkgusr

      RETURN
      END SUBROUTINE EIRENE_INIUSR
