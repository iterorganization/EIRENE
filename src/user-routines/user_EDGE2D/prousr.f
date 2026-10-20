CDK USER
C
C   USER-SUPPLIED SUBROUTINES
C
C           ************
C           *  EDGE2D  *  (fem-interface)
C           ************
C
C s.wiesen@fz-juelich.de (mar06)
C
C get plasma values for bulk species IPLS (module comprt, set by subr. plasma)
c
c reads casename.zplasma
c
      SUBROUTINE EIRENE_PROUSR (PRO,INDX,P0,P1,P2,P3,P4,P5,PROVAC,N)
C
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CGRID
      USE EIRMOD_CGEOM
      USE EIRMOD_CINIT
      USE EIRMOD_CCONA
      USE EIRMOD_COMPRT
      USE EIRMOD_USRDATA

      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: P0, P1, P2, P3, P4, P5, PROVAC
      REAL(DP), INTENT(OUT) :: PRO(*)
      INTEGER, INTENT(IN) :: INDX, N
      CHARACTER(256) :: FILENAME, line,sstr
      INTEGER :: NTR, I, J, LL,ier
      REAL(DP), ALLOCATABLE, SAVE :: PLAS(:,:,:)
      REAL(DP) :: bnorm
      INTEGER, SAVE :: INDAR(11)=(/ (0, i=1,11) /)
      integer, parameter :: fp=31
      character(2) :: cstr2
      integer :: idum,idum2,idum3,idum4,idum5
      real(dp) :: rdum
      INTEGER, SAVE :: IPLSTI, IPLSV
csw
c     reset?
      if(indx < 0) then
         if(allocated(plas)) deallocate(plas)
         indar(1:11) = 0
c just in case...:
         close(fp+ifoff)
         return
      endif
csw
       npls_fix = nplsi
       if(nbkgusr>0) npls_fix = nplsi - nbkgusr

c read in plasma data from fort.31 ?
      if (.not.allocated(plas)) then
        allocate(plas(12,nrad,0:npls))
        plas=0.

        ll=len_trim(casename)
        filename=casename(1:ll) // '.zplasma'
        open (unit=fp+ifoff,file=filename,access='sequential',
     .        form='formatted')


c misc plasma data:
        write(sstr,'(a20)')
     .          '*** MISC PLASMA DATA'
        CALL EIRENE_locstr_usr(fp+ifoff,sstr,ier)
        if(ier /=0) then
           write(*,*) 'PROUSR: ',sstr,' not found'
           close(fp+ifoff)
           call EIRENE_exit_own(1)
        endif
        do j=1,2
           read(fp+ifoff,'(a)') line
        enddo

        read(fp+ifoff,'(i8)') ntr
        if (ntr /= nr1st-1) then
           write (*,*) 'PROUSR:', sstr
           write (*,*) ' wrong number of triangles in plasma file'
           write (*,*) ' check for correct number in file ',filename,
     .                      ntr,nr1st-1
           call EIRENE_exit_own(1)
        endif

        do j=1,ntr
           read(fp+ifoff,'(i7,7(1x,e14.7))') idum,
c     .          te,ne,bx,by,bz,pot,psi
     .          plas(1,j,0),
     .          plas(12,j,0),
     .          plas(7,j,0),
     .          plas(8,j,0),
     .          plas(9,j,0),
     .          rdum,
!pb     .          rdum
     .          plas(6,j,0)

           bnorm=sqrt(plas(7,j,0)**2 + plas(8,j,0)**2 + plas(9,j,0)**2)
           if (bnorm < eps12) then
              plas(7,j,0) = 0._dp
              plas(8,j,0) = 0._dp
              plas(9,j,0) = 1._dp
              plas(10,j,0) = 0._dp
           else
              plas(7:9,j,0) = plas(7:9,j,0)/bnorm
              plas(10,j,0) =  bnorm
           endif
        enddo
        plas(7,ntr+1:nrad,0) = 0._dp
        plas(8,ntr+1:nrad,0) = 0._dp
        plas(9,ntr+1:nrad,0) = 1._dp
        plas(10,ntr+1:nrad,0) = 1._dp

c loop over species:
        do i=1,npls_fix
           write(cstr2,'(i2.2)') i
           write(sstr,'(a23)')
     .          '*** ION #'//cstr2//' PLASMA DATA'
           CALL EIRENE_locstr_usr(fp+ifoff,sstr,ier)
           if(ier /=0) then
              write(*,*) 'PROUSR: tag ',sstr,'not found'
              close(fp+ifoff)
              call EIRENE_exit_own(1)
           endif

           do j=1,3
              read(fp+ifoff,'(a)') line
           enddo

           read(fp+ifoff,'(i8)') ntr
           if (ntr /= nr1st-1) then
              write (*,*) 'PROUSR:', sstr
              write (*,*) ' wrong number of triangles in plasma file'
              write (*,*) ' check for correct number in file ',filename,
     .                      ntr,nr1st-1
              call EIRENE_exit_own(1)
           endif

           do j=1,ntr
              read(fp+ifoff,'(i7,6(1x,e14.7))') idum,
c     .             ti,ni,vx,vy,vz,zi
     .             plas(2,j,i),
     .             plas(3,j,i),
     .             plas(4,j,i),
     .             plas(5,j,i),
     .             plas(6,j,i),
     .             plas(11,j,i)
           enddo
        enddo

        close(fp+ifoff)
      endif

c
c fort.31 read in now
c
c sort in values
c
      if(indx == 0) then
c te
         pro(1:n) = plas(1,1:n,0)
c ne (fill using global). Temporary hack, better with another indx
         if(eirene_use_e2d_ne) dein(1:n) = plas(12,1:n,0)
         indar(1) = indar(1)+1
c     WARNING!: Assume te is always called before Ti an Vi
         IPLSTI = 0
         IPLSV  = 0
      elseif (indx == 1) then
c ti
         IPLSTI = IPLSTI+1
         IF (IPLSTI .GT. NPLSI) THEN
            write(IUNOUT,*) "PROUSR CALLED FOR TI WITH IPLSTI > NPLS"
            write(IUNOUT,*) "IPLSTI = ",IPLSTI," NPLS=",NPLS
            call EIRENE_exit_own(1)
         ENDIF
         IF(IPLSTI .GT. npls_fix) then
             pro(1:n) = tiin(iplsti,1:n)
         ELSE
             pro(1:n) = plas(2,1:n,iplsti)
         ENDIF
         indar(2) = indar(2)+1

      elseif (indx == 1+1*npls) then
c ni
         IF(ipls .GT. npls_fix) then
             pro(1:n) = diin(ipls,1:n)
         ELSE
             pro(1:n) = plas(3,1:n,ipls)
         ENDIF
         indar(3) = indar(3) + 1

      elseif (indx == 1+2*npls) then
c vx
         IPLSV = IPLSV+1
         IF (IPLSV .GT. NPLSI) THEN
            write(IUNOUT,*) "PROUSR CALLED FOR VX WITH IPLSV > NPLS"
            write(IUNOUT,*) "IPLSV = ",IPLSV," NPLS=",NPLS
            call EIRENE_exit_own(1)
         ENDIF
         IF(ipls .GT. npls_fix) then
             pro(1:n) = vxin(ipls,1:n)
         ELSE
             pro(1:n) = plas(4,1:n,iplsv)
         ENDIF
         indar(4) = indar(4) + 1

      elseif (indx == 1+3*npls) then
c vy
         IF(ipls .GT. npls_fix) then
             pro(1:n) = vyin(ipls,1:n)
         ELSE
             pro(1:n) = plas(5,1:n,iplsv)
         ENDIF
         indar(5) = indar(5) + 1

      elseif (indx == 1+4*npls) then
c vz
         IF(ipls .GT. npls_fix) then
             pro(1:n) = vzin(ipls,1:n)
         ELSE
             pro(1:n) = plas(6,1:n,iplsv)
         ENDIF
         indar(6) = indar(6) + 1

      elseif (indx == 1+1*npls+NPLSTI+3*NPLSV) then
c bx
         pro(1:n) = plas(7,1:n,0)
         indar(7) = indar(7) + 1

      elseif (indx == 2+1*npls+NPLSTI+3*NPLSV) then
c by
         pro(1:n) = plas(8,1:n,0)
         indar(8) = indar(8) + 1

      elseif (indx == 3+1*npls+NPLSTI+3*NPLSV) then
c bz
         pro(1:n) = plas(9,1:n,0)
         indar(9) = indar(9) + 1

      elseif (indx == 4+1*NPLS+NPLSTI+3*NPLSV) then
c bf
         pro(1:n) = plas(10,1:n,0)
         indar(10) = indar(10) + 1

c      reserved for electric field:
c      EXIN: indx=7+1*NPLS+NPLSTI+3*NPLSV
c      EYIN: indx=8+1*NPLS+NPLSTI+3*NPLSV
c      EZIN: indx=9+1*NPLS+NPLSTI+3*NPLSV
c      EFIN: indx=10+1*NPLS+NPLSTI+3*NPLSV

      elseif (indx == 11+1*NPLS+NPLSTI+3*NPLSV) then
! psi
         pro(1:n) = plas(6,1:n,0)

c      ZI for charge bundling
      elseif (indx == 12+1*NPLS+NPLSTI+3*NPLSV) then
c! zi
         IF(ipls .GT. npls_fix) then
             pro(1:n) = ziin(ipls,1:n)
         ELSE
             pro(1:n) = plas(11,1:n,ipls)
         ENDIF

      else
         write (iunout,*) ' prousr: no data provided for index ',indx
         pro(1:n) = 0._dp
      endif

      return
      end subroutine eirene_prousr
