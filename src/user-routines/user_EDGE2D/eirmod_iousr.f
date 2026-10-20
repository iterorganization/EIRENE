      module eirmod_iousr

C  User specific read routines

      USE EIRMOD_PRECISION
      use eirmod_parmmod
      IMPLICIT NONE

      private

      public :: eirene_read_block_11_usr, eirene_write_block_11_usr,
     .          eirene_read_block_11_usr_fixed

      interface eirene_read_block_11_usr
        procedure :: eirene_read_block_11_usr_fixed
#ifndef NO_JSON
        procedure :: eirene_read_block_11_usr_json
#endif
      end interface eirene_read_block_11_usr

#ifndef NO_JSON
      interface eirene_write_block_11_usr
        procedure :: eirene_read_block_11_usr_json
      end interface eirene_write_block_11_usr
#endif

      contains

      subroutine eirene_read_block_11_usr_fixed

      use eirmod_comprt, only: iunin, iunout

      implicit none
      integer :: j, idumm(9)

      return
      end subroutine eirene_read_block_11_usr_fixed

#ifndef NO_JSON
      subroutine eirene_read_block_11_usr_json(json, me)

      use eirmod_json

      integer :: json
      integer :: me

      return
      end subroutine eirene_read_block_11_usr_json


      subroutine eirene_write_block_11_usr_json(json, me)

      use eirmod_json

      integer :: json
      integer :: me

      return
      end subroutine eirene_write_block_11_usr_json
#else
      subroutine eirene_write_block_11_usr
      return
      end subroutine eirene_write_block_11_usr
#endif

      end module eirmod_iousr


