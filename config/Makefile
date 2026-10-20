# program name
EIRLIB = libeirene.a
PROG   = eirobj
PROGX  = eirobjx

# Test whether necessary environment variables are defined; if not, exit
ifndef HOST_NAME
  $(error HOST_NAME not defined)
endif
ifndef COMPILER
  $(error COMPILER not defined)
endif

MAKETAGS ?= ctags -e -f

SRCDIR=${PWD}/src
SRCEIR=${PWD}
SRCLOCAL=${SRCEIR}/src.local
ifdef USE_B25
B2SRC=${SOLPSTOP}/modules/B2.5/src
endif
ifndef INTDIR
 ifdef SOLPSTOP
  SOLPS_MINOR_VERSION = $(shell ( cd ${SOLPSTOP} ; git describe --tags --dirty --always | cut -d '.' -f 2 ) )
  ifeq ($(shell test ${SOLPS_MINOR_VERSION} -ge 2; echo $$?),0)
   INTDIR ?= ${SRCDIR}/interfaces/couple_SOLPS_WG
   $(warning Assuming SOLPS_WG interface.)
  else
   INTDIR ?= ${SRCDIR}/interfaces/couple_SOLPS-ITER
   $(warning Assuming SOLPS-ITER interface.)
  endif
 else
  $(warning Need to choose an interface directory!)
  $(warning Set it with the command:)
  $(warning setenv INTDIR ${SRCDIR}/interfaces/couple_[Name-of-interface])
  $(error Interface directory not defined!)
 endif
else
 $(warning Using ${INTDIR} for the interface directory.)
endif

MAKES = ${SRCEIR}/config/Makefile ${SRCEIR}/config/compile ${SRCEIR}/config/config.${HOST_NAME}.${COMPILER}
DEFINES = ${EIRENE_DEFINES} ${SOLPS_CPP}
# Include global SOLPS compiler settings
ifndef SOLPS_CPP
  NODENAME = $(shell echo `hostname`)
  ifeq ($(shell [ -e ${SOLPSTOP}/SETUP/config.${HOST_NAME}.${COMPILER} ] && echo yes || echo no ),yes)
    include ${SOLPSTOP}/SETUP/config.${HOST_NAME}.${COMPILER}
    MAKES += ${SOLPSTOP}/SETUP/setup.csh.${HOST_NAME}.${COMPILER} ${SOLPSTOP}/SETUP/config.${HOST_NAME}.${COMPILER}
  else
    $(warning ${SOLPSTOP}/SETUP/config.${HOST_NAME}.${COMPILER} not found. Assuming stand-alone compilation.)
  endif
  ifeq ($(shell [ -e ${SOLPSTOP}/SETUP/config.common.${COMPILER} ] && echo yes || echo no ),yes)
    include ${SOLPSTOP}/SETUP/config.common.${COMPILER}
    MAKES += ${SOLPSTOP}/SETUP/config.common.${COMPILER}
  endif
else
  MAKES += ${SOLPSTOP}/Makefile ${SOLPSTOP}/SETUP/setup.csh.${HOST_NAME}.${COMPILER} ${SOLPSTOP}/SETUP/config.${HOST_NAME}.${COMPILER}
  ifeq ($(shell [ -e ${SOLPSTOP}/SETUP/config.common.${COMPILER} ] && echo yes || echo no ),yes)
    MAKES += ${SOLPSTOP}/SETUP/config.common.${COMPILER}
  endif
endif
ifeq ($(shell [ -e ${SOLPSTOP}/SETUP/setup.csh.${HOST_NAME}.${COMPILER}.local ] && echo yes || echo no ),yes)
  MAKES += ${SOLPSTOP}/SETUP/setup.csh.${HOST_NAME}.${COMPILER}.local
endif
ifeq ($(shell [ -e ${SOLPSTOP}/SETUP/config.${HOST_NAME}.${COMPILER}.local ] && echo yes || echo no ),yes)
  include ${SOLPSTOP}/SETUP/config.${HOST_NAME}.${COMPILER}.local
  MAKES += ${SOLPSTOP}/SETUP/config.${HOST_NAME}.${COMPILER}.local
endif

ifdef USE_B25
  ifndef SOLPSTOP
    $(error SOLPSTOP not defined, but trying to compile with B2.5)
  endif
endif
ifdef SOLPS_MPI
  USE_MPI ?= -DUSE_MPI
else
  ifdef USE_MPI
    SOLPS_MPI = yes
  endif
endif
ifdef SOLPS_OPENMP
  ifndef USE_OPENMP
    USE_OPENMP = -DUSE_OPENMP
  endif
endif

# Default prefix for OBJDIR: standalone
PREF_OBJDIR = standalone
ifdef USE_B25
PREF_OBJDIR = couple_SOLPS-ITER
endif

# Extensions for object directories when various options are used
ifdef USE_IMPGYRO
EXT_IMPGYRO = .ig
else
ifdef USE_MPI
EXT_MPI = .mpi
endif
endif
ifdef USE_OPENMP
EXT_OPENMP = .openmp
endif
ifdef SOLPS_DEBUG
DBG       = 1
EXT_DEBUG = .debug
endif

# Directory where objectcode/binaries will be produced
OBJDIR = ${SRCEIR}/builds/${PREF_OBJDIR}.${HOST_NAME}.${COMPILER}${EXT_OPENMP}${EXT_MPI}${EXT_IMPGYRO}${EXT_DEBUG}

INCLUDE =
# Defines MPI_VERSION, which is the MPI version number
# MPI_VERSION is used in eirmod_mpi to check whether the MPI 3 subroutines are available
ifdef USE_IMPGYRO
  USE_MPI ?= -DUSE_MPI
  INCLUDE = ${MPI_CPP}
  include ${OBJDIR}/mpiversion.mk
else
 ifdef USE_MPI
  INCLUDE = ${MPI_CPP}
  include ${OBJDIR}/mpiversion.mk
 endif
endif
INCLUDE += -I${SRCLOCAL} -I${SRCDIR}
USE_DIMENSIONS = 0
ifdef USE_B25
  ifeq ($(shell [ -s ${B2SRC}/modules/b2mod_dimensions.F ] && echo yes || echo no ),yes)
    USE_DIMENSIONS = 1
    DIMSDIR = ${B2SRC}/modules
    ifeq ($(shell [ -s ${B2SRC}/modules.local/b2mod_dimensions.F ] && echo yes || echo no ),yes)
      DIMSDIR = ${B2SRC}/modules.local
    endif
    DEFINES += -DDIMENSIONS_MODULE
  else
    ifeq ($(shell [ -d ${B2SRC}/include.local ] && echo yes || echo no ),yes)
      INCLUDE += -I${B2SRC}/include.local
    endif
    INCLUDE += -I${B2SRC}/include
  endif
endif
GIT_DESCRIPTION=$(shell git describe --tags --dirty --always | cut -c 1-40 )
GIT_BRANCH=$(shell git branch --show-current )
GIT_DATE=$(shell git show -s --format="%ci" )
GIT_HASH=$(shell git show -s --format="%H" )

ifdef USE_B25
EIRENE_DEFINES += ${USE_B25}
endif
ifdef USE_MPI
EIRENE_DEFINES += ${USE_MPI} -DMPI_VERSION=${MPI_VERSION}
ifdef OPEN_MPI
EIRENE_DEFINES += -DOPEN_MPI=${OPEN_MPI}
endif
LDLIBS += ${LD_MPI}
endif
ifdef USE_OPENMP
EIRENE_DEFINES += ${USE_OPENMP}
endif
ifdef USE_IMPGYRO
EIRENE_DEFINES += ${USE_IMPGYRO}
endif
ifdef SOLPS_DEBUG
EIRENE_DEFINES += -DDBG -DTRACE
endif

#.SUFFIXES: .o .f .f90

ifeq ("$(DBG)","")
     CFLAGS= -O2
else
     SOLPS_DEBUG= yes
     CFLAGS= -g -fbounds-check
endif

ifdef USE_MPI
FPATH=${SRCDIR}/broadcast
else
FPATH=${SRCDIR}/broadcast_dummy
INCLUDE+=-I${SRCDIR}/broadcast_dummy
endif
FPATH+=:${SRCLOCAL}:${SRCDIR}/modules:${SRCDIR}/photons:${SRCDIR}/assistant:${INTDIR}:${SRCDIR}/diagno:${SRCDIR}/file-handling:${SRCDIR}/geometry:${SRCDIR}/iterate:${SRCDIR}/main-routines:${SRCDIR}/mathematics:${SRCDIR}/output:${SRCDIR}/particle-tracing:${SRCDIR}/plotting:${SRCDIR}/sampling:${SRCDIR}/scoring:${SRCDIR}/startup-routines:${SRCDIR}/surface-processes:${SRCDIR}/user-routines/user_iter:${SRCDIR}/volume-processes:${SRCDIR}/tetrahedra
FFPATH=${SRCDIR}/geometry/time-routines
NOGPATH:=${FPATH}:${FFPATH}
ifdef USE_B25
FPATH+=:${SRCDIR}/grsoft-dummy
else
ifndef LD_GR
FPATH+=:${SRCDIR}/grsoft-dummy
endif
endif
F90PATH=${INTDIR}:${SRCDIR}/broadcast:${SRCDIR}/modules
FF90PATH=${INTDIR}:${SRCDIR}/broadcast:${SRCDIR}/modules
VPATH=${FPATH}:${FFPATH}:${F90PATH}

$(shell awk 'FNR==1{if(!/^OBJS *=/){e=1;exit}} END{exit e}' \
        ${OBJDIR}/LISTOBJ 2>/dev/null || rm -f ${OBJDIR}/LISTOBJ)
include ${OBJDIR}/LISTOBJ
include ${SRCEIR}/config/compile
ifeq ($(shell [ -e ${SRCEIR}/config/config.${HOST_NAME}.${COMPILER}.local ] && echo yes || echo no ),yes)
include ${SRCEIR}/config/config.${HOST_NAME}.${COMPILER}.local
MAKES+= ${SRCEIR}/config/config.${HOST_NAME}.${COMPILER}.local
endif

MODULES := ${patsubst ${SRCDIR}/*/eirmod_%.f,%,${shell echo ${SRCDIR}/*/eirmod_*.f}}
MODLOCASE:=$(shell echo $(MODULES)| tr A-Z a-z)
ALLOBJS = ${OBJNOG:%.o=${OBJDIR}/%.o}
ALLOBJS+= ${OBJSF90:%.o=${OBJDIR}/%.o}
ALLOBJS+= ${OBJSFF90:%.o=${OBJDIR}/%.o}
GRSOFT := ${patsubst ${SRCDIR}/grsoft-dummy/%.f,${OBJDIR}/%.o,${shell echo ${SRCDIR}/grsoft-dummy/*.f}}
GRSOFT_EXCLUDE := ${patsubst ${SRCDIR}/grsoft-dummy/%.f,%.o,${shell echo ${SRCDIR}/grsoft-dummy/*.f}}
MPILIST = ${patsubst ${SRCDIR}/broadcast_dummy/%.f,%.o,${shell echo ${SRCDIR}/broadcast_dummy/mpi_*.f}}
EXELIST = eirene_main.o ${GRSOFT_EXCLUDE} ${MPILIST}
MAINLIST = ${addprefix $(OBJDIR)/,${EXELIST}}
EXCLUDELIS = ${EXELIST:.o=\\.o}

testvar = ${addprefix $(OBJDIR)/,${OBJS}}

all: dirs local links VERSION $(OBJDIR)/${PROG} ${OBJDIR}/${PROGX}

.PHONY: dirs clean depend listobj tags echo local links VERSION force

ifeq (${USE_DIMENSIONS},1)
DIMENSIONS:
	@ln -sf ${DIMSDIR}/b2mod_dimensions.F ${INTDIR}
ifneq (${MOD},o)
${OBJDIR}/eirmod_extrab25.${MOD}: ${OBJDIR}/b2mod_dimensions.${MOD}
endif
${OBJDIR}/couple_param.o: ${OBJDIR}/b2mod_dimensions.o
${OBJDIR}/eirmod_extrab25.o: ${OBJDIR}/b2mod_dimensions.o
else
DIMENSIONS:
endif

${OBJDIR}/${PROGX}: ${OBJDIR}/eirene_main.o ${OBJDIR}/$(EIRLIB) ${OBJDIR}/libgr_dummy.a ${MAKES}
	${FC} ${FCOPTS} ${FFLAGSEXTRA} -o ${OBJDIR}/${PROGX} ${OBJDIR}/eirene_main.o ${OBJDIR}/$(EIRLIB) ${OBJDIR}/libgr_dummy.a ${LDLIBS} ${LD_NAG} ${LDEXTRA}

ifdef LD_GR
ifdef LD_GKS
ifeq ($(strip ${GLI_HOME}),)
$(warning Eirene graphics may not work because GLI_HOME is not defined.)
endif
${OBJDIR}/${PROG}: ${OBJDIR}/eirene_main.o ${OBJDIR}/${EIRLIB} ${MAKES}
	${FC} ${FCOPTS} ${FFLAGSEXTRA} -o ${OBJDIR}/${PROG} ${OBJDIR}/eirene_main.o ${OBJDIR}/$(EIRLIB) ${LDLIBS} ${GKS} ${LD_NAG} ${LDEXTRA}
else
$(warning Eirene graphics will be turned off because LD_GKS is not defined.)
${OBJDIR}/${PROG}: ${OBJDIR}/eirene_main.o ${OBJDIR}/$(EIRLIB) ${OBJDIR}/libgr_dummy.a ${MAKES}
	${FC} ${FCOPTS} ${FFLAGSEXTRA} -o ${OBJDIR}/${PROG} ${OBJDIR}/eirene_main.o ${OBJDIR}/$(EIRLIB) ${OBJDIR}/libgr_dummy.a ${LDLIBS} ${LD_NAG} ${LDEXTRA}
endif
else
$(warning Eirene graphics will be turned off because LD_GR is not defined.)
${OBJDIR}/${PROG}: ${OBJDIR}/eirene_main.o ${OBJDIR}/$(EIRLIB) ${OBJDIR}/libgr_dummy.a ${MAKES}
	${FC} ${FCOPTS} ${FFLAGSEXTRA} -o ${OBJDIR}/${PROG} ${OBJDIR}/eirene_main.o ${OBJDIR}/$(EIRLIB) ${OBJDIR}/libgr_dummy.a ${LDLIBS} ${LD_NAG} ${LDEXTRA}
endif

ifdef USE_OPENMP
ifeq ($(COMPILER),ifort64)
${OBJDIR}/read_fixform.o : read_fixform.F
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} $< $*.F
	${FC} ${FCOPTS} -qoverride-limits ${FFLAGSEXTRA} -c ${INCLUDE} ${INCMOD}${OBJDIR} ${OBJDEST} $*.F
	@if [ -f $*.${MOD} ] ; then touch $*.${MOD} ; fi
endif
endif

ifdef USE_MPI
ifeq ($(COMPILER),gfortran)
ifeq ($(shell test ${GFORTRAN_MAJOR_VERSION} -ge 10; echo $$?),0)
${OBJDIR}/broadcast.o : broadcast.F
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} $< $*.F
	${FC} ${FCOPTS} -fallow-argument-mismatch ${FFLAGSEXTRA} -c ${INCLUDE} ${INCMOD}${OBJDIR} ${OBJDEST} $*.F
	@if [ -f $*.${MOD} ] ; then touch $*.${MOD} ; fi

${OBJDIR}/calstr.o : calstr.F
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} $< $*.F
	${FC} ${FCOPTS} -fallow-argument-mismatch ${FFLAGSEXTRA} -c ${INCLUDE} ${INCMOD}${OBJDIR} ${OBJDEST} $*.F
	@if [ -f $*.${MOD} ] ; then touch $*.${MOD} ; fi

${OBJDIR}/collect_census.o : collect_census.F
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} $< $*.F
	${FC} ${FCOPTS} -fallow-argument-mismatch ${FFLAGSEXTRA} -c ${INCLUDE} ${INCMOD}${OBJDIR} ${OBJDEST} $*.F
	@if [ -f $*.${MOD} ] ; then touch $*.${MOD} ; fi

${OBJDIR}/collect_coutau.o : collect_coutau.F
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} $< $*.F
	${FC} ${FCOPTS} -fallow-argument-mismatch ${FFLAGSEXTRA} -c ${INCLUDE} ${INCMOD}${OBJDIR} ${OBJDEST} $*.F
	@if [ -f $*.${MOD} ] ; then touch $*.${MOD} ; fi

${OBJDIR}/eirmod_balanced_strategy.o : eirmod_balanced_strategy.F90
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} $< $*.F90
	${FC} ${FCOPTS} -fallow-argument-mismatch ${FFLAGSEXTRA} -c ${INCLUDE} ${INCMOD}${OBJDIR} ${OBJDEST} $*.F90
	@if [ -f $*.${MOD} ] ; then touch $*.${MOD} ; fi

${OBJDIR}/eirmod_cadgeo.o : eirmod_cadgeo.F
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} $< $*.F
	${FC} ${FCOPTS} -fallow-argument-mismatch ${FFLAGSEXTRA} -c ${INCLUDE} ${INCMOD}${OBJDIR} ${OBJDEST} $*.F
	@if [ -f $*.${MOD} ] ; then touch $*.${MOD} ; fi

${OBJDIR}/eirmod_cai.o : eirmod_cai.F
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} $< $*.F
	${FC} ${FCOPTS} -fallow-argument-mismatch ${FFLAGSEXTRA} -c ${INCLUDE} ${INCMOD}${OBJDIR} ${OBJDEST} $*.F
	@if [ -f $*.${MOD} ] ; then touch $*.${MOD} ; fi

${OBJDIR}/eirmod_calstr_buffered.o : eirmod_calstr_buffered.F90
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} $< $*.F90
	${FC} ${FCOPTS} -fallow-argument-mismatch ${FFLAGSEXTRA} -c ${INCLUDE} ${INCMOD}${OBJDIR} ${OBJDEST} $*.F90
	@if [ -f $*.${MOD} ] ; then touch $*.${MOD} ; fi

${OBJDIR}/eirmod_ccoupl.o : eirmod_ccoupl.F
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} $< $*.F
	${FC} ${FCOPTS} -fallow-argument-mismatch ${FFLAGSEXTRA} -c ${INCLUDE} ${INCMOD}${OBJDIR} ${OBJDEST} $*.F
	@if [ -f $*.${MOD} ] ; then touch $*.${MOD} ; fi

${OBJDIR}/eirmod_cestim.o : eirmod_cestim.F
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} $< $*.F
	${FC} ${FCOPTS} -fallow-argument-mismatch ${FFLAGSEXTRA} -c ${INCLUDE} ${INCMOD}${OBJDIR} ${OBJDEST} $*.F
	@if [ -f $*.${MOD} ] ; then touch $*.${MOD} ; fi

${OBJDIR}/eirmod_cgeom.o : eirmod_cgeom.F
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} $< $*.F
	${FC} ${FCOPTS} -fallow-argument-mismatch ${FFLAGSEXTRA} -c ${INCLUDE} ${INCMOD}${OBJDIR} ${OBJDEST} $*.F
	@if [ -f $*.${MOD} ] ; then touch $*.${MOD} ; fi

${OBJDIR}/eirmod_cgrid.o : eirmod_cgrid.F
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} $< $*.F
	${FC} ${FCOPTS} -fallow-argument-mismatch ${FFLAGSEXTRA} -c ${INCLUDE} ${INCMOD}${OBJDIR} ${OBJDEST} $*.F
	@if [ -f $*.${MOD} ] ; then touch $*.${MOD} ; fi

${OBJDIR}/eirmod_cinit.o : eirmod_cinit.F
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} $< $*.F
	${FC} ${FCOPTS} -fallow-argument-mismatch ${FFLAGSEXTRA} -c ${INCLUDE} ${INCMOD}${OBJDIR} ${OBJDEST} $*.F
	@if [ -f $*.${MOD} ] ; then touch $*.${MOD} ; fi

${OBJDIR}/eirmod_clgin.o : eirmod_clgin.F
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} $< $*.F
	${FC} ${FCOPTS} -fallow-argument-mismatch ${FFLAGSEXTRA} -c ${INCLUDE} ${INCMOD}${OBJDIR} ${OBJDEST} $*.F
	@if [ -f $*.${MOD} ] ; then touch $*.${MOD} ; fi

${OBJDIR}/eirmod_comnnl.o : eirmod_comnnl.F
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} $< $*.F
	${FC} ${FCOPTS} -fallow-argument-mismatch ${FFLAGSEXTRA} -c ${INCLUDE} ${INCMOD}${OBJDIR} ${OBJDEST} $*.F
	@if [ -f $*.${MOD} ] ; then touch $*.${MOD} ; fi

${OBJDIR}/eirmod_comsig.o : eirmod_comsig.F
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} $< $*.F
	${FC} ${FCOPTS} -fallow-argument-mismatch ${FFLAGSEXTRA} -c ${INCLUDE} ${INCMOD}${OBJDIR} ${OBJDEST} $*.F
	@if [ -f $*.${MOD} ] ; then touch $*.${MOD} ; fi

${OBJDIR}/eirmod_comsou.o : eirmod_comsou.F
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} $< $*.F
	${FC} ${FCOPTS} -fallow-argument-mismatch ${FFLAGSEXTRA} -c ${INCLUDE} ${INCMOD}${OBJDIR} ${OBJDEST} $*.F
	@if [ -f $*.${MOD} ] ; then touch $*.${MOD} ; fi

${OBJDIR}/eirmod_comspl.o : eirmod_comspl.F
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} $< $*.F
	${FC} ${FCOPTS} -fallow-argument-mismatch ${FFLAGSEXTRA} -c ${INCLUDE} ${INCMOD}${OBJDIR} ${OBJDEST} $*.F
	@if [ -f $*.${MOD} ] ; then touch $*.${MOD} ; fi

${OBJDIR}/eirmod_comusr.o : eirmod_comusr.F
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} $< $*.F
	${FC} ${FCOPTS} -fallow-argument-mismatch ${FFLAGSEXTRA} -c ${INCLUDE} ${INCMOD}${OBJDIR} ${OBJDEST} $*.F
	@if [ -f $*.${MOD} ] ; then touch $*.${MOD} ; fi

${OBJDIR}/eirmod_comxs.o : eirmod_comxs.F
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} $< $*.F
	${FC} ${FCOPTS} -fallow-argument-mismatch ${FFLAGSEXTRA} -c ${INCLUDE} ${INCMOD}${OBJDIR} ${OBJDEST} $*.F
	@if [ -f $*.${MOD} ] ; then touch $*.${MOD} ; fi

${OBJDIR}/eirmod_cpes.o : eirmod_cpes.F
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} $< $*.F
	${FC} ${FCOPTS} -fallow-argument-mismatch ${FFLAGSEXTRA} -c ${INCLUDE} ${INCMOD}${OBJDIR} ${OBJDEST} $*.F
	@if [ -f $*.${MOD} ] ; then touch $*.${MOD} ; fi

${OBJDIR}/eirmod_cpolyg.o : eirmod_cpolyg.F
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} $< $*.F
	${FC} ${FCOPTS} -fallow-argument-mismatch ${FFLAGSEXTRA} -c ${INCLUDE} ${INCMOD}${OBJDIR} ${OBJDEST} $*.F
	@if [ -f $*.${MOD} ] ; then touch $*.${MOD} ; fi

${OBJDIR}/eirmod_cref.o : eirmod_cref.F
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} $< $*.F
	${FC} ${FCOPTS} -fallow-argument-mismatch ${FFLAGSEXTRA} -c ${INCLUDE} ${INCMOD}${OBJDIR} ${OBJDEST} $*.F
	@if [ -f $*.${MOD} ] ; then touch $*.${MOD} ; fi

${OBJDIR}/eirmod_cstep.o : eirmod_cstep.F
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} $< $*.F
	${FC} ${FCOPTS} -fallow-argument-mismatch ${FFLAGSEXTRA} -c ${INCLUDE} ${INCMOD}${OBJDIR} ${OBJDEST} $*.F
	@if [ -f $*.${MOD} ] ; then touch $*.${MOD} ; fi

${OBJDIR}/eirmod_ctetra.o : eirmod_ctetra.F
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} $< $*.F
	${FC} ${FCOPTS} -fallow-argument-mismatch ${FFLAGSEXTRA} -c ${INCLUDE} ${INCMOD}${OBJDIR} ${OBJDEST} $*.F
	@if [ -f $*.${MOD} ] ; then touch $*.${MOD} ; fi

${OBJDIR}/eirmod_ctrcei.o : eirmod_ctrcei.F
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} $< $*.F
	${FC} ${FCOPTS} -fallow-argument-mismatch ${FFLAGSEXTRA} -c ${INCLUDE} ${INCMOD}${OBJDIR} ${OBJDEST} $*.F
	@if [ -f $*.${MOD} ] ; then touch $*.${MOD} ; fi

${OBJDIR}/eirmod_ctrig.o : eirmod_ctrig.F
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} $< $*.F
	${FC} ${FCOPTS} -fallow-argument-mismatch ${FFLAGSEXTRA} -c ${INCLUDE} ${INCMOD}${OBJDIR} ${OBJDEST} $*.F
	@if [ -f $*.${MOD} ] ; then touch $*.${MOD} ; fi

${OBJDIR}/eirmod_infcop.o : eirmod_infcop.F
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} $< $*.F
	${FC} ${FCOPTS} -fallow-argument-mismatch ${FFLAGSEXTRA} -c ${INCLUDE} ${INCMOD}${OBJDIR} ${OBJDEST} $*.F
	@if [ -f $*.${MOD} ] ; then touch $*.${MOD} ; fi

${OBJDIR}/eirmod_photon.o : eirmod_photon.F
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} $< $*.F
	${FC} ${FCOPTS} -fallow-argument-mismatch ${FFLAGSEXTRA} -c ${INCLUDE} ${INCMOD}${OBJDIR} ${OBJDEST} $*.F
	@if [ -f $*.${MOD} ] ; then touch $*.${MOD} ; fi

${OBJDIR}/eirmod_wneutrals.o : eirmod_wneutrals.F
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} $< $*.F
	${FC} ${FCOPTS} -fallow-argument-mismatch ${FFLAGSEXTRA} -c ${INCLUDE} ${INCMOD}${OBJDIR} ${OBJDEST} $*.F
	@if [ -f $*.${MOD} ] ; then touch $*.${MOD} ; fi

${OBJDIR}/eirsrt.o : eirsrt.F
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} $< $*.F
	${FC} ${FCOPTS} -fallow-argument-mismatch ${FFLAGSEXTRA} -c ${INCLUDE} ${INCMOD}${OBJDIR} ${OBJDEST} $*.F
	@if [ -f $*.${MOD} ] ; then touch $*.${MOD} ; fi

ifneq ($(MOD),o)
${OBJDIR}/eirmod_balanced_strategy.{MOD} : eirmod_balanced_strategy.F90
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} $< $*.F90
	${FC} ${FCOPTS} -fallow-argument-mismatch ${FFLAGSEXTRA} -c ${INCLUDE} ${INCMOD}${OBJDIR} ${OBJDEST} $*.F90

${OBJDIR}/eirmod_cadgeo.{MOD} : eirmod_cadgeo.F
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} $< $*.F
	${FC} ${FCOPTS} -fallow-argument-mismatch ${FFLAGSEXTRA} -c ${INCLUDE} ${INCMOD}${OBJDIR} ${OBJDEST} $*.F

${OBJDIR}/eirmod_cai.{MOD} : eirmod_cai.F
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} $< $*.F
	${FC} ${FCOPTS} -fallow-argument-mismatch ${FFLAGSEXTRA} -c ${INCLUDE} ${INCMOD}${OBJDIR} ${OBJDEST} $*.F

${OBJDIR}/eirmod_calstr_buffered.{MOD} : eirmod_calstr_buffered.F90
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} $< $*.F90
	${FC} ${FCOPTS} -fallow-argument-mismatch ${FFLAGSEXTRA} -c ${INCLUDE} ${INCMOD}${OBJDIR} ${OBJDEST} $*.F90

${OBJDIR}/eirmod_ccoupl.{MOD} : eirmod_ccoupl.F
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} $< $*.F
	${FC} ${FCOPTS} -fallow-argument-mismatch ${FFLAGSEXTRA} -c ${INCLUDE} ${INCMOD}${OBJDIR} ${OBJDEST} $*.F

${OBJDIR}/eirmod_cestim.{MOD} : eirmod_cestim.F
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} $< $*.F
	${FC} ${FCOPTS} -fallow-argument-mismatch ${FFLAGSEXTRA} -c ${INCLUDE} ${INCMOD}${OBJDIR} ${OBJDEST} $*.F

${OBJDIR}/eirmod_cgeom.{MOD} : eirmod_cgeom.F
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} $< $*.F
	${FC} ${FCOPTS} -fallow-argument-mismatch ${FFLAGSEXTRA} -c ${INCLUDE} ${INCMOD}${OBJDIR} ${OBJDEST} $*.F

${OBJDIR}/eirmod_cgrid.{MOD} : eirmod_cgrid.F
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} $< $*.F
	${FC} ${FCOPTS} -fallow-argument-mismatch ${FFLAGSEXTRA} -c ${INCLUDE} ${INCMOD}${OBJDIR} ${OBJDEST} $*.F

${OBJDIR}/eirmod_cinit.{MOD} : eirmod_cinit.F
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} $< $*.F
	${FC} ${FCOPTS} -fallow-argument-mismatch ${FFLAGSEXTRA} -c ${INCLUDE} ${INCMOD}${OBJDIR} ${OBJDEST} $*.F

${OBJDIR}/eirmod_clgin.{MOD} : eirmod_clgin.F
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} $< $*.F
	${FC} ${FCOPTS} -fallow-argument-mismatch ${FFLAGSEXTRA} -c ${INCLUDE} ${INCMOD}${OBJDIR} ${OBJDEST} $*.F

${OBJDIR}/eirmod_comnnl.{MOD} : eirmod_comnnl.F
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} $< $*.F
	${FC} ${FCOPTS} -fallow-argument-mismatch ${FFLAGSEXTRA} -c ${INCLUDE} ${INCMOD}${OBJDIR} ${OBJDEST} $*.F

${OBJDIR}/eirmod_comsig.{MOD} : eirmod_comsig.F
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} $< $*.F
	${FC} ${FCOPTS} -fallow-argument-mismatch ${FFLAGSEXTRA} -c ${INCLUDE} ${INCMOD}${OBJDIR} ${OBJDEST} $*.F

${OBJDIR}/eirmod_comsou.{MOD} : eirmod_comsou.F
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} $< $*.F
	${FC} ${FCOPTS} -fallow-argument-mismatch ${FFLAGSEXTRA} -c ${INCLUDE} ${INCMOD}${OBJDIR} ${OBJDEST} $*.F

${OBJDIR}/eirmod_comspl.{MOD} : eirmod_comspl.F
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} $< $*.F
	${FC} ${FCOPTS} -fallow-argument-mismatch ${FFLAGSEXTRA} -c ${INCLUDE} ${INCMOD}${OBJDIR} ${OBJDEST} $*.F

${OBJDIR}/eirmod_comusr.{MOD} : eirmod_comusr.F
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} $< $*.F
	${FC} ${FCOPTS} -fallow-argument-mismatch ${FFLAGSEXTRA} -c ${INCLUDE} ${INCMOD}${OBJDIR} ${OBJDEST} $*.F

${OBJDIR}/eirmod_comxs.{MOD} : eirmod_comxs.F
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} $< $*.F
	${FC} ${FCOPTS} -fallow-argument-mismatch ${FFLAGSEXTRA} -c ${INCLUDE} ${INCMOD}${OBJDIR} ${OBJDEST} $*.F

${OBJDIR}/eirmod_cpes.{MOD} : eirmod_cpes.F
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} $< $*.F
	${FC} ${FCOPTS} -fallow-argument-mismatch ${FFLAGSEXTRA} -c ${INCLUDE} ${INCMOD}${OBJDIR} ${OBJDEST} $*.F

${OBJDIR}/eirmod_cpolyg.{MOD} : eirmod_cpolyg.F
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} $< $*.F
	${FC} ${FCOPTS} -fallow-argument-mismatch ${FFLAGSEXTRA} -c ${INCLUDE} ${INCMOD}${OBJDIR} ${OBJDEST} $*.F

${OBJDIR}/eirmod_cref.{MOD} : eirmod_cref.F
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} $< $*.F
	${FC} ${FCOPTS} -fallow-argument-mismatch ${FFLAGSEXTRA} -c ${INCLUDE} ${INCMOD}${OBJDIR} ${OBJDEST} $*.F

${OBJDIR}/eirmod_cstep.{MOD} : eirmod_cstep.F
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} $< $*.F
	${FC} ${FCOPTS} -fallow-argument-mismatch ${FFLAGSEXTRA} -c ${INCLUDE} ${INCMOD}${OBJDIR} ${OBJDEST} $*.F

${OBJDIR}/eirmod_ctetra.{MOD} : eirmod_ctetra.F
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} $< $*.F
	${FC} ${FCOPTS} -fallow-argument-mismatch ${FFLAGSEXTRA} -c ${INCLUDE} ${INCMOD}${OBJDIR} ${OBJDEST} $*.F

${OBJDIR}/eirmod_ctrcei.{MOD} : eirmod_ctrcei.F
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} $< $*.F
	${FC} ${FCOPTS} -fallow-argument-mismatch ${FFLAGSEXTRA} -c ${INCLUDE} ${INCMOD}${OBJDIR} ${OBJDEST} $*.F

${OBJDIR}/eirmod_ctrig.{MOD} : eirmod_ctrig.F
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} $< $*.F
	${FC} ${FCOPTS} -fallow-argument-mismatch ${FFLAGSEXTRA} -c ${INCLUDE} ${INCMOD}${OBJDIR} ${OBJDEST} $*.F

${OBJDIR}/eirmod_infcop.{MOD} : eirmod_infcop.F
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} $< $*.F
	${FC} ${FCOPTS} -fallow-argument-mismatch ${FFLAGSEXTRA} -c ${INCLUDE} ${INCMOD}${OBJDIR} ${OBJDEST} $*.F

${OBJDIR}/eirmod_photon.{MOD} : eirmod_photon.F
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} $< $*.F
	${FC} ${FCOPTS} -fallow-argument-mismatch ${FFLAGSEXTRA} -c ${INCLUDE} ${INCMOD}${OBJDIR} ${OBJDEST} $*.F

${OBJDIR}/eirmod_wneutrals.{MOD} : eirmod_wneutrals.F
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} $< $*.F
	${FC} ${FCOPTS} -fallow-argument-mismatch ${FFLAGSEXTRA} -c ${INCLUDE} ${INCMOD}${OBJDIR} ${OBJDEST} $*.F
endif
endif
endif
endif

$(OBJDIR)/$(EIRLIB): ${ALLOBJS} ${SRCDIR}/git_version_Eirene.h
	@ar uvr $@ ${ALLOBJS}
	ranlib $@

${OBJDIR}/libgr_dummy.a: ${GRSOFT}
	ar r $@ $?
	ranlib $@

clean:
	rm -f $(OBJDIR)/*.f $(OBJDIR)/*.f90 $(OBJDIR)/*.o $(OBJDIR)/*.txt $(OBJDIR)/*.in $(OBJDIR)/*.cmake $(OBJDIR)/*.a $(OBJDIR)/${PROG} ${OBJDIR}/${PROGX} ${SRCDIR}/git_version_Eirene.h ${OBJDIR}/dependencies ${OBJDIR}/LISTOBJ ${OBJDIR}/Makefile ${OBJDIR}/mpiversion.mk ${OBJDIR}/mpi_version ${SRCEIR}/builds/binRelease/mpi_version.mk ${SRCEIR}/builds/binRelease/mpi_version.f90 ${SRCEIR}/builds/binRelease/mpi_version
ifneq (${MOD},o)
	rm -f ${OBJDIR}/*.${MOD}
endif
	rm -Rf ${OBJDIR}/CMakeFiles
ifndef NO_MANUAL
	cd ${SRCEIR}/Manual ; rm -f *.aux *.log *.dvi *.pdf *.ps *.out *.toc *.bbl *.blg
endif

clean_manual:
	cd ${SRCEIR}/Manual ; rm -f *.aux *.log *.dvi *.pdf *.ps *.out *.toc *.bbl *.blg

depend: ${SRCLOCAL}/eirlocal.f ${OBJSF:.o=.f} ${OBJSFF:.o=.F} ${OBJSF90:.o=.f90} ${OBJSFF90:.o=.F90} eirene_main.F
	@`which makedepend` -p'$${OBJDIR}/' ${DEFINES} -f- ${INCLUDE} $^ | \
	sed 's,^$${OBJDIR}/[^ ][^ ]*/,\$${OBJDIR}/,' | \
        sed 's,: ${SOLPSTOP},: $${SOLPSTOP},' > ${OBJDIR}/dependencies
	@echo '# 1' >> ${OBJDIR}/dependencies
ifneq (${MOD},o)
	@`which makedepend` -p'$${OBJDIR}/' ${DEFINES} -f- ${INCLUDE} ${SRCDIR}/*/eirmod_*.f ${SRCDIR}/*/eirmod_*.[fF]90 ${INTDIR}/eirmod_*.[fF] ${INTDIR}/eirmod_*.[fF]90 -o.${MOD} | \
	sed 's,^$${OBJDIR}/[^ ][^ ]*/,\$${OBJDIR}/,' | \
        sed 's,: ${SOLPSTOP},: $${SOLPSTOP},' >> ${OBJDIR}/dependencies
	@echo '# 2' >> ${OBJDIR}/dependencies
endif
	@egrep -aiH '^[ 	]+use ' ${SRCDIR}/*/*.[fF] ${SRCDIR}/*/*.[fF]90 | grep -v 'IGNORE' | tr , ' '  |  awk '{sub(/\.f:|\.F:|\.f90:|\.F90:/,".o:",$$1); sub("^.*/","$${OBJDIR}/",$$1); print $$1,"$${OBJDIR}/"tolower($$3)".${MOD}"}' >> ${OBJDIR}/dependencies
	@echo '# 3a' >> ${OBJDIR}/dependencies
	@egrep -aiH '^[ 	]+use ' ${INTDIR}/*.[fF] ${INTDIR}/*.[fF]90 ${SRCDIR}/user-routines/user_iter/*.[fF] ${SRCDIR}/geometry/time-routines/*.F | grep -v 'IGNORE' | tr , ' '  | awk '{sub(/\.f:|\.F:|\.f90:|\.F90:/,".o:",$$1);sub("^.*/","$${OBJDIR}/",$$1); print $$1,"$${OBJDIR}/"tolower($$3)".${MOD}"}' >> ${OBJDIR}/dependencies
	@echo '# 4a' >> ${OBJDIR}/dependencies
	@egrep -aiH -m 1 '^[ 	]*module ' ${SRCDIR}/*/*.[fF] ${SRCDIR}/*/*.[fF]90 ${INTDIR}/*.[Ff] ${INTDIR}/*.[Ff]90 ${SRCDIR}/user-routines/user_iter/*.[fF] ${SRCDIR}/geometry/time-routines/*.F | sed 's/:module/: module/' | tr , ' ' | awk '{sub(/\.f:|\.F:|\.f90:|\.F90:/,".o", $$1); sub("^.*/","$${OBJDIR}/",$$1); sub(",","",$$3); print "$${OBJDIR}/" tolower($$3) ".${MOD}:", $$1;}' >> ${OBJDIR}/dependencies
	@echo '# 5a' >> ${OBJDIR}/dependencies
ifneq (${MOD},o)
	@egrep -aiH '^[ 	]+use ' ${SRCDIR}/*/*.[fF] ${SRCDIR}/*/*.[fF]90 | grep -v 'IGNORE' | tr , ' '  |  awk '{sub(/\.f:|\.F:|\.f90:|\.F90:/,".o:",$$1); sub("^.*/","$${OBJDIR}/",$$1); print $$1,"$${OBJDIR}/"tolower($$3)".o"}' >> ${OBJDIR}/dependencies
	@echo '# 3b' >> ${OBJDIR}/dependencies
	@egrep -aiH '^[ 	]+use ' ${INTDIR}/*.[fF] ${INTDIR}/*.[fF]90 ${SRCDIR}/user-routines/user_iter/*.[fF] ${SRCDIR}/geometry/time-routines/*.F | grep -v 'IGNORE' | tr , ' '  | awk '{sub(/\.f:|\.F:|\.f90:|\.F90:/,".o:",$$1);sub("^.*/","$${OBJDIR}/",$$1); print $$1,"$${OBJDIR}/"tolower($$3)".o"}' >> ${OBJDIR}/dependencies
	@echo '# 4b' >> ${OBJDIR}/dependencies
	@egrep -aiH '^[ 	]+use ' ${INTDIR}/eirmod_*.[fF] ${INTDIR}/eirmod_*.[fF]90 ${SRCDIR}/*/eirmod_*.[fF] ${SRCDIR}/*/eirmod_*.[Ff]90 | grep -v 'IGNORE' | tr , ' '  | awk '{sub(/\.f:|\.F:|\.f90:|\.F90:/,".${MOD}:",$$1);sub("^.*/","$${OBJDIR}/",$$1); print $$1,"$${OBJDIR}/"tolower($$3)".${MOD}"}' >> ${OBJDIR}/dependencies
	@echo '# 5b' >> ${OBJDIR}/dependencies
endif

listobj: DIMENSIONS ${SRCDIR}/assistant/version.f local | $(OBJDIR)
	@rm -f ${OBJDIR}/LISTOBJ ; touch ${OBJDIR}/LISTOBJ ; l="OBJS =" ; \
	for d in `echo "$(FPATH)" | tr : \ `; do \
		l="$$l `(cd $$d > /dev/null; echo *.f)`"; \
	done; \
	for d in `echo "$(FFPATH)" | tr : \ `; do \
		l="$$l `(cd $$d > /dev/null; echo *.F)`"; \
	done; \
	E="-e 's/ \*\.f//g' -e 's/ \*\.F//g' -e 's/\.f/\.o/g' -e 's/\.F/\.o/g'" ; for f in $(EXCLUDELIS); do \
		E="$$E -e 's/ $$f//'"; \
	done; \
	echo "$$l" | eval sed "$$E" > ${OBJDIR}/LISTOBJ
	@lf="OBJSF =" ; \
	for d in `echo "$(FPATH)" | tr : \ `; do \
		lf="$$lf `(cd $$d > /dev/null; echo *.f)`"; \
	done; \
	E="-e 's/ \*\.f//g' -e 's/\.f/\.o/g'" ; for f in $(EXCLUDELIS); do \
		E="$$E -e 's/ $$f//'"; \
	done; \
	echo "$$lf" | eval sed "$$E" >> ${OBJDIR}/LISTOBJ
	@lff="OBJSFF =" ; \
	for d in `echo "$(FFPATH)" | tr : \ `; do \
		lff="$$lff `(cd $$d > /dev/null; echo *.F)`"; \
	done; \
	E="-e 's/ \*\.F//g' -e 's/\.F/\.o/g'" ; for f in $(EXCLUDELIS); do \
		E="$$E -e 's/ $$f//'"; \
	done; \
	echo "$$lff" | eval sed "$$E" >> ${OBJDIR}/LISTOBJ
	@lg="OBJNOG =" ; \
	for d in `echo "$(NOGPATH)" | tr : \ `; do \
		lg="$$lg `(cd $$d > /dev/null; echo *.[fF] )`"; \
	done; \
	E="-e 's/ \*\.f//g' -e 's/ \*\.F//g' -e 's/\.f/\.o/g' -e 's/\.F/\.o/g'" ; \
	echo "$$lg" | eval sed "$$E" >> ${OBJDIR}/LISTOBJ
	@ll="OBJSF90 =" ; \
	for d in `echo "$(F90PATH)" | tr : \ `; do \
		ll="$$ll `(cd $$d > /dev/null; echo *.f90)`"; \
	done; \
	E="-e 's/ \*\.f90//g' -e 's/\.f90/\.o/g'" ; \
	echo "$$ll" | eval sed "$$E" >> ${OBJDIR}/LISTOBJ
	@lll="OBJSFF90 =" ; \
	for d in `echo "$(FF90PATH)" | tr : \ `; do \
		lll="$$lll `(cd $$d > /dev/null; echo *.F90)`"; \
	done; \
	E="-e 's/ \*\.F90//g' -e 's/\.F90/\.o/g'" ; \
	echo "$$lll" | eval sed "$$E" >> ${OBJDIR}/LISTOBJ

tags:
	rm -f ${SRCEIR}/TAGS ; ${MAKETAGS} ${SRCEIR}/TAGS ${SRCLOCAL}/*.[fF] ${SRCDIR}/*/*.[Ff] ${SRCDIR}/*/*.[Ff]90 ${SRCDIR}/*/*/*.[Ff]90 ${INTDIR}/*.[fF] ${SRCDIR}/user-routines/user_iter/*.[fF] ${SRCDIR}/geometry/time-routines/*.F ${SRCEIR}/Manual/eirene.tex ${SRCEIR}/Manual/tex/*.tex || touch ${SRCEIR}/TAGS

${OBJDIR}/LISTOBJ: listobj

VERSION: ${SRCDIR}/git_version_Eirene.h ${SRCDIR}/assistant/version.f

${SRCDIR}/assistant/version.f: ${SRCDIR}/cmake/version.f.in ${MAKES}
	cat ${SRCDIR}/cmake/version.f.in | sed -e 's|@git_description@|${GIT_DESCRIPTION}|' -e 's|@git_date@|${GIT_DATE}|' -e 's|@git_hash@|${GIT_HASH}|' -e 's|@git_branch@|${GIT_BRANCH}|' > ${SRCDIR}/assistant/version.f

${SRCDIR}/git_version_Eirene.h: force
	@echo "      character*32 :: git_version_Eirene =" > ${SRCDIR}/git_version_new.h
	@echo "     . '`git describe --tags --dirty --always | cut -c 1-32`'" >> ${SRCDIR}/git_version_new.h
	@if cmp -s ${SRCDIR}/git_version_new.h ${SRCDIR}/git_version_Eirene.h; then rm ${SRCDIR}/git_version_new.h; else mv ${SRCDIR}/git_version_new.h ${SRCDIR}/git_version_Eirene.h; fi

${OBJDIR}/dependencies: ${SRCDIR}/modules/.new_modules
	-mkdir -p ${OBJDIR}
	printf '# Dummy dependencies file for Eirene\n' > ${OBJDIR}/dependencies
	${MAKE} -f ${SRCEIR}/config/Makefile tags
	${MAKE} -f ${SRCEIR}/config/Makefile VERSION
	${MAKE} -f ${SRCEIR}/config/Makefile links
	${MAKE} -f ${SRCEIR}/config/Makefile local
	${MAKE} -f ${SRCEIR}/config/Makefile listobj
	${MAKE} -f ${SRCEIR}/config/Makefile depend

$(shell [ -s ${OBJDIR}/dependencies ] || rm -f ${OBJDIR}/dependencies)
include ${OBJDIR}/dependencies
ifeq ($(shell [ -e ${SRCEIR}/config/dependencies.local ] && echo yes || echo no ),yes)
include ${SRCEIR}/config/dependencies.local
endif

${SRCEIR}/TAGS: tags

dirs: $(OBJDIR)

$(OBJDIR):
	mkdir -p $(OBJDIR)

local: ${SRCLOCAL}/eirlocal.f

${SRCLOCAL}/eirlocal.f:
	mkdir -p ${SRCLOCAL}
	echo "      subroutine eirlocal" > ${SRCLOCAL}/eirlocal.f
	echo "c" >> ${SRCLOCAL}/eirlocal.f
	echo "c store local or locally modified subroutines in this directory" >> ${SRCLOCAL}/eirlocal.f
	echo "c" >> ${SRCLOCAL}/eirlocal.f
	echo "      end" >> ${SRCLOCAL}/eirlocal.f

links: ${SRCEIR}/Database/AMJUEL ${SRCEIR}/Database/AMMONX ${SRCEIR}/Database/H2VIBR ${SRCEIR}/Database/HYDHEL ${SRCEIR}/Database/METHANE ${SRCEIR}/Database/METHAN

${SRCEIR}/Database/AMJUEL:
	-ln -s ${SRCEIR}/Database/AMdata/amjuel.tex ${SRCEIR}/Database/AMJUEL

${SRCEIR}/Database/AMMONX:
	-ln -s ${SRCEIR}/Database/AMdata/AMMONX_Arrh-elast.tex ${SRCEIR}/Database/AMMONX

${SRCEIR}/Database/H2VIBR:
	-ln -s ${SRCEIR}/Database/AMdata/h2vibr.tex ${SRCEIR}/Database/H2VIBR

${SRCEIR}/Database/HYDHEL:
	-ln -s ${SRCEIR}/Database/AMdata/hydhel.tex ${SRCEIR}/Database/HYDHEL

${SRCEIR}/Database/METHANE:
	-ln -s ${SRCEIR}/Database/AMdata/methane.tex ${SRCEIR}/Database/METHANE

${SRCEIR}/Database/METHAN: ${SRCEIR}/Database/METHANE
	-ln -s ${SRCEIR}/Database/METHANE ${SRCEIR}/Database/METHAN

${OBJDIR}/mpiversion.mk: ${MAKES} | $(OBJDIR)
ifdef NO_MPI
	echo 'MPI_VERSION=0' > ${OBJDIR}/mpiversion.mk
else
	printf "use mpi\ninteger OPEN_MPI\nWRITE(*,fmt='(A12,I1)') 'MPI_VERSION=', MPI_VERSION\nWRITE(*,fmt='(A9)') 'MPI_MOD=1'\nif (OPEN_MPI.ne.0) WRITE(*,fmt='(A10)') 'OPEN_MPI=1'\nEND\n" > ${OBJDIR}/mpi_version.f90
	( ${MPI_FC} ${FCOPTS} ${INCLUDE} -o ${OBJDIR}/mpi_version ${OBJDIR}/mpi_version.f90 ${LD_MPI} && ( ${OBJDIR}/mpi_version | tail -n3 ) || \
	( printf "include 'mpif.h'\nWRITE(*,fmt='(A12,I1)') 'MPI_VERSION=', MPI_VERSION\nWRITE(*,fmt='(A9)') 'MPI_MOD=0'\nEND\n" > ${OBJDIR}/mpi_version.f90 ; \
	( ${MPI_FC} ${FCOPTS} ${INCLUDE} -o ${OBJDIR}/mpi_version ${OBJDIR}/mpi_version.f90 ${LD_MPI} && ( ${OBJDIR}/mpi_version | tail -n2 ) || ( echo MPI_VERSION=0 ; echo MPI_MOD=0 ) ) ) ) > ${OBJDIR}/mpiversion.mk 2> /dev/null
endif

echo:
	@echo GRSOFT=${GRSOFT}
#	@echo GRSOFT_EXCLUDE=${GRSOFT_EXCLUDE}
#	@echo PLOT=${PLOT}
#	@echo PLOT_EXCLUDE=${PLOT_EXCLUDE}
#	@echo ${SRCDIR}/grsoft-dummy/*.F
#	@echo SRCLOCAL=${SRCLOCAL}:
#	@echo LOCAL=${LOCAL}
	@echo FPATH=${FPATH}
	@echo NOGPATH=${NOGPATH}
	@echo ALLOBJS=${ALLOBJS}
	@echo OBJS=${OBJS}
	@echo OBJNOG=${OBJNOG}
	@echo OBJSF90=${OBJSF90}
	@echo OBJSFF90=${OBJSFF90}
#	@echo DEST=${DEST}
#	@echo MODULES=${MODULES}
#	@echo MODLOCASE=${MODLOCASE}
	@echo EXCLUDELIS=${EXCLUDELIS}
	@echo MAINLIST=${MAINLIST}
#	@echo SRCF=${SRCF}
#	@echo OBJDEST=${OBJDEST}
	@echo LD_GR=${LD_GR}
	@echo FCOPTS=${FCOPTS}
