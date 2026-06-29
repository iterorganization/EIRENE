# CMake project file for EIRENE

cmake_minimum_required(VERSION 3.10)

project(EIRENE Fortran)

# Set a default build type if none was specified.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE "Release" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Release" "Debug" "Coverage")
  # "Coverage" only available for GNU compiler
endif()

# Adds compiler flags to all build types.
if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
  set(CMAKE_Fortran_FLAGS_RELEASE "-fp-model source -fpp")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -fpp -gen-interfaces -warn all,interfaces -check all -traceback -fpe0 -diag-disable:10182")
elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "IntelLLVM")
  set(CMAKE_Fortran_FLAGS_RELEASE "-fp-model source -fpp")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -fpp -gen-interfaces -warn all,interfaces -check all,nouninit -traceback -fpe0 -diag-disable:10182,10440")
  string(APPEND CMAKE_EXE_LINKER_FLAGS_DEBUG "-check all,nouninit")
elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -cpp -frecursive")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -cpp -frecursive -Wall -fcheck=all -Wtarget-lifetime")
  if (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER 8.1)
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -Wdo-subscript")
  endif()
  set(CMAKE_Fortran_FLAGS_COVERAGE "-g -cpp -frecursive -fprofile-arcs -ftest-coverage")
  add_definitions(-DGFORTRAN)
elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "PGI" OR "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "NVHPC")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -cpp")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -cpp -C -Mchkptr -Mchkstk")
elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "NAG")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O -x77 -compatible -fpp")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g90 -fpp -x77 -compatible")
elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Cray")
  set(CMAKE_Fortran_FLAGS_RELEASE "-e0 -ez -ef -em -O3 -hfp3")
  set(CMAKE_Fortran_FLAGS_DEBUG "-e0 -ez -ef -em -O0 -g -hfp0 -K trap=divz -K trap=fp -G0 -Rb")
elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "G95")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O -cpp")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -cpp -C -fbounds-check -ftrace=full")
endif()
if (FLAGS_EXTRA)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FLAGS_EXTRA}")
endif()
message(STATUS "Using ${CMAKE_Fortran_COMPILER_ID} compiler with version ${CMAKE_Fortran_COMPILER_VERSION}")

# Output names and postfix.
if(DEFINED ENV{SOLPS_CPP})
set(EIRENE_LIBRARY_NAME "eirene" CACHE STRING "EIRENE library name")
set(EIRENE_EXECUTABLE_NAME "eirobj" CACHE STRING "EIRENE executable name")
set(SOLPS_POSTFIX "" CACHE STRING "SOLPS executable names postfix")
else()
set(EIRENE_LIBRARY_NAME "EIRENE" CACHE STRING "EIRENE library name")
set(EIRENE_EXECUTABLE_NAME "eirene" CACHE STRING "EIRENE executable name")
endif()
set(EIRENE_POSTFIX "" CACHE STRING "EIRENE output names postfix")

# Check for modules and library for JSON
# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# This is a bit convoluted but allows for some user error in the JSON library location
# First check for user defined paths
if(LibJSON)
  message(STATUS " Searching for user defined LibJSON. Found " ${LibJSON} )
  if(IS_DIRECTORY ${LibJSON})
    if((EXISTS "${LibJSON}/libjsonfortran.a" OR EXISTS "${LibJSON}/libjsonfortran.so"))
      message(STATUS " Found JSON library in user defined LibJSON path. ")
      if(NOT JSON_MODULES)
	set(JSON_MODULES ${LibJSON})
        if(NOT (EXISTS "${JSON_MODULES}/json_module.mod"))
          string(CONCAT JSON_MODULES ${LibJSON} "/../include")
        endif()
      endif()
      set(LibJSON "${LibJSON}/libjsonfortran.a")
      set(jsonfortran_FOUND TRUE)
    else()
      message(WARNING "JSON Library libjsonfortran.a NOT found at " ${LibJSON})
    endif()
  else()  # Assume the library is included in the path
    if(EXISTS ${LibJSON})
      message(STATUS " Found user defined LibJSON. ")
      set(jsonfortran_FOUND TRUE)
    endif()
    if(NOT JSON_MODULES)
      string(REPLACE "libjsonfortran.a" "" MOD_TRY ${LibJSON})
      if(EXISTS "${MOD_TRY}json_module.mod")
        set(JSON_MODULES ${MOD_TRY})
      else()
        string(CONCAT MOD_TRY2 ${MOD_TRY} "../include")
        if(EXISTS "${MOD_TRY2}json_module.mod")
          set(JSON_MODULES ${MOD_TRY2})
        endif()
      endif()
      if(NOT JSON_MODULES)
        string(REPLACE "libjsonfortran.so" "" MOD_TRY ${LibJSON})
        if(EXISTS "${MOD_TRY}json_module.mod")
          set(JSON_MODULES ${MOD_TRY})
        else()
          string(CONCAT MOD_TRY2 ${MOD_TRY} "../include")
          if(EXISTS "${MOD_TRY2}json_module.mod")
            set(JSON_MODULES ${MOD_TRY2})
          endif()
        endif()
      endif()
    endif()
  endif()
  if(EXISTS "${JSON_MODULES}/json_module.mod")
    message(STATUS " Found user defined JSON modules.")
    set(jsonfortran_INCLUDE_DIRS ${JSON_MODULES})
  else()
    message(STATUS "JSON Modules not found at " ${JSON_MODULES})
    message(FATAL_ERROR "The -DJSON_MODULES variable must be set to the directory containing JSON Modules.")
  endif()
else()
  if(DEFINED ENV{JSONF_LIB_DIR})
    message(STATUS " JSON library location defined by environment variable JSONF_LIB_DIR: " $ENV{JSONF_LIB_DIR})
    if(EXISTS $ENV{JSONF_LIB_DIR}/libjsonfortran.a)
      set(LibJSON $ENV{JSONF_LIB_DIR}/libjsonfortran.a)
      set(jsonfortran_FOUND TRUE)
    elseif(EXISTS $ENV{JSONF_LIB_DIR}/libjsonfortran.so)
      set(LibJSON $ENV{JSONF_LIB_DIR}/libjsonfortran.so)
      set(jsonfortran_FOUND TRUE)
    else()
      message(FATAL_ERROR "JSON library libjsonfortran.a or libjsonfortran.so not found in location specified by environment variable JSONF_LIB_DIR")
    endif()
    message(STATUS " LibJSON: " ${LibJSON})
    if(DEFINED ENV{JSONF_MOD_DIR})
      message(STATUS " JSON library modules location defined by environment variable JSONF_MOD_DIR: " $ENV{JSONF_MOD_DIR})
      set(JSON_MODULES $ENV{JSONF_MOD_DIR})
    else()
      set(JSON_MODULES $ENV{JSONF_LIB_DIR})
    endif()
    if(EXISTS "${JSON_MODULES}/json_module.mod")
      message(STATUS " Found JSON modules defined by environment variable.")
      set(jsonfortran_INCLUDE_DIRS ${JSON_MODULES})
    endif()
  else()
    message(STATUS " Detecting LibJSON. By default the static library will be used.")
    if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
      find_package(jsonfortran NAMES jsonfortran-intel)
    elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "IntelLLVM")
      find_package(jsonfortran NAMES jsonfortran-intelllvm)
    elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
      find_package(jsonfortran NAMES jsonfortran-gnu)
    else()
      message(STATUS " Compiler version not known for JSON library")
    endif()
    if(jsonfortran_FOUND)
      if(EXISTS "${jsonfortran_INCLUDE_DIRS}/libjsonfortran.a")
        set(LibJSON "${jsonfortran_INCLUDE_DIRS}/libjsonfortran.a")
      elseif(EXISTS "${jsonfortran_INCLUDE_DIRS}/libjsonfortran.so")
        set(LibJSON "${jsonfortran_INCLUDE_DIRS}/libjsonfortran.so")
      endif()
      if( NOT LibJSON)
        string(REPLACE "include" "lib" LibJSON ${jsonfortran_INCLUDE_DIRS})
        if(EXISTS "${LibJSON}/libjsonfortran.a")
          set(LibJSON "${LibJSON}/libjsonfortran.a")
        elseif(EXISTS "${LibJSON}/libjsonfortran.so")
          set(LibJSON "${LibJSON}/libjsonfortran.so")
        else()
          message(FATAL_ERROR "JSON Library not found at " ${LibJSON})
        endif()
      endif()
      message(STATUS " Found LibJSON. Version: " ${jsonfortran_VERSION})
    endif()
  endif()
endif()

if(jsonfortran_FOUND)
  message(STATUS " JSON_MODULES path " ${jsonfortran_INCLUDE_DIRS})
  include_directories(${jsonfortran_INCLUDE_DIRS})
else()
  message(STATUS "JSON library not found.")
  message(STATUS "Use -DLibJSON=\"<directory>\" to set it if available.")
  message(STATUS "Code will be compiled with NO_JSON option.")
  add_definitions(-DNO_JSON)
endif()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Control of interface and user-routines to compile.
if( NOT EIRENE_INTERFACE OR NOT EIRENE_USER-ROUTINES)
  message(STATUS " --  Setting interface and user-routines to compile:")
endif()
if( NOT EIRENE_INTERFACE)
  message(STATUS "       If flag -DEIRENE_INTERFACE is unspecified, \"Dummy\"")
  message(STATUS "       will be set. Try e.g. -DEIRENE_INTERFACE=SOLPS-ITER.  ")
endif()
if( NOT EIRENE_USER-ROUTINES)
  message(STATUS "       If flag -DEIRENE_USER-ROUTINES is unspecified, \"default\"")
  message(STATUS "       will be set. Try e.g. -DEIRENE_USER-ROUTINES=iter.  ")
endif()
set(EIRENE_INTERFACE "Dummy" CACHE STRING "EIRENE interface to compile")
set(EIRENE_USER-ROUTINES "default" CACHE STRING "EIRENE user-routines to compile")

# Add option to control whether trace output will be compiled.
option(TRACE "Build with trace output" OFF)
if (TRACE)
  message(STATUS "-- Adding TRACE flag...")
  add_definitions(-DTRACE)
endif()

# Add option to control whether binary files are also written in ASCII.
option(CHECKBIN "Check binary files" OFF)
if (CHECKBIN)
  message(STATUS "-- Adding CHECKBIN flag...")
  add_definitions(-DCHECKBIN)
endif()

# Add option to control whether compiling coupled to B2.5
option(B25_EIRENE "Build coupled B2.5-Eirene executable" OFF)
if (B25_EIRENE)
  message(STATUS "-- Adding B25_EIRENE flag...")
  add_definitions(-DB25_EIRENE)
  if (IS_DIRECTORY $ENV{SOLPSTOP}/modules/B2.5/src/include.local)
    include_directories($ENV{SOLPSTOP}/modules/B2.5/src/include.local)
  endif()
  include_directories($ENV{SOLPSTOP}/modules/B2.5/src/include)
  if (DIMENSIONS_MODULE)
    add_definitions(-DDIMENSIONS_MODULE)
    if(EXISTS "$ENV{SOLPSTOP}/modules/B2.5/src/modules.local/b2mod_dimensions.F")
      include_directories($ENV{SOLPSTOP}/modules/B2.5/src/modules.local)
      set(B2MOD_DIMENSIONS_SRC "$ENV{SOLPSTOP}/modules/B2.5/src/modules.local/b2mod_dimensions.F")
    else()
      include_directories($ENV{SOLPSTOP}/modules/B2.5/src/modules)
      set(B2MOD_DIMENSIONS_SRC "$ENV{SOLPSTOP}/modules/B2.5/src/modules/b2mod_dimensions.F")
    endif()
  endif()
  set(GR-DUMMY_LIBRARY_NAME "gr_dummy" CACHE STRING "GR-DUMMY library name")
endif()
if(DEFINED ENV{SOLPS_CPP})
  if( NOT "$ENV{SOLPS_CPP}" STREQUAL "")
    message(STATUS "-- Adding SOLPS-ITER environment flags...")
    string(REPLACE "-DGFORTRAN" "" SOLPS_DEFS $ENV{SOLPS_CPP})
    add_definitions(${SOLPS_DEFS})
  else()
    message(STATUS "-- Adding Fortran 2003 flag...")
    add_definitions(-DF2003)
  endif()
else()
  message(STATUS "-- Adding Fortran 2003 flag...")
  add_definitions(-DF2003)
endif()
if(DEFINED ENV{EIRENE_DEFINES})
  if( NOT "$ENV{EIRENE_DEFINES}" STREQUAL "")
    message(STATUS "-- Adding specific Eirene definitions...")
    add_definitions($ENV{EIRENE_DEFINES})
  endif()
endif()

# Add option to control whether MPI will be used.
option(MPI "Build with MPI" ON)
if (MPI)
  message(STATUS "-- Adding USE_MPI and MPI_VERSION flags...")
  add_definitions(-DUSE_MPI)
  if (MPI_VERSION)
    add_definitions(-DMPI_VERSION=${MPI_VERSION})
  endif()
  if (OPEN_MPI)
    message(STATUS "-- Adding OPEN_MPI flag...")
    add_definitions(-DOPEN_MPI=${OPEN_MPI})
  endif()
  find_package(MPI REQUIRED)
  if ( NOT MPI_FOUND )

  endif()
  include_directories(${MPI_Fortran_INCLUDE_PATH})

  if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
    if (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER 9.5)
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch")
    endif()
  endif()
endif()

# Add option to enable OpenMP with the parallel region defined in eirmod_mcarlo.
option(OPENMP "Build for use with internal parallel region" OFF)
if (OPENMP)
  cmake_minimum_required(VERSION 3.1...3.10 FATAL_ERROR)
  message(STATUS "-- Adding USE_OPENMP flag...")
  add_definitions(-DUSE_OPENMP)
  find_package(OpenMP REQUIRED)
  if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "IntelLLVM")
    message(STATUS "-- OpenMP is unsafe with ifx compiler! Do not use. --")
  endif()
endif()

# Add option to enable OpenMP with the parallel region defined in the calling program.
option(EXT_OPENMP "Build for use with external parallel region" OFF)
if (EXT_OPENMP)
  cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
  message(STATUS "-- Adding USE_EXT_OPENMP flag...")
  add_definitions(-DUSE_EXT_OPENMP)
  find_package(OpenMP REQUIRED)
endif()

# Adds openMP using the old style CMake approach to ensure backwards compatibility
# This should be changed at some point to use the newer OpenMP Fortran capabilities
if (OPENMP_FOUND)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
  if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qoverride-limits")
# elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Cray")
#   set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -M 3140")
  endif()
endif()

# Finds all fortran files (*.[fF]/*.[fF]90) in the source directory and extract the
# relative path. Relative path is stored in variable 'fortran_dirs'.
message(STATUS "Detecting fortran files (*.[fF]/*.[fF]90) in subdirectories of "
               ${PROJECT_SOURCE_DIR})
# The git-version source (version.f) is generated per build into this build's
# own binary directory (see the 'version' custom target below) rather than into
# the shared source tree.  Generating it under PROJECT_SOURCE_DIR caused a
# parallel-build race: 'make -jN all' configures two CMake build trees that
# share this same source directory (the standalone Eirene executable and the
# B2.5-coupled EIRENE library), and both would write and compile
# assistant/version.f concurrently, corrupting the build.
set(EIRENE_VERSION_F ${PROJECT_BINARY_DIR}/version.f)
if(NOT EXISTS "${EIRENE_VERSION_F}")
  file(TOUCH ${EIRENE_VERSION_F})
endif()
set_source_files_properties(${EIRENE_VERSION_F} PROPERTIES GENERATED TRUE)
file(GLOB_RECURSE fortran_files RELATIVE ${PROJECT_SOURCE_DIR}
     ${PROJECT_SOURCE_DIR}/*.[fF]
     ${PROJECT_SOURCE_DIR}/*.[fF]90)
set(fortran_dirs "")
foreach(fortran_file ${fortran_files})
  get_filename_component(fortran_dir ${fortran_file} PATH)
  list(APPEND fortran_dirs ${fortran_dir})
endforeach()
list(REMOVE_DUPLICATES fortran_dirs)

# Sets variables 'src_[relative path to subdirectory]' which contain a list of
# all fortran files (*.[fF]/*.[fF]90) in this subdirectory.
message(STATUS "Configuring variables src_[relative path to subdirectory]")
foreach(fortran_dir ${fortran_dirs})
  file(GLOB src_${fortran_dir} ${PROJECT_SOURCE_DIR}/${fortran_dir}/*.[fF]
                               ${PROJECT_SOURCE_DIR}/${fortran_dir}/*.[fF]90)
endforeach()

# version.f is now generated into the build dir (EIRENE_VERSION_F) and added to
# the EIRENE target explicitly below; drop any stale source-tree copy left by an
# older layout or by a traditional (NO_CMAKE) build so it is not compiled twice.
list(REMOVE_ITEM src_assistant ${PROJECT_SOURCE_DIR}/assistant/version.f)

# Remove program file to be able to add it separately if required.
list(GET src_main-routines 0 fortran_file)
get_filename_component(fortran_dir ${fortran_file} PATH)
list(REMOVE_ITEM src_main-routines ${fortran_dir}/eirene_main.F)

# Specify interface and user-routines to use.
set(src_interface/x ${src_interfaces/couple_${EIRENE_INTERFACE}})
set(src_user-routines/x ${src_user-routines/user_${EIRENE_USER-ROUTINES}})
if( NOT src_interface/x )
  message(FATAL_ERROR "No interface named \"${EIRENE_INTERFACE}\"")
endif()
if( NOT src_user-routines/x )
  message(FATAL_ERROR "No user-routines named \"${EIRENE_USER-ROUTINES}\"")
endif()

# Add option to control whether graphics will be enabled.
option(GRAPHICS "Build with graphics enabled" OFF)

# Add linking to the plotting libraries
set(plotpath "")
if( NOT B25_EIRENE )

  add_executable(eirene main-routines/eirene_main.F)

  if(GRAPHICS)

    message(STATUS " --  Setting libraries for internal plotting:  --  --")

    if(LibGKS)
      if (EXISTS ${LibGKS})
        get_filename_component(LibGKS ${LibGKS} ABSOLUTE)
        set(LibGKS_local "yes")
      else()
        message(STATUS "Specified LibGKS was not found")
        unset(LibGKS)
      endif()
    else()
      find_library(LibGKS gks)
      if(LibGKS)
        message(STATUS "Using system GKS library " ${LibGKS})
        message(STATUS "To change it, use: -DLibGKS=$PATH_TO_GKS_LIB/libgks.a")
      else()
        message(STATUS "Path to LibGKS unset")
        message(STATUS "You may want to try -DLibGKS=$PATH_TO_GKS_LIB/libgks.a")
      endif()
    endif()

    if(LibGRS)
      if (EXISTS ${LibGRS})
        get_filename_component(LibGRS ${LibGRS} ABSOLUTE)
        set(LibGRS_local "yes")
      else()
        message(STATUS "Specified LibGRS was not found")
        unset(LibGRS)
      endif()
    else()
      find_library(LibGRS gr)
      if(LibGRS)
        message(STATUS "Using system GR library " ${LibGRS})
        message(STATUS "To change it, use: -DLibGRS=$PATH_TO_GRSOFT_LIB/libgr.a")
      else()
        message(STATUS "Path to LibGRS unset")
        message(STATUS "You may want to try -DLibGRS=$PATH_TO_GRSOFT_LIB/libgr.a")
      endif()
    endif()

    if(LibX11)
      if (EXISTS ${LibX11})
        get_filename_component(LibX11 ${LibX11} ABSOLUTE)
        set(LibX11_local "yes")
      else()
        message(STATUS "Specified LibX11 was not found")
        unset(LibX11)
      endif()
    else()
      find_library(LibX11 X11)
      if(LibX11)
        message(STATUS "Using system X11 library " ${LibX11})
        message(STATUS "To change it, use: -DLibX11=$PATH_TO_X11_LIB/libX11.so")
      else()
        message(STATUS "Path to LibX11 unset")
        message(STATUS "You may want to try -DLibX11=$PATH_TO_X11_LIB/libX11.so")
      endif()
    endif()

    if(LibXt)
      if (EXISTS ${LibXt})
        get_filename_component(LibXt ${LibXt} ABSOLUTE)
        set(LibXt_local "yes")
      else()
        message(STATUS "Specified LibXt was not found")
        unset(LibXt)
      endif()
    else()
      find_library(LibXt Xt)
      if(LibXt)
        message(STATUS "Using system libXt library " ${LibXt})
        message(STATUS "To change it, use: -DLibXt=$PATH_TO_Xt_LIB/libXt.so")
      else()
        message(STATUS "Path to LibXt unset")
        message(STATUS "You may want to try -DLibXt=$PATH_TO_Xt_LIB/libXt.so")
      endif()
    endif()

    if(Libz)
      if (EXISTS ${Libz})
        get_filename_component(Libz ${Libz} ABSOLUTE)
        set(Libz_local "yes")
      else()
        message(STATUS "Specified Libz was not found")
        unset(Libz)
      endif()
    else()
      find_library(Libz z)
      if(Libz)
        message(STATUS "Using system libz library " ${Libz})
        message(STATUS "To change it, use: -DLibz=$PATH_TO_z_LIB/libz.so")
      else()
        message(STATUS "Path to Libz unset")
        message(STATUS "You may want to try -DLibz=$PATH_TO_z_LIB/libz.so")
      endif()
    endif()

    if(LibGKS AND LibGRS AND LibX11 AND Libz AND LibXt)
      message(STATUS "Found the following libraries used for plotting:")
      message(STATUS "  Lib GKS: " ${LibGKS})
      message(STATUS "  Lib GRS: " ${LibGRS})
      message(STATUS "  Lib X11: " ${LibX11})
      message(STATUS "  Lib Xt:  " ${LibXt})
      message(STATUS "  Lib z:   " ${Libz})

      if(DEFINED ENV{GRSOFT_DEVICE})
        string(FIND $ENV{GRSOFT_DEVICE} "62" GRSOFT_OK)
      endif()
      if( NOT ( DEFINED ENV{GLI_HOME} AND ( GRSOFT_OK GREATER -1 ) ) )
        message(STATUS "For running, please export the following statements "
                       "into your console:")
        if(NOT ( DEFINED ENV{GLI_HOME} ) )
          message(STATUS "export GLI_HOME=$PATH_TO_GKS_LIB")
        endif()
        if( GRSOFT_OK EQUAL -1 )
          message(STATUS "export GRSOFT_DEVICE=62")
        endif()
      endif()
      message(STATUS " --  Plotting will be enabled.")
    endif()

  else()
    message(STATUS " --  --  Set -DGRAPHICS=ON to enable plotting.")
  endif()

endif()

# Decide whether to use plotting or plot_dummy.
if(LibGRS AND LibGKS AND LibX11 AND Libz AND LibXt AND GRAPHICS)
  list(APPEND plotpath ${src_plotting})
else()
  if(DEFINED ENV{SOLPS_CPP})
    list(APPEND plotpath ${src_plotting} ${src_grsoft-dummy})
  else()
    list(APPEND plotpath ${src_plot_dummy} ${src_grsoft-dummy})
  endif()
  message(STATUS " --  Plotting will be disabled.")
  set(SOLPS_POSTFIX "x")
endif()

# Adds EIRENE library to the project.
add_library(EIRENE ${src_assistant} ${EIRENE_VERSION_F} ${src_broadcast} ${src_diagno}
                   ${src_file-handling} ${src_function_modules} ${src_geometry}
                   ${src_geometry/time-routines} ${src_interface/x}
                   ${src_iterate} ${src_main-routines} ${src_mathematics}
                   ${src_modules} ${src_output} ${src_particle-tracing}
                   ${src_photons} ${plotpath} ${src_sampling}
                   ${src_scoring} ${src_startup-routines}
                   ${src_surface-processes} ${src_tetrahedra}
                   ${src_user-routines/x} ${src_volume-processes})

if(B25_EIRENE)
  add_library(GR-DUMMY ${src_grsoft-dummy})
  # When coupling to B2.5 with the DIMENSIONS_MODULE pragma, b2mod_dimensions.F
  # (owned by B2.5) must be compiled as part of the EIRENE library so that the
  # .mod file is available to eirmod_extrab25.F90 in this same build directory.
  if(DIMENSIONS_MODULE AND B2MOD_DIMENSIONS_SRC)
    get_target_property(_eirene_sources EIRENE SOURCES)
    list(REMOVE_ITEM _eirene_sources
      ${PROJECT_SOURCE_DIR}/interfaces/couple_${EIRENE_INTERFACE}/b2mod_dimensions.F)
    set_target_properties(EIRENE PROPERTIES SOURCES "${_eirene_sources}")
    target_sources(EIRENE PRIVATE ${B2MOD_DIMENSIONS_SRC})
    set_source_files_properties(
      ${PROJECT_SOURCE_DIR}/interfaces/couple_${EIRENE_INTERFACE}/eirmod_extrab25.F90
      PROPERTIES OBJECT_DEPENDS ${B2MOD_DIMENSIONS_SRC})
  endif()
endif()

if( NOT B25_EIRENE )
  target_link_libraries(eirene EIRENE)
# Link it to the JSON library.
  if (LibJSON)
    target_link_libraries(eirene ${LibJSON})
  endif()
# Link it to the MPI library.
  if (MPI)
    target_link_libraries(eirene ${MPI_Fortran_LIBRARIES})
  endif()
# Link it to the graphics libraries.
  if(LibGRS_local)
    target_link_libraries(eirene ${LibGRS})
  endif()
  if(LibGKS_local)
    target_link_libraries(eirene ${LibGKS})
  endif()
  if(LibX11_local)
    target_link_libraries(eirene ${LibX11})
  endif()
  if(LibXt_local)
    target_link_libraries(eirene ${LibXt})
  endif()
  if(Libz_local)
    target_link_libraries(eirene ${Libz})
  endif()
# Link to any other additional libraries.
  if (ADD_LIBS)
    target_link_libraries(eirene ${ADD_LIBS})
  endif()
# Include any additional linker options
  if (LINK_OPTIONS)
    target_link_options(eirene PUBLIC "LINKER:" ${LINK_OPTIONS})
  endif()

endif()

# This is the newer approach for adding OpenMP - ver >= 3.9
#if(OPENMP OR EXT_OPENMP)
#  target_compile_options(EIRENE PRIVATE ${OpenMP_Fortran_FLAGS})
#  target_compile_options(eirene PRIVATE ${OpenMP_Fortran_FLAGS})
#  target_link_libraries(eirene EIRENE ${OpenMP_Fortran_LIBRARIES})
#endif()

# Sets output directory of the library and executable and output name
# of the executable.
if(DEFINED ENV{SOLPS_CPP})
  set_target_properties(EIRENE PROPERTIES
                      ARCHIVE_OUTPUT_DIRECTORY
                      ${PROJECT_BINARY_DIR}
                      OUTPUT_NAME ${EIRENE_LIBRARY_NAME})
  if ( NOT B25_EIRENE )
    set_target_properties(eirene PROPERTIES
                      RUNTIME_OUTPUT_DIRECTORY
                      ${PROJECT_BINARY_DIR}
                      OUTPUT_NAME ${EIRENE_EXECUTABLE_NAME}${SOLPS_POSTFIX})
  else()
    set_target_properties(GR-DUMMY PROPERTIES
                      ARCHIVE_OUTPUT_DIRECTORY
                      ${PROJECT_BINARY_DIR}
                      OUTPUT_NAME ${GR-DUMMY_LIBRARY_NAME})
    add_dependencies(GR-DUMMY EIRENE)
  endif()
else()
  set_target_properties(EIRENE PROPERTIES
                      ARCHIVE_OUTPUT_DIRECTORY
                      ${PROJECT_BINARY_DIR}/../lib${CMAKE_BUILD_TYPE}
                      OUTPUT_NAME ${EIRENE_LIBRARY_NAME}${EIRENE_POSTFIX})
  if ( NOT B25_EIRENE )
    set_target_properties(eirene PROPERTIES
                      RUNTIME_OUTPUT_DIRECTORY
                      ${PROJECT_BINARY_DIR}/../bin${CMAKE_BUILD_TYPE}
                      OUTPUT_NAME ${EIRENE_EXECUTABLE_NAME}${EIRENE_POSTFIX})
  endif()
endif()

# Adds optional local dependencies
if(EXISTS "${PROJECT_SOURCE_DIR}/../src.local/CMakeLists.txt")
  string(APPEND LOCAL_OBJDIR "${PROJECT_BINARY_DIR}/CMakeFiles/EIRENE.dir")
  include(${PROJECT_SOURCE_DIR}/../src.local/CMakeLists.txt)
endif()

# Add dependencies on local environment settings
set_property(TARGET EIRENE PROPERTY OBJECT_DEPENDS ${PROJECT_SOURCE_DIR}/CMakeLists.txt)
if(DEFINED ENV{SOLPS_CPP})
  string(REPLACE " " ";" SOLPS_Makefiles "${MAKES}")
  set_property(TARGET EIRENE PROPERTY
               OBJECT_DEPENDS ${SOLPS_Makefiles})
endif()

# Gets the git revision hash and stores it in the variable 'git_hash'.
find_package(Git)

if(GIT_FOUND)
  if(EXISTS ${PROJECT_SOURCE_DIR}/../.git)
    set(git_hash_command log --pretty=format:"%ci" --abbrev-commit -1)
    execute_process(COMMAND ${GIT_EXECUTABLE} ${git_hash_command}
                    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/..
                    OUTPUT_VARIABLE git_date
                    OUTPUT_STRIP_TRAILING_WHITESPACE)

    set(git_hash_command log --pretty=format:"%H" --abbrev-commit -1)
    execute_process(COMMAND ${GIT_EXECUTABLE} ${git_hash_command}
                    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/..
                    OUTPUT_VARIABLE git_hash
                    OUTPUT_STRIP_TRAILING_WHITESPACE)

    set(git_hash_command log --pretty=format:"%d" --abbrev-commit -1)
    execute_process(COMMAND ${GIT_EXECUTABLE} ${git_hash_command}
                    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/..
                    OUTPUT_VARIABLE git_branch
                    OUTPUT_STRIP_TRAILING_WHITESPACE)

    set(git_hash_command log --pretty=format:"%ci %H %d" --abbrev-commit -1)
    execute_process(COMMAND ${GIT_EXECUTABLE} ${git_hash_command}
                    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/..
                    OUTPUT_VARIABLE git_full_hash
                    OUTPUT_STRIP_TRAILING_WHITESPACE)

    set(git_hash_command describe --tags --dirty --always)
    execute_process(COMMAND ${GIT_EXECUTABLE} ${git_hash_command}
                    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/..
                    OUTPUT_VARIABLE git_description
                    OUTPUT_STRIP_TRAILING_WHITESPACE)

    # Configures the version CMakeLists file.
    configure_file(${PROJECT_SOURCE_DIR}/cmake/version.cmake.in
                   version.cmake
                   @ONLY)

    # Trigger the EIRENE version output.
    add_custom_target(version
                      ${CMAKE_COMMAND} -D PROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}
                                       -P version.cmake)
    add_dependencies(EIRENE version)
    if ( NOT B25_EIRENE )
      add_dependencies(eirene version)
    endif()
  endif()
endif()

# Adds a target to generate a documentation with Doxygen.
find_package(Doxygen)

if(DOXYGEN_FOUND)
  configure_file(${PROJECT_SOURCE_DIR}/cmake/Doxyfile.in
                 ${PROJECT_BINARY_DIR}/../doc/Doxyfile @ONLY)

  add_custom_target(doc
                   ${DOXYGEN_EXECUTABLE} Doxyfile
                   WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/../doc
                   COMMENT "Generating documentation with Doxygen" VERBATIM)

  message(STATUS "-----------------------------------------------------")
  message(STATUS "To generate the documentation, type \"make doc\"")
  message(STATUS "-----------------------------------------------------")
endif()

