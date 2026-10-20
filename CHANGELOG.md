# Rules for the Changelog

- Add new entries at the top of the file
- Include Author, branch and last commit sha-1 hash
- Include the last EIRENE release version the development is based on
- The commit sha-1 hash is obviously not known at the time of the commit, so best practice is to add an additional commit just including an update to the CHANGELOG.md file. This commit can then reference the previous commits sha-1 hash.
- Include list of changes with short description
- If test cases were updated describe the necessity and the amount by which results have changed

---

---

# Changelog of EIRENE repository:
---

## R. Futtersack - feature/JINTRAC_userEDGE2D - a561ec0839c9af6dc3140e5099388d1a9881c9f3

** Changes **

- Take changes from SOLPS_develop
- Update EDGE2D specific users files
- to be merged to SOLPS_develop as base of JINTRAC_develop branch

## X. Bonnin - SOLPS_develop - bc2a2055bc9e2086dec38e452ccde59f31c62a4a

### Based on SOLPS_wide_grids - fed40a626b2730ec9a0fbf987833905f70237b8f

** Changes **

- Fixed X11 library location in PITAGORA host configuration files (with contribution from E. Tonello)
- Pre-evaluated pkg-config statements to reduce compilation time
- Cleaned up CMake configuration files
- Fixed bug preventing correct compilation in debug mode
- Made CMake build self-contained; improved compiler detection logic (ifx vs. ifort fallback)
- Compile b2mod_dimensions.F as part of EIRENE library to ensure module availability (with contribution from S. Pinches)
- Added IntelLLVM compiler and JSON-Fortran CMake package detection; added safeties for mis-specified graphical library paths
- Added error message for OpenMP incompatibility with ifx compiler; made error non-fatal to allow 2025b toolchain builds to proceed
- Fixed parallel-build race condition: generate version.f into build directory instead of source tree (with contribution from S. Pinches)
- Increased minimum CMake version requirement to 3.10; removed CMake deprecation warnings
- Corrected LaTeX syntax error in Manual
- Fixed fort.44 wall-loading output title field sizes (with contribution from J. Lore)
- Added JINTRAC coupling interface to the SOLPS-ITER/Eirene integration (with contribution from R. Futtersack)
- Fixed W-0279 "Cannot inline - data type mismatch" compiler warnings with NVHPC 24.1: replaced AIMAG() intrinsic calls with Fortran 2003 %IM component selectors for complex variables in eirmod_photon.F
- Eliminated remaining NVHPC 24.1 inliner warnings by moving EIRENE_ph_faddeeva2 from module-level to internal procedure of EIRENE_ph_lorvdw

No CI cases updated.

---

## X. Bonnin - SOLPS_wide_grids

### Based on SOLPS_event_limit - 99e131dfef3a448d0ccb512459f1a55bf211c6fd

** Changes **

- Added SOLPS_WG interface directory for SOLPS-ITER Wide Grids code version
- Automatic detection of SOLPS-ITER version for traditional Makefile
- Automatic build of version.f file
- Made output from JSON and fixed format input file more consistent
- Bug fix for write-up of user_data.input file when there are trailing lines in the input.dat file
- Bug fix for reading LKINDP in JSON formatted input
- Added SNE_PIEL output for SOLPS-ITER balances
- Added L_MACRO quantity in COMXS module to govern condensation for hybrid fluid-kinetic model
- Added IFCTRI array in CTRIG module to store triangle faces for SOLPS-ITER Wide Grids interface
- Removed misleading warning in LEARC1 when particles are on surface
- Made file extensions consistent across user and interface directories
- Updated configuration files for traditional compilation
- Expanded comments
- Output prettifying
- Minor code clean-up

---

## X. Bonnin and A. Pshenov - SOLPS_event_limit

### Based on SOLPS_nprll_wrstrt - f5dcd68af8af6fdc736b93f159a1a9cf71d75a6d

** Changes **

- Implemented an event number limit for individual trajectories

---

## X. Bonnin - SOLPS_nprll_wrstrt

### Based on SOLPS_vector_tallies - bd537444943622e1e37d50acb55e18aadb999e01

** Changes **

- Removed condition to not write WRSTRT files for EMBARRASSINGLY PARALLEL option

---

## X. Bonnin and D. Reiter - SOLPS_vector_tallies

### Based on SOLPS_DR_changes - 7f59f4c18f601a0a35a33d85866f2d5a486f4ef1

** Changes **

- Implemented vectorial momentum tallies

---

## X. Bonnin - SOLPS_DR_changes - 0638d7e91f10d5ce09f561a0e9991df19b6f43e8

### Based on SOLPS_coupling - db157c759dbe59d340df6a0aca48200c1c024fe4

** Changes **

- Other minor changes Detlev wanted
- Modified treatment for generation limit
- Preparatory work for new photon transport variables
- Improvements to BGK background treatment
- Removed polarization reflection file usage
- Removing obsolete Helium emissivity routines
- Code clean-up from D. Reiter and new arguments in EIRENE_SLOPE

---

## X. Bonnin - SOLPS_coupling - 37fb1f2f592d4d15ccdfe3ec68894c09b4d7e5dc

### Based on SOLPS_bug_fixes - 765b12f7f7ab978e8340d215d330ec49c6c5417b

** Changes **

- Implemented new sputtering tallies resolved by incident species causing the sputtering
- Added support for legacy compilers (-DLEGACYCOMP) and G95
- Bug fix for data mapping in snowflake geometries
- New IMPROVED MPI parallelization strategy
- Bug fixes for time-dependent runs and ensuring we always have a filled-in background for coupled runs
- Debugged hybrid OpenMP+MPI parallelization
- Added passing of ZI to SAMSRF
- Output prettifying
- Added 1 more significant digit to fort.44 and fort.46 output files
- Modifed CMakeLists.txt to handle compilation within and for SOLPS-ITER coupling environment
- Added configuration files for CMAKE-less compilation
- Added NO_JSON pragma to allow for compilation without JSON library
- Avoid duplicated write-up of user_data.input JSON file
- Updated version file to include SOLPS-ITER commit information
- Introduced DIMENSIONS_MODULE pragma to grab problem dimensions from SOLPS-ITER environment
- Added NFLA to fort.44 output
- Correction to post-couple test to better handle time census stratum

---

## X. Bonnin and D. Reiter - SOLPS_bug_fixes - 8f586da88c5113442f53dc42c865d8e40daa8e6a

### Based on V1.2.0 (develop - b0027581a7c95e4fc92671f7a354b6c556eacb23)

** Changes **

- Bug fix for computation of momentum source associated with inhomogeneous heavy-heavy collisions
- OpenMP fix for MARLOWE reflection model
- Do not close the master log output file to avoid problems with coupled runs
- Added some input safeties
- Corrected output in SOLPS-ITER interface for EMISS, EMISSMOL, EDISSML, and SRCML arrays
- Improved treatment of time census stratum (only use DELTAT if non-zero, correct broadcast of census stratum quantities, enforced consistency of related inputs)
- Bug fixes to tally print-outs for special cases
- Performance improvements for MPI parallelization
- Added OPEN_MPI pragma to better manage compiler warnings
- Variable type corrections
- Renamed some internal variables to avoid clashes
- Added IOLD argument to EIRENE_UPDATE_SPTFLX
- Added ISTRA argument to EIRENE_SAMVL1
- Corrected units for spectra output + added message about NFILEL and fort.13 format + adapted output formats
- Corrected comments and typos in printout statements
- Removed trailing whitespaces and alignment corrections
- Minor code clean-up

---

## D. Harting - New EIRENE version: Release-V1.2.0 - 1ae241d96213ee7f22f83da44a15d1a561d2f476

#### Based on V1.1.2 and develop 580f5ab43fe927c23fb42daee40cd6f3026c4719

New_EIRENE_Version=1.2.0

```
+===============================================+
|                                               |
|   New EIRENE version: Release-V1.2.0          |
|   Date: 23th May 2025                         |
|   MsV:  MsV-JuelichOrigin_1.0-Oct2024-V1.1.0  |
|                                               |
+===============================================+
```
**Most relevant changes since V1.1.0**

1. EIRENE is now under the dual license:
   - EIRENE Public Licence - "EPL" (see EPL.md) provided by the Forschungszentrum Juelich GmbH (FZJ)
   - and Creative Commons Public License "CC BY-NC-ND 4.0" (see CC_BY-NC-ND.md file and its pdf-converted version)
   The user is bound by at least one of those licenses - at their own choice. The EPL allows research purposed development including creating an own developer community and even providing limited commercial services within the boudaries regulated by the EPL, however it requires registration at www.eirene.de/EPL
2. Minor bugfixes to EIRENE Database.
   - amjuel.tex (Reaction 2.7.14:  e + N_2^+ -> N + N)
   - AMMONX_Arrh-elast.tex (Reaction 14s: e + N2 -> N + N + e)
3. New option NFLAGS > 10:
   - order the triangle faces belonging to one spatially resolved surface into a continuous surface for the output (EIRENE input block 11).
   - See EIRENE manual chapter 2.11.
4. Added rescaling of flux stored on the time boundary when changing the timestep in time depending calculations.
5. Bug fixes related bring EIRENE version inline with SOLPS-ITER version.
   - For details see X. Bonnin's list of changes for commits c3b449d6e828 and 781960eea7b5 in this CHANGELOG.md.
6. Fixed some CI issues.

---

## D. Harting - New EIRENE version: Release-V1.1.2 - 5b961fd9f8718634046f55f828a21fc66d220cb8

#### Based on V1.1.1 - Hotfix for CI

New_EIRENE_Version=1.1.2

```
+===============================================+
|                                               |
|   New EIRENE version: Release-V1.1.2          |
|   Date: 7th April 2025                        |
|   Hotfix                                      |
|                                               |
+===============================================+
```
**Changes**

1. Scheduled CI pipelines (run daily) to update GIT forks at CEA, ITER and JET
was triggering also the 'coverage_report' step. Reordering the rules for the
'coverage_report' step should fix this problem.


---

## D. Harting - New EIRENE version: Release-V1.1.1 - a3df2da6f46bfca0c28cd69ab771145a888a74b5

#### Based on V1.1.0 - Hotfix for CI

New_EIRENE_Version=1.1.1

```
+===============================================+
|                                               |
|   New EIRENE version: Release-V1.1.1          |
|   Date: 27th November 2024                    |
|   Hotfix                                      |
|                                               |
+===============================================+
```
**Changes**

1. Scheduled CI pipelines (run daily) to update GIT forks at CEA, ITER and JET were failing (yaml-error) due to wrong rules setup for some jobs in .gitlab-ci.yml file. This is now fixed.

---
## P. Boerner - PB/bugfixes_and_improvements - 9744e33a0f57c55dc11725a6df83013e909e16bb

### Based on V1.1.0 (develop - 549f0660027fa64b28fb59d11ba951d6759f5f5e)

** Changes **

- In subroutine DBG_PRINTOUT use correct dimensions NPLSTI and NPLSV when calculating averages of TIIN, VXIN, VYIN, VZIN as those are the species limits of these tallies.
- In time-dependent calculations it is necessary to rescale the flux stored on the time boundary when changing the timestep. The rescaled flux needs to be used for FLXCEN as well.
- Set FLXCEN to rescaled flux when changing timestep size in time-dependent mode.
- In EIRENE_OUTFLX for triangular meshes not based on an underlying structured grid printing of spatially resolved surface tallies was not done correctly. This issue has been fixed.
- Additionally an option to order the triangle faces belonging to one spatially resolved surface into a continuous surface has been added (NFLAGS > 10). A table holding information such as coordinates of the start and endpoints of the triangle faces, the arc length along the surface and the surface tally for all requested species will be printed. Produce a sorted list of triangle edges for spatially resolved surfaces and use it in the output
- For triangular meshes the NFLAGS flag has been extended. Here NFLAGS > 10 triggers the production of an additional table holding information about the spatially resolved surface namely start and end points of the triangle edges forming the surface, the arc length along the surface and the values of the surface tally for the requested species.

---

## X. Bonnin - feature/SOLPS_push_bug-fixes - c3b449d6e828da5e17a0d88cf981f7fbd2d84b8c

#### Based on V1.1.0 (develop - bbf7bd63111955e513bedd25188e3a008e9b4928)

** Changes **

- Bug fixes in AMJUEL and AMMONX data files
- Bug fixes for JSON linking and I/O in SOLPS cases
- Added LGEXIT variable in OpenMP treatment for error handling
- Added 0.1 eV minimum value for input electron temperature in PLASMA_DERIV
- Corrected handling of multiple output files for parallel runs
- Corrected array bounds for data read loop in couple_B2 interface
- Corrected location of barrier in BROADCAST routine
- Changed handling of electromagnetic field input tally flags
- Corrected treatment of transparent surfaces in SOLPS-ITER interface
- Corrected charge assignment for bremsstrahlung calculation in EIRENE_SAMVOL
- Added array size check in EIRMOD_COUTAU
- Removed NLSRF condition to calls to EIRENE_UPSUSR and EIRENE_UPDATE_SPECTRUM in EIRMOD_LOCATE
- Bug fix for reading pressure loop feedback data
- Bug fix due to typo in SGNAL routine + associated correction to ITER CI serial test cases
- Bug fixes for species-specific rescaling scheme
- Corrected printout of ESHEATH in SOLPS-ITER cases
- Removed resetting of IFIRST to 0 after PREPARE step in EIRMOD_UPDLIN
- Added some missing USE_OPENMP pragmas
- Corrected pragmas in CMake version file
- Added definition of MPI_VERSION pragma in CMakeLists.txt
- Corrected handling of NPESTA and NPESTR arrays
- Added broadcasting of NPLRPI in EIRMOD_COMXS
- Use of NFSTPI instead of NFRSTP in EIRMOD_PLTEIR
- Added deallocation call of WNEUTRAL_FLUXES
- Corrected EIRDIAG_NDS size for SOLPS-ITER checking procedures
- Corrected procedure names in IOUSR interface files
- Removed a premature RETURN in WRPLAM_LONG
- Conditioned a gfortran compiler flag to minimal version
- Corrected some reading formats
- Corrected some comments
- Alignment corrections

---

## X. Bonnin - feature/SOLPS_push_clean-up - 781960eea7b56a886fccecd4ad98e4002fbc5292

#### Based on develop - 1138b60ed32d3337e392f716ece854bfc89675a3

** Changes **

This commit is meant to apply good coding practices with minimal consequences on actual results. It addresses warnings from many different compilers to make the Eirene compilation as clean as reasonably achievable.

- Applied CONTRIBUTING rule that Fortran file extensions should be .F instead of .f when the file must be pre-processed by the Makefile
- Pointing to modified set of examples with .F file extensions
- Removed forced usage of OpenMP compiler options
- Added GRAPHICS option in CMakeLists.txt for compilation with or without plotting enabled
- Added EIRENE_ALLOC_CGRID call in COUPLE_PARAM from SOLPS-ITER routine
- Used the module provided NREACI, NATMI, NMOLI, NIONI, and NPHOTI in EIRENE_BROWSE_BLOCK_5 inside FIND_PARAM_JSON
- Corrected initialization of WEISPZ in EIRMOD_LOCAT0
- Added ITYP_OLD assignments in EIRMOD_LOCAT1
- Added IFIRST initialization to zero in READ_COLRAD
- Adding COUPLE_INIT routine in B2.5 interface
- Moved EIRENE_EXPINT, EIRENE_LAX_M and EIRENE_MMDEI to their own files
- Apply separate scaling to additional, algebraic and pumping tallies
- Added TRCMELD output switch
- Refactored computation of random seeds to be compatible with both OpenMP and MPI parallelization schemes
- Do not use IUNOUT in EIRENE_MAIN before it is defined or after it has been deallocated
- Added MPI_VERSION and MPI_MOD pragmas to handle older MPI implementations
- Ensuring MPI_STATUSES has consistent size for MPI runs
- Added a safety and error message in FIND_PARAM if the input file is found empty
- Added a safety to check the fort.34 file is present before attempting to read it
- Added safety for reading correct number of integers in block 10F of input file in READ_FIXFORM
- Added safety against IPRNL=0 in SETUP_TIME_SURFACE
- Adding safety against TEIL=0 and an error message in EIRMOD_ALGEBRA
- Adding safety for incorrect algebraic manipulation of intensive tallies in ALGTAL
- Added calls to flush out buffers
- Added allocated status safeties
- Used NPLSI instead of NPLS in READ_FIXFORM
- Replaced NPHOT by NPHOTI in SCALE_TALLIES
- Making NDXD, etc... PUBLIC in EIRMOD_EIRBRA
- Removed NDXP, NDYP, NFL from EIRMOD_BRAEIR
- Removed MAXPOIN from EIRMOD_SOLPS
- Making LMETSP2 and LMETSPW2 threadprivate
- Rearranged OpenMP critical region in EIRMOD_UPTBGK
- Added OpenMP compliant error handling
- Added information about thread number in error state in EIRENE_EXIT_OWN
- Added census score output in TMSTEP
- Added SHSTEP output in EIRMOD_INFCOP from SOLPS-ITER interface
- Made FLUXES and FLUXS arrays in EIRSRT allocatable to avoid a segmentation fault
- Made DUMMY allocatable in INTEGRATE_TALLIES
- Made arrays in EIRMOD_PLTEIR allocatable
- Added SAVE statement in EIRMOD_LOCAT0
- Removed misleading superfluous output from CMakeLists.txt
- Added standard deviation output in OUTIDLTAL
- Removed redundant output in OUTIDLPLA
- Corrected sizes of tally output arrays in SETTXT_INTAL
- Conditioned some debugging output
- Rewrote some output formats to avoid creation of temporary arrays
- Clarified some code output
- Some aesthetic output changes
- Corrected format statements
- Refactored reading of fort.13 file
- Rewrote EIRENE_FTHOMP to avoid a fatal round-off error
- Use of EVKEL instead of approximate value
- Generalized usaged of FORT and FORT_LC variables
- Introduced usage of EIRENE_OPENFILE in all interface routines
- Avoiding obsolescent shared labels END DO statements, replaced with CONTINUE
- Added _DP type declarations for floating-point constants
- Added type assignments when using integers inside floating-point operations
- Added TRIM operations and character ranges to avoid string length overflows
- Added USE_OPENMP and USE_EXT_OPENMP pragmas to avoid compiler warnings from !$OMP instructions when compiling without OpenMP
- Removed double definition of USE_OPENMP and USE_EXT_OPENMP in CMakeLists.txt
- Added B25_EIRENE pre-processor pragmas in SOLPS-ITER interface routines
- Added LEGACYCOMP alternate code
- Added definition of GFORTRAN pragma in CMakeLists.txt
- Added G95 pre-processor pragmas to deal with g95 compiler errors
- Added NAGFOR pre-processor pragmas to deal with NAG compiler specific warnings
- Replaced F2003 pre-processor (the compiler supports the Fortran 2003 standard) pragma with LEGACYCOMP (the compiler does not support the Fortran 2003 standard)
- Removed unused IFOFF argument in EIRENE_FIND_TRIANG_DIM
- Removed usage of obsolete variables NITER0 in INPUT and NTIME0 in READ_FIXFORM
- Removed conflicting DIMENSION statements
- Avoided use of module variables as loop indices
- Added EXTERNAL definitions for called subroutines and functions
- Using RETURN and END SUBROUTINE <SUB_NAME> statements to close subroutines
- Moved statement functions to become contained subroutines
- Added IMPLICIT NONE statements
- Added some INTENT(IN) declarations
- Added pre-processor pragmas around relevant variables only used if the pragma is declared
- Moved declarations of loop index variables inside contained routines when appropriate
- Add IGNORE comments to MPI and JSON modules use statements for SOLPS-style dependencies build
- Removed superfluous MPI interface routines
- Removed unused variables and modules
- Removed unused labels, formats and error statements
- Making use of lower/upper case more consistent
- Removed tabulations, alignment corrections
- Added some more clarifying comments
- Removed trailing whitespaces and empty lines
- Have figures in Manual listed by full relative path and file extension
- Added missing comment about density validity range in AMJUEL file
- Added some more produced file types to .gitignore file
- Improved handling of external libraries in CMakeLists.txt
- Added compilation dependency on local environment setup files
- Added -fpe0 debug Intel compilation flag
- Added -allow-argument-mismatch argument to GNU compilation above version 9.5 to avoid errors with MPI routines
- Added Cray compiler support
- Implemented new logic in test-case scripts for input_fixed, input_json, output_fixed and output_json directories
- Adapted test-cases for new outtal_* file format.
- Adjusted openMP test-case Makefiles for easy CI data creations.

---

## D. Harting - New EIRENE version: Release-V1.1.0 - a339c4b70032409286e399b59f01e6bcb6940f06

New_EIRENE_Version=1.1.0

```
+===============================================+
|                                               |
|   New EIRENE version: Release-V1.1.0          |
|   Date: 31st October 2024                     |
|   MsV:  MsV-JuelichOrigin_1.0-Oct2024-V1.1.0  |
|                                               |
+===============================================+
```
**Changes**

1. Added new JSON input file format fro EIRENE.
   - See EIRENE manual for more details.
2. Added OpenMP capability for EIRENE.
3. Updates to EIRENE database
   - Changes to coefficients
      - amjuel.tex: Reactions 2.2.9 (e + H_2 -> 2e + H_2^+), 2.3.6A0 (C^+ + e -> C), 2.6A0 (e + C -> C^+ + 2e (ADAS 93 -> ADAS 96)), 2.8A0 (e + O -> O^+ + 2e)
      - hydhel.tex: Reaction 3.2.2 (p + H_2(v=0) -> p + H_2(v > 0))
      - methane.tex: Reaction 3.2 (p + C -> C^+ + H)
      - h2vibr.tex: Reaction 2.4l1T
   - Changes in several reactions to T1MIN, T1MAX, N2MIN, N2MAX and Eth
4. Unification of internal EIRENE versions.
5. Preparations for unification of EIRENE with SOLPS-ITER version.
6. Internal code restructuring with some minor impacts on coupling/user-routines.
7. Restructuring of GitLab CI.

---

## D. Harting - DMH/fix_ci_version_check - b6e4b2df6d71f4195faa1e23fa1dd79a0f56a5c8

#### Based on Vx.y.z(undefined before this commit - use instead develop - 1138b60ed32d3337e392f716ece854bfc89675a3)

**Changes**

- Stop the CI pipeline if the compliance_check failed
- Corrected logic for version check in compliance_check
   - compliance_check failed when e.g. minor version number decreased and patch number decreased. This is now fixed
   - Check not only the last commit if version number increased but check the last two commits which changed the first line of version.txt if the version number increased

## H. Leggate - HJL/add_changelog_and_coding_rules_version_check - 5cc860bae983754e7aaabb8797f8b340e80e2ac7

#### Based on Vx.y.z(undefined before this commit - use instead develop - 375c8f50e4e5641fed3654c9b59a41e429f590f9)

New_EIRENE_Version=1.0.10

**Changes**

This adds further checks on versions at various points in the code and adds a coding rules document

- Added coding_rules.md imported from ModCR
- Added check on Version in coding_rules.md to pre-commit git hook
- Added check for New_EIRENE_Version in CHANGELOG.md in ci_version_check.sh
- Added short guide to merging in CONTRIBUTING.md
- Changed OpenMP Sample Case version to add in 2 exclusions from 2D-poly

---

## D. Harting - DMH/streamline_ci - 54b400477165ea16b146ca28b0361cf4ff916a48

#### Based on on Vx.y.z (undefined before this commit - use instead develop - 375c8f50e4e5641fed3654c9b59a41e429f590f9)

**Changes**

- All the ITER test cases were building the whole EIRENE library again

  - Build now only once the ITER EIRENE executable at preparation stage to speedup the CI pipeline

  - Reuse the ITER EIRENE executable at run stage

  - No changes in referece outputs of ITER test cases


- Reduced number of MC particles in EMC3 test cases to speedup CI pipeline

  - Reference outputs of EMC3 test cases were updated due to new results with reduced number of MC particles


- Removed -fopenmp compile switch from non openmp test cases

  - Removed -fopenmp switch from CMakeList.txt of EIRENE repository and from the Makefiles in the sample case repository

  - Updated test case ITER_1573_scaling due to small changes from removal of -fopenmp switch

  - Removed ouput files fort.111 and fort.113 generated if CHECKBIN environment variable is set during compile time as these files are now inconsistent with updated test cases and not used during CI checks

---

## X. Bonnin - feature/SOLPS_push_whitespace-comments - 66c03ee2edeeb58def8b6a4ef7b6aefd02ae797b

#### Based on Vx.y.z(undefined before this commit - use instead develop - d3033c4626f220acfb060aecc8c29c523b8f4d4c)

**Changes**

This commit is meant is decrease the distance between the SOLPS Eirene branch and the MsV reference develop branch. It is not meant to contain any functional changes.

- Updated version number to 1.0.9.
- Streamlined detection of graphical library dependencies for compilation.
- Some formatting corrections to avoid stars in output.
- Added a more informative error message when EWALL=0 is set.
- Added informative message from CMake as to what compiler version is being used.
- Added flushing of buffers on error exit to make sure all output is printed.
- Bug fix to OUTAU array size to prevent a CI test failure.
- Corrected generalized interface UPTUSR routine to do the same as its counterparts.
- Added SAVE and THREADPRIVATE statements for local arrays in SUMOSTRA, EIRMOD_COUTAU and EIRMOD_LOCATE routines.
- Alignment bug fix in EIRMOD_COMXS to ensure RP%POLY is allocated.
- Keeping function and routine name declarations on a single line to facilitate code searches.
- Changed CI rules to always produce artifacts.
- Cleaned up debugging output.
- Corrected explanation of SRCML array (folllowing bug fix from Niels Horsten).
- Corrected default value of DTIMV in Manual.
- Added ISPCOPT variable description in Manual.
- Corrected typos.
- Added clarification comments (mostly from Detlev Reiter).
- Right-justifying line labels to satisfy the NAG compiler preprocessor.
- Some aethestic changes in the positioning of empty lines (in both the code source and the output).
- Some minor additional output to match that from the SOLPS branch.
- Removing trailing whitespaces as per the CONTRIBUTING rules.
- Removing tabs as indentation characters.
- Removed some superfluous commented out code lines.
- Removing some superfluous parentheses.
- Truncating long comment lines to avoid compiler warnings.
- Enforced proper indentation when found misaligned.
- Added some spaces to improve code readability.
- Changes of character case to match the SOLPS branch.
- Removed duplicate prb89.dat file.
- Added local user source files to .gitignore list.

---

## D. Harting - DMH/develop_fix_compare_scripts - 0d9077f21ce0c99835d44719821f7e00ad5073f3

#### Based on Vx.y.z (undefined before this commit - use instead develop - d3033c4626f220acfb060aecc8c29c523b8f4d4c)

**Changes**

- Only changes to CI logic, no code changes

  - Removed fort.37 from reference outputs of ITER cases. Output to fort.37 was abandoned and is now included in the standard output of EIRENE

- Removed dependency of OPENMP cases on previous cases

  - OPENMP library and test case preparation downloaded all previous artifacts. This is not necessary as the OPENMP cases should be fresh compilations

- Removed some cells in the 2D-D_polygon_openmp  test cases from the comparison as they contain small numbers and relative change is large due to statistical noise

---

## D. Harting - DMH/develop_fix_compare_scripts - 7cf909853e2065c8213459946d13e63a9e9a427d

#### Based on Vx.y.z (undefined before this commit - use instead develop - d3033c4626f220acfb060aecc8c29c523b8f4d4c)

**Description**

- CI checks (compare scrips) were switched off and CI signaled 'green' even though the output has changed ([issue #82](https://jugit.fz-juelich.de/eirene/eirene/-/issues/82))
- The checks were switched off by commits 40a86ae38b0eeb860f7f4553a15ade890749f3f3 (2021-10-14) and 28ba4aa745452a67b75d6836c1e60aefa061af20 (2022-12-14)
- Fixed compare scripts to flag differences to reference output
- Updated standard test cases
- **IMPORTANT!!! NEW TEST CASES WERE NOT YET VALIDATED !!!**

**Changes in EIRENE Database since compare error was introduced in  40a86ae38b0eeb860f7f4553a15ade890749f3f3**

- amjuel.tex

  - Changes in coefficients:
    - Amjuel reaction 2.2.9  e + H_2  -> 2e + H_2^+
    - Amjuel Reaction 2.3.6A0  C^+ + e -> C
    - Amjuel Reaction 2.6A0  e + C  -> C^+   + 2e (ADAS 93 -> ADAS 96)
    - Amjuel Reaction 2.8A0  e + O  -> O^+   + 2e
  - Changes to T1MIN, T1MAX, N2MIN, N2MAX
    - Amjuel Reaction 2.7A0r  e + N  -> N^+   + 2e
    - Amjuel Reaction 2.6A0 C  + e -> C^+   + 2e
  - Changes to N2MIN, N2MAX
    - Amjuel reaction 2.1.5JH  e + H -> H^+ + 2e
    - Amjuel reaction 2.1.5o  e + H -> H^+ + 2e   Ly-opaque
    - Amjuel Reaction 2.1.8o H^+ + e -> H(1s)  Ly-opaque
    - Amjuel Reaction 2.2.h2c H_2 + e ->  ....
    - Amjuel Reaction 2.2.h2r  H_2 + e ->  ....
    - Amjuel Reaction 2.3.9a  e + He(1s^21S) -> He^+(1s) + 2e
    - Amjuel Reaction 2.3.13a  e+ He^+(1s) ->   He(1s^21S)
  - Changes to T1MIN
    - Amjuel reaction 7.2.3a    p + H^{-} ->  H + H (for cold H^-)
    - Amjuel reaction 7.2.3b   p + H^{-} ->  H + H^+ + 2e (for cold H^-)
  - Changes to P2MIN, P2MAX
    - Amjuel Reaction 2.1.5   H + e -> H^+ +2e,  Ratio H^+/H(1)
    - Amjuel Reaction 2.1.5a  H + e -> H^+ + 2e , Ratio H(3)/H(1)
    - Amjuel Reaction Reaction 2.1.5b  H + e -> H^+ + 2e, Ratio H(2)/H(1)
    - Amjuel eaction 2.1.5c  H + e -> H^+ + 2e, Ratio H(4)/H(1)
    - Amjuel Reaction 2.1.5d  H + e -> H^+ + 2e, Ratio H(5)/H(1)
    - Amjuel Reaction 2.1.5e  H + e -> H^+ + 2e, Ratio H(6)/H(1)
    - Amjuel Reaction 2.1.5tot  H + e -> H^+ + 2e, Ratio H(tot)/H(1)
    - Amjuel Reaction 2.1.5de   H + e -> H^+ + 2e
    - Amjuel Reaction 2.1.5o    H + e -> H^+ + 2e  Ly-opaque
    - Amjuel Reaction 2.1.8   H^+ + e -> H(1s), Ratio $H(1)/H^+
    - Amjuel Reaction 2.1.8a  H^+ + e -> H(1s), Ratio H(3)/H^+
    - Amjuel Reaction 2.1.8b  H^+ + e -> H(1s) ,  Ratio H(2)/H^+
    - Amjuel Reaction 2.1.8c  H^+ + e -> H(1s) ,  Ratio H(4)/H^+
    - Amjuel Reaction 2.1.8d  H^+ + e -> H(1s) ,  Ratio H(5)/H^+
    - Amjuel Reaction 2.1.8e  H^+ + e -> H(1s) ,  Ratio H(6)/H^+
    - Amjuel Reaction 2.1.8tot  H^+ + e -> H(1s) , Ratio H(tot)/H^+
    - Amjuel Reaction 2.1.8de H^+  + e -> H(1s) , + 13.6 eV
    - Amjuel Reaction 2.1.8o H^+ + e -> H(1s) , 13.6 eV
    - Plus others......

- h2vibr.tex

  - H.1 :  Fits for sigma(E)
    - Nearly all Eth changed
  - Changes in coefficients
    - Reaction 2.4l1T
  - A lot of reactions changed T1MIN T1MAX
  - New Reactions added

- hydhel.tex

  - Changes of Eth
    - Reaction 3.2.2    p + H_2(v=0) -> p + H_2(v > 0)

- methane.tex

  - Changes in coefficients
    - Reaction 3.2      p + C -> C^+ + H

**Differences of new test cases generated**

- fort.37 generated by eirmod.infcop was abandoned (ecdedfcbb9d2f9a105462ead067205e38b5bc421).

  - The output to fort.37 was moved to the general EIRENE output file (IUNOUT)
  - Removed fort.37 from reference cases

- Differences in standard output file of EIRENE

  - The maximum number of ATOMIC and molecular reaction has increased by one (NREAC)
    - Affected test cases: 1D-H_slab, 1D_cylindric, 1D-elliptic, 2D-D_slab, 2D-D_slab_finite_cylinder_11x11,  2D-D_triang, 2D_cylindric, 2D_elliptic, 3D_slab_11x11x11,  ALTS2,  ALTS2_bits, cylinder, cylinder_triang_no_octree,  cylinder_triang_octree
    - Test cases not affected:  2D-D_polygon, 2D-D_polygon_NLPLG
  - Some test cases are signalling IEEE_DIVIDE_BY_ZERO
    - 2D-D_polygon, 2D-D_polygon_NLPLG
  - Some test cases are **NOT** signalling any more IEEE_INVALID_FLAG
    - cylinder
  - CREACD changed significantly
    - 2D-D_polygon,  2D-D_slab,  2D-D_slab_finite_cylinder_11x11,  2D-D_triang,  3D_slab_11x11x11,  3D-D_tetra,  3D-D_tetra
  - TABEI1 (AVR) changed significantly
    - 2D-D_polygon
  - EELEI1 (AVR) changes slightly
    - 2D-D_polygon_NLPLG
  - More warnings from SLREAC extrapolation
    - 2D-D_slab_finite_cylinder_11x11, 2D-D_triang, 3D_slab_11x11x11,  3D-D_tetra,

- **All EMC3 test cases are identical**

  - They have limited precision in their ASCII outputs

- Even though these testcases show differences, their global balances in standard output of EIRENE show no differences -> **OK to update these test cases**

  - 1D_cylindric, 1D_elliptic, 1D-H_slab, 2D_cylindric, cylinder, ITER_1573, ITER_1573_def_lines, ITER_1573_dens_model, ITER_1573_scaling, ALT2S, 2D-D_polygon_NLPLG, 2D_elliptic, ALT2S_bits

  - ITER test cases show differences due to neutral-neutral collisions in the IN-TALLIES

- These test cases show differences with the new EIRENE database, but they show no differences in their global balances with the old EIRENE database -> **OK to update test cases**

  - 2D-D_slab_finite_cylinder_11x11, 3D-D_tetra, 2D-D_polygon, 3D-D_slab_11x11x11

- These test cases show differences with the new and the old EIRENE Database. **Reasson for changes needs still to be undestood, not OK to update**

  - 2D-D_triang, cylinder_triang_no_octree, cylinder_triang_octree

  - 2D-D_slab: Shows only some small differences in global balances with old database !!!

---
