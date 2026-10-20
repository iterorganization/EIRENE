cdr  activate laguerre test (temporarily)
cdr  logicals for lralch, lvainst, ...
cdr  collisional_bray --> collisional_vainst, but what is collisional2?
cdr  return atomic structure: aik, energy levels, stat. weights.

      subroutine eirene_he_colrad (te, ne, bfield, iform,
     .                             pop_esc, q_ext, l_ext,
     .                             aikeinst,energlev,statwght)

cdr  derived and adapted from original "Goto code":
cdr  M. Goto,
cdr  Collisional-radiative model for neutral helium in plasma revisited,
cdr  Journal of Quantitative Spectroscopy and Radiative Transfer, 76(3-4), 331 (2003)

cdr  Further adaptations:
cdr  W. Zholobenko, M. Rack, D. Reiter, et al., Synthetic helium beam diagnostic and underlying atomic data
cdr  Journal Nuclear Fusion 58 (2018) 126006
cdr  W. Zholobenko et al., Report FZ Juelich, JUEL-4407, Feb. 2018, ISSN 0944-2952

cdr  translated back from C++ to Fortran: P.Boerner, FZ Juelich, 2018
cdr  final adaptations for direct use inside eirene: D.Reiter, 2018-2021
cdr  2021: comments, sanity checks, and form 1 vs. form 2 structures.
cdr        pop-coeffs. made dimensionless (valid also for corona ne=0.0 limit)

      use eirmod_precision
      use eirmod_hecr
      USE EIRMOD_COMPRT, ONLY: IUNOUT

      implicit none
      real*8, intent(in) :: te, ne, bfield
      integer, intent(in) :: iform
      real*8, intent(in) :: pop_esc(ULH,ULH), q_ext(ULH,ULH)
      logical, intent(in) :: l_ext

      real*8, intent(out) :: aikeinst(ULH,ULH), energlev(ULH),
     .                        statwght(ULH)

      aikeinst = 0.
      energlev = 0.
      statwght = 0.

      return
      end subroutine eirene_he_colrad

