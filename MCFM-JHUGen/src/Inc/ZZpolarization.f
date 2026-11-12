c--- Helicity selection for the two decaying Z/gamma* bosons in the
c--- WBFZZ (qq_ZZqq) process; see getVVZZamps_pol.f.
c---
c--- ZZpol_mode  = 0 : off (default). Original amplitudes, unchanged.
c---             = 1 : on. Only diagrams in which each lepton pair
c---                   comes from its own Z/gamma* propagator are kept,
c---                   and each pair is projected onto the helicities
c---                   chosen in ZZpol_hel.
c--- ZZpol_hel(1) selects the (3,4) pair, ZZpol_hel(2) the (5,6) pair:
c---             = 0 : all helicities (no projection)
c---             = 1 : longitudinal, lambda = 0
c---             = 2 : transverse, lambda = +1 and -1 (coherent sum)
c---             = 3 : lambda = +1 only
c---             = 4 : lambda = -1 only
c---             = 5 : all three via the projection (check; same as 0)
c--- ZZpol_frame = 0 : helicity axes defined in the frame of the
c---                   momenta given to MCFM
c---             = 1 : four-lepton (ZZ) rest frame
c---             = 2 : partonic centre-of-mass frame, -(p1+p2)
      integer ZZpol_mode,ZZpol_hel(2),ZZpol_frame
      common/ZZpolarization/ZZpol_mode,ZZpol_hel,ZZpol_frame
c--- internal flag, set only by getVVZZamps_pol while it computes
c--- polarized amplitudes; switches off the non-resonant diagrams
      logical ZZpol_resonly
      common/ZZpolinternal/ZZpol_resonly
!$omp threadprivate(/ZZpolinternal/)
