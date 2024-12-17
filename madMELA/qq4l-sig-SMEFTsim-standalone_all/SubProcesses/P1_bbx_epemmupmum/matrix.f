      SUBROUTINE M27_SMATRIXHEL(P,HEL,ANS)
      IMPLICIT NONE
C     
C     CONSTANT
C     
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=6)
      INTEGER                 NCOMB
      PARAMETER (             NCOMB=64)
CF2PY INTENT(OUT) :: ANS
CF2PY INTENT(IN) :: HEL
CF2PY INTENT(IN) :: P(0:3,NEXTERNAL)

C     
C     ARGUMENTS 
C     
      REAL*8 P(0:3,NEXTERNAL),ANS
      INTEGER HEL
C     
C     GLOBAL VARIABLES
C     
      INTEGER USERHEL
      COMMON/M27_HELUSERCHOICE/USERHEL
C     ----------
C     BEGIN CODE
C     ----------
      USERHEL=HEL
      CALL M27_SMATRIX(P,ANS)
      USERHEL=-1

      END

      SUBROUTINE M27_SMATRIX(P,ANS)
C     
C     Generated by MadGraph5_aMC@NLO v. 2.9.21, 2024-09-26
C     By the MadGraph5_aMC@NLO Development Team
C     Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
C     
C     MadGraph5_aMC@NLO StandAlone Version
C     
C     Returns amplitude squared summed/avg over colors
C     and helicities
C     for the point in phase space P(0:3,NEXTERNAL)
C     
C     Process: b b~ > e+ e- mu+ mu- NP<=1 NPprop=0 SMHLOOP=0 @1
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=6)
      INTEGER    NINITIAL
      PARAMETER (NINITIAL=2)
      INTEGER NPOLENTRIES
      PARAMETER (NPOLENTRIES=(NEXTERNAL+1)*6)
      INTEGER                 NCOMB
      PARAMETER (             NCOMB=64)
      INTEGER HELAVGFACTOR
      PARAMETER (HELAVGFACTOR=4)
C     
C     ARGUMENTS 
C     
      REAL*8 P(0:3,NEXTERNAL),ANS
CF2PY INTENT(OUT) :: ANS
CF2PY INTENT(IN) :: P(0:3,NEXTERNAL)
C     
C     LOCAL VARIABLES 
C     
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
C     put in common block to expose this variable to python interface
      COMMON/M27_PROCESS_NHEL/NHEL
      REAL*8 T
      REAL*8 M27_MATRIX
      INTEGER IHEL,IDEN, I, J
C     For a 1>N process, them BEAMTWO_HELAVGFACTOR would be set to 1.
      INTEGER BEAMS_HELAVGFACTOR(2)
      DATA (BEAMS_HELAVGFACTOR(I),I=1,2)/2,2/
      INTEGER JC(NEXTERNAL)
      LOGICAL GOODHEL(NCOMB)
      DATA NTRY/0/
      DATA GOODHEL/NCOMB*.FALSE./

C     
C     GLOBAL VARIABLES
C     
      INTEGER USERHEL
      COMMON/M27_HELUSERCHOICE/USERHEL
      DATA USERHEL/-1/
      LOGICAL HELRESET
      COMMON/M27_HELRESET/HELRESET
      DATA HELRESET/.TRUE./

      DATA (NHEL(I,   1),I=1,6) / 1,-1, 1,-1, 1,-1/
      DATA (NHEL(I,   2),I=1,6) / 1,-1, 1,-1, 1, 1/
      DATA (NHEL(I,   3),I=1,6) / 1,-1, 1,-1,-1,-1/
      DATA (NHEL(I,   4),I=1,6) / 1,-1, 1,-1,-1, 1/
      DATA (NHEL(I,   5),I=1,6) / 1,-1, 1, 1, 1,-1/
      DATA (NHEL(I,   6),I=1,6) / 1,-1, 1, 1, 1, 1/
      DATA (NHEL(I,   7),I=1,6) / 1,-1, 1, 1,-1,-1/
      DATA (NHEL(I,   8),I=1,6) / 1,-1, 1, 1,-1, 1/
      DATA (NHEL(I,   9),I=1,6) / 1,-1,-1,-1, 1,-1/
      DATA (NHEL(I,  10),I=1,6) / 1,-1,-1,-1, 1, 1/
      DATA (NHEL(I,  11),I=1,6) / 1,-1,-1,-1,-1,-1/
      DATA (NHEL(I,  12),I=1,6) / 1,-1,-1,-1,-1, 1/
      DATA (NHEL(I,  13),I=1,6) / 1,-1,-1, 1, 1,-1/
      DATA (NHEL(I,  14),I=1,6) / 1,-1,-1, 1, 1, 1/
      DATA (NHEL(I,  15),I=1,6) / 1,-1,-1, 1,-1,-1/
      DATA (NHEL(I,  16),I=1,6) / 1,-1,-1, 1,-1, 1/
      DATA (NHEL(I,  17),I=1,6) / 1, 1, 1,-1, 1,-1/
      DATA (NHEL(I,  18),I=1,6) / 1, 1, 1,-1, 1, 1/
      DATA (NHEL(I,  19),I=1,6) / 1, 1, 1,-1,-1,-1/
      DATA (NHEL(I,  20),I=1,6) / 1, 1, 1,-1,-1, 1/
      DATA (NHEL(I,  21),I=1,6) / 1, 1, 1, 1, 1,-1/
      DATA (NHEL(I,  22),I=1,6) / 1, 1, 1, 1, 1, 1/
      DATA (NHEL(I,  23),I=1,6) / 1, 1, 1, 1,-1,-1/
      DATA (NHEL(I,  24),I=1,6) / 1, 1, 1, 1,-1, 1/
      DATA (NHEL(I,  25),I=1,6) / 1, 1,-1,-1, 1,-1/
      DATA (NHEL(I,  26),I=1,6) / 1, 1,-1,-1, 1, 1/
      DATA (NHEL(I,  27),I=1,6) / 1, 1,-1,-1,-1,-1/
      DATA (NHEL(I,  28),I=1,6) / 1, 1,-1,-1,-1, 1/
      DATA (NHEL(I,  29),I=1,6) / 1, 1,-1, 1, 1,-1/
      DATA (NHEL(I,  30),I=1,6) / 1, 1,-1, 1, 1, 1/
      DATA (NHEL(I,  31),I=1,6) / 1, 1,-1, 1,-1,-1/
      DATA (NHEL(I,  32),I=1,6) / 1, 1,-1, 1,-1, 1/
      DATA (NHEL(I,  33),I=1,6) /-1,-1, 1,-1, 1,-1/
      DATA (NHEL(I,  34),I=1,6) /-1,-1, 1,-1, 1, 1/
      DATA (NHEL(I,  35),I=1,6) /-1,-1, 1,-1,-1,-1/
      DATA (NHEL(I,  36),I=1,6) /-1,-1, 1,-1,-1, 1/
      DATA (NHEL(I,  37),I=1,6) /-1,-1, 1, 1, 1,-1/
      DATA (NHEL(I,  38),I=1,6) /-1,-1, 1, 1, 1, 1/
      DATA (NHEL(I,  39),I=1,6) /-1,-1, 1, 1,-1,-1/
      DATA (NHEL(I,  40),I=1,6) /-1,-1, 1, 1,-1, 1/
      DATA (NHEL(I,  41),I=1,6) /-1,-1,-1,-1, 1,-1/
      DATA (NHEL(I,  42),I=1,6) /-1,-1,-1,-1, 1, 1/
      DATA (NHEL(I,  43),I=1,6) /-1,-1,-1,-1,-1,-1/
      DATA (NHEL(I,  44),I=1,6) /-1,-1,-1,-1,-1, 1/
      DATA (NHEL(I,  45),I=1,6) /-1,-1,-1, 1, 1,-1/
      DATA (NHEL(I,  46),I=1,6) /-1,-1,-1, 1, 1, 1/
      DATA (NHEL(I,  47),I=1,6) /-1,-1,-1, 1,-1,-1/
      DATA (NHEL(I,  48),I=1,6) /-1,-1,-1, 1,-1, 1/
      DATA (NHEL(I,  49),I=1,6) /-1, 1, 1,-1, 1,-1/
      DATA (NHEL(I,  50),I=1,6) /-1, 1, 1,-1, 1, 1/
      DATA (NHEL(I,  51),I=1,6) /-1, 1, 1,-1,-1,-1/
      DATA (NHEL(I,  52),I=1,6) /-1, 1, 1,-1,-1, 1/
      DATA (NHEL(I,  53),I=1,6) /-1, 1, 1, 1, 1,-1/
      DATA (NHEL(I,  54),I=1,6) /-1, 1, 1, 1, 1, 1/
      DATA (NHEL(I,  55),I=1,6) /-1, 1, 1, 1,-1,-1/
      DATA (NHEL(I,  56),I=1,6) /-1, 1, 1, 1,-1, 1/
      DATA (NHEL(I,  57),I=1,6) /-1, 1,-1,-1, 1,-1/
      DATA (NHEL(I,  58),I=1,6) /-1, 1,-1,-1, 1, 1/
      DATA (NHEL(I,  59),I=1,6) /-1, 1,-1,-1,-1,-1/
      DATA (NHEL(I,  60),I=1,6) /-1, 1,-1,-1,-1, 1/
      DATA (NHEL(I,  61),I=1,6) /-1, 1,-1, 1, 1,-1/
      DATA (NHEL(I,  62),I=1,6) /-1, 1,-1, 1, 1, 1/
      DATA (NHEL(I,  63),I=1,6) /-1, 1,-1, 1,-1,-1/
      DATA (NHEL(I,  64),I=1,6) /-1, 1,-1, 1,-1, 1/
      DATA IDEN/36/

      INTEGER POLARIZATIONS(0:NEXTERNAL,0:5)
      COMMON/M27_BORN_BEAM_POL/POLARIZATIONS
      DATA ((POLARIZATIONS(I,J),I=0,NEXTERNAL),J=0,5)/NPOLENTRIES*-1/

C     
C     FUNCTIONS
C     
      LOGICAL M27_IS_BORN_HEL_SELECTED

C     ----------
C     Check if helreset mode is on
C     ---------
      IF (HELRESET) THEN
        NTRY = 0
        DO I=1,NCOMB
          GOODHEL(I) = .FALSE.
        ENDDO
        HELRESET = .FALSE.
      ENDIF

C     ----------
C     BEGIN CODE
C     ----------
      IF(USERHEL.EQ.-1) NTRY=NTRY+1
      DO IHEL=1,NEXTERNAL
        JC(IHEL) = +1
      ENDDO
C     When spin-2 particles are involved, the Helicity filtering is
C      dangerous for the 2->1 topology.
C     This is because depending on the MC setup the initial PS points
C      have back-to-back initial states
C     for which some of the spin-2 helicity configurations are zero.
C      But they are no longer zero
C     if the point is boosted on the z-axis. Remember that HELAS
C      helicity amplitudes are no longer
C     lorentz invariant with expternal spin-2 particles (only the
C      helicity sum is).
C     For this reason, we simply remove the filterin when there is
C      only three external particles.
      IF (NEXTERNAL.LE.3) THEN
        DO IHEL=1,NCOMB
          GOODHEL(IHEL)=.TRUE.
        ENDDO
      ENDIF
      ANS = 0D0
      DO IHEL=1,NCOMB
        IF (USERHEL.EQ.-1.OR.USERHEL.EQ.IHEL) THEN
          IF (GOODHEL(IHEL) .OR. NTRY .LT. 20.OR.USERHEL.NE.-1) THEN
            IF(NTRY.GE.2.AND.POLARIZATIONS(0,0).NE.
     $       -1.AND.(.NOT.M27_IS_BORN_HEL_SELECTED(IHEL))) THEN
              CYCLE
            ENDIF
            T=M27_MATRIX(P ,NHEL(1,IHEL),JC(1))
            IF(POLARIZATIONS(0,0).EQ.
     $       -1.OR.M27_IS_BORN_HEL_SELECTED(IHEL)) THEN
              ANS=ANS+T
            ENDIF
            IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL)) THEN
              GOODHEL(IHEL)=.TRUE.
            ENDIF
          ENDIF
        ENDIF
      ENDDO
      ANS=ANS/DBLE(IDEN)
      IF(USERHEL.NE.-1) THEN
        ANS=ANS*HELAVGFACTOR
      ELSE
        DO J=1,NINITIAL
          IF (POLARIZATIONS(J,0).NE.-1) THEN
            ANS=ANS*BEAMS_HELAVGFACTOR(J)
            ANS=ANS/POLARIZATIONS(J,0)
          ENDIF
        ENDDO
      ENDIF
      END


      REAL*8 FUNCTION M27_MATRIX(P,NHEL,IC)
C     
C     Generated by MadGraph5_aMC@NLO v. 2.9.21, 2024-09-26
C     By the MadGraph5_aMC@NLO Development Team
C     Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
C     
C     Returns amplitude squared -- no average over initial
C      state/symmetry factor
C     for the point with external lines W(0:6,NEXTERNAL)
C     
C     Process: b b~ > e+ e- mu+ mu- NP<=1 NPprop=0 SMHLOOP=0 @1
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=60)
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=6)
      INTEGER    NWAVEFUNCS, NCOLOR
      PARAMETER (NWAVEFUNCS=11, NCOLOR=1)
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
      COMPLEX*16 IMAG1
      PARAMETER (IMAG1=(0D0,1D0))
C     
C     ARGUMENTS 
C     
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C     
C     LOCAL VARIABLES 
C     
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR), TMP_JAMP(0)
      COMPLEX*16 W(20,NWAVEFUNCS)
      COMPLEX*16 DUM0,DUM1
      DATA DUM0, DUM1/(0D0, 0D0), (1D0, 0D0)/
C     
C     GLOBAL VARIABLES
C     
      INCLUDE 'coupl.inc'

C     
C     COLOR DATA
C     
      DATA (CF(I,  1),I=  1,  1) /3.000000000000000D+00/
C     1 T(2,1)
C     ----------
C     BEGIN CODE
C     ----------
      CALL IXXXXX(P(0,1),MDL_MB,NHEL(1),+1*IC(1),W(1,1))
      CALL OXXXXX(P(0,2),MDL_MB,NHEL(2),-1*IC(2),W(1,2))
      CALL IXXXXX(P(0,3),MDL_ME,NHEL(3),-1*IC(3),W(1,3))
      CALL OXXXXX(P(0,4),MDL_ME,NHEL(4),+1*IC(4),W(1,4))
      CALL IXXXXX(P(0,5),MDL_MMU,NHEL(5),-1*IC(5),W(1,5))
      CALL OXXXXX(P(0,6),MDL_MMU,NHEL(6),+1*IC(6),W(1,6))
      CALL FFS1_3(W(1,1),W(1,2),GC_520,MDL_MH,MDL_WH,W(1,7))
      CALL FFV1_3_3(W(1,3),W(1,4),GC_197,GC_136,MDL_MZ,MDL_WZ,W(1,8))
      CALL FFV1_3_3(W(1,5),W(1,6),GC_197,GC_136,MDL_MZ,MDL_WZ,W(1,9))
C     Amplitude(s) for diagram number 1
      CALL VVS4_0(W(1,8),W(1,9),W(1,7),GC_285,AMP(1))
      CALL FFS2_3(W(1,1),W(1,2),GC_468,MDL_MH,MDL_WH,W(1,7))
      CALL FFV1P0_3(W(1,3),W(1,4),GC_4,ZERO,ZERO,W(1,10))
      CALL FFV1P0_3(W(1,5),W(1,6),GC_4,ZERO,ZERO,W(1,11))
C     Amplitude(s) for diagram number 2
      CALL VVS3_0(W(1,10),W(1,11),W(1,7),GC_270,AMP(2))
C     Amplitude(s) for diagram number 3
      CALL VVS3_0(W(1,10),W(1,11),W(1,7),GC_325,AMP(3))
C     Amplitude(s) for diagram number 4
      CALL VVS6_0(W(1,10),W(1,11),W(1,7),GC_269,AMP(4))
C     Amplitude(s) for diagram number 5
      CALL VVS3_0(W(1,10),W(1,11),W(1,7),GC_330,AMP(5))
C     Amplitude(s) for diagram number 6
      CALL VVS6_0(W(1,10),W(1,11),W(1,7),GC_322,AMP(6))
C     Amplitude(s) for diagram number 7
      CALL VVS6_0(W(1,10),W(1,11),W(1,7),GC_329,AMP(7))
C     Amplitude(s) for diagram number 8
      CALL VVS2_0(W(1,10),W(1,9),W(1,7),GC_437,AMP(8))
C     Amplitude(s) for diagram number 9
      CALL VVS2_0(W(1,10),W(1,9),W(1,7),GC_320,AMP(9))
C     Amplitude(s) for diagram number 10
      CALL VVS6_0(W(1,10),W(1,9),W(1,7),GC_436,AMP(10))
C     Amplitude(s) for diagram number 11
      CALL VVS2_0(W(1,10),W(1,9),W(1,7),GC_326,AMP(11))
C     Amplitude(s) for diagram number 12
      CALL VVS6_0(W(1,10),W(1,9),W(1,7),GC_319,AMP(12))
C     Amplitude(s) for diagram number 13
      CALL VVS6_0(W(1,10),W(1,9),W(1,7),GC_321,AMP(13))
C     Amplitude(s) for diagram number 14
      CALL VVS2_0(W(1,11),W(1,8),W(1,7),GC_437,AMP(14))
C     Amplitude(s) for diagram number 15
      CALL VVS2_0(W(1,11),W(1,8),W(1,7),GC_320,AMP(15))
C     Amplitude(s) for diagram number 16
      CALL VVS6_0(W(1,11),W(1,8),W(1,7),GC_436,AMP(16))
C     Amplitude(s) for diagram number 17
      CALL VVS2_0(W(1,11),W(1,8),W(1,7),GC_326,AMP(17))
C     Amplitude(s) for diagram number 18
      CALL VVS6_0(W(1,11),W(1,8),W(1,7),GC_319,AMP(18))
C     Amplitude(s) for diagram number 19
      CALL VVS6_0(W(1,11),W(1,8),W(1,7),GC_321,AMP(19))
C     Amplitude(s) for diagram number 20
      CALL VVS3_0(W(1,8),W(1,9),W(1,7),GC_272,AMP(20))
C     Amplitude(s) for diagram number 21
      CALL VVS6_0(W(1,8),W(1,9),W(1,7),GC_271,AMP(21))
C     Amplitude(s) for diagram number 22
      CALL VVS1_5_0(W(1,8),W(1,9),W(1,7),GC_323,GC_438,AMP(22))
C     Amplitude(s) for diagram number 23
      CALL VVS4_0(W(1,8),W(1,9),W(1,7),GC_285,AMP(23))
C     Amplitude(s) for diagram number 24
      CALL VVS3_0(W(1,8),W(1,9),W(1,7),GC_324,AMP(24))
C     Amplitude(s) for diagram number 25
      CALL VVS6_0(W(1,8),W(1,9),W(1,7),GC_327,AMP(25))
C     Amplitude(s) for diagram number 26
      CALL VVS4_0(W(1,8),W(1,9),W(1,7),GC_435,AMP(26))
C     Amplitude(s) for diagram number 27
      CALL VVS3_0(W(1,8),W(1,9),W(1,7),GC_328,AMP(27))
C     Amplitude(s) for diagram number 28
      CALL VVS4_0(W(1,8),W(1,9),W(1,7),GC_464,AMP(28))
C     Amplitude(s) for diagram number 29
      CALL VVS4_0(W(1,8),W(1,9),W(1,7),GC_465,AMP(29))
C     Amplitude(s) for diagram number 30
      CALL VVS4_0(W(1,8),W(1,9),W(1,7),GC_466,AMP(30))
C     Amplitude(s) for diagram number 31
      CALL VVS4_0(W(1,8),W(1,9),W(1,7),GC_467,AMP(31))
      CALL FFV5_3(W(1,5),W(1,6),GC_352,MDL_MZ,MDL_WZ,W(1,11))
C     Amplitude(s) for diagram number 32
      CALL VVS4_0(W(1,8),W(1,11),W(1,7),GC_285,AMP(32))
      CALL FFV2_3(W(1,5),W(1,6),GC_803,MDL_MZ,MDL_WZ,W(1,11))
C     Amplitude(s) for diagram number 33
      CALL VVS4_0(W(1,8),W(1,11),W(1,7),GC_285,AMP(33))
      CALL FFV8_3(W(1,5),W(1,6),GC_804,MDL_MZ,MDL_WZ,W(1,11))
C     Amplitude(s) for diagram number 34
      CALL VVS4_0(W(1,8),W(1,11),W(1,7),GC_285,AMP(34))
      CALL FFV1_3(W(1,5),W(1,6),GC_364,MDL_MZ,MDL_WZ,W(1,11))
C     Amplitude(s) for diagram number 35
      CALL VVS4_0(W(1,8),W(1,11),W(1,7),GC_285,AMP(35))
      CALL FFV3_3(W(1,5),W(1,6),GC_353,MDL_MZ,MDL_WZ,W(1,11))
C     Amplitude(s) for diagram number 36
      CALL VVS4_0(W(1,8),W(1,11),W(1,7),GC_285,AMP(36))
      CALL FFV5_3(W(1,5),W(1,6),GC_392,MDL_MZ,MDL_WZ,W(1,11))
C     Amplitude(s) for diagram number 37
      CALL VVS4_0(W(1,8),W(1,11),W(1,7),GC_285,AMP(37))
      CALL FFV2_3(W(1,5),W(1,6),GC_813,MDL_MZ,MDL_WZ,W(1,11))
C     Amplitude(s) for diagram number 38
      CALL VVS4_0(W(1,8),W(1,11),W(1,7),GC_285,AMP(38))
      CALL FFV8_3(W(1,5),W(1,6),GC_814,MDL_MZ,MDL_WZ,W(1,11))
C     Amplitude(s) for diagram number 39
      CALL VVS4_0(W(1,8),W(1,11),W(1,7),GC_285,AMP(39))
      CALL FFV1_3(W(1,5),W(1,6),GC_395,MDL_MZ,MDL_WZ,W(1,11))
C     Amplitude(s) for diagram number 40
      CALL VVS4_0(W(1,8),W(1,11),W(1,7),GC_285,AMP(40))
      CALL FFV3_3(W(1,5),W(1,6),GC_384,MDL_MZ,MDL_WZ,W(1,11))
C     Amplitude(s) for diagram number 41
      CALL VVS4_0(W(1,8),W(1,11),W(1,7),GC_285,AMP(41))
      CALL FFV5_3(W(1,5),W(1,6),GC_398,MDL_MZ,MDL_WZ,W(1,11))
C     Amplitude(s) for diagram number 42
      CALL VVS4_0(W(1,8),W(1,11),W(1,7),GC_285,AMP(42))
      CALL FFV3_3(W(1,5),W(1,6),GC_387,MDL_MZ,MDL_WZ,W(1,11))
C     Amplitude(s) for diagram number 43
      CALL VVS4_0(W(1,8),W(1,11),W(1,7),GC_285,AMP(43))
      CALL FFV5_3(W(1,3),W(1,4),GC_352,MDL_MZ,MDL_WZ,W(1,11))
C     Amplitude(s) for diagram number 44
      CALL VVS4_0(W(1,11),W(1,9),W(1,7),GC_285,AMP(44))
      CALL FFV8_3(W(1,3),W(1,4),GC_709,MDL_MZ,MDL_WZ,W(1,11))
C     Amplitude(s) for diagram number 45
      CALL VVS4_0(W(1,11),W(1,9),W(1,7),GC_285,AMP(45))
      CALL FFV2_3(W(1,3),W(1,4),GC_708,MDL_MZ,MDL_WZ,W(1,11))
C     Amplitude(s) for diagram number 46
      CALL VVS4_0(W(1,11),W(1,9),W(1,7),GC_285,AMP(46))
      CALL FFV1_3(W(1,3),W(1,4),GC_364,MDL_MZ,MDL_WZ,W(1,11))
C     Amplitude(s) for diagram number 47
      CALL VVS4_0(W(1,11),W(1,9),W(1,7),GC_285,AMP(47))
      CALL FFV3_3(W(1,3),W(1,4),GC_353,MDL_MZ,MDL_WZ,W(1,11))
C     Amplitude(s) for diagram number 48
      CALL VVS4_0(W(1,11),W(1,9),W(1,7),GC_285,AMP(48))
      CALL FFV5_3(W(1,3),W(1,4),GC_392,MDL_MZ,MDL_WZ,W(1,11))
C     Amplitude(s) for diagram number 49
      CALL VVS4_0(W(1,11),W(1,9),W(1,7),GC_285,AMP(49))
      CALL FFV8_3(W(1,3),W(1,4),GC_719,MDL_MZ,MDL_WZ,W(1,11))
C     Amplitude(s) for diagram number 50
      CALL VVS4_0(W(1,11),W(1,9),W(1,7),GC_285,AMP(50))
      CALL FFV2_3(W(1,3),W(1,4),GC_718,MDL_MZ,MDL_WZ,W(1,11))
C     Amplitude(s) for diagram number 51
      CALL VVS4_0(W(1,11),W(1,9),W(1,7),GC_285,AMP(51))
      CALL FFV1_3(W(1,3),W(1,4),GC_395,MDL_MZ,MDL_WZ,W(1,11))
C     Amplitude(s) for diagram number 52
      CALL VVS4_0(W(1,11),W(1,9),W(1,7),GC_285,AMP(52))
      CALL FFV3_3(W(1,3),W(1,4),GC_384,MDL_MZ,MDL_WZ,W(1,11))
C     Amplitude(s) for diagram number 53
      CALL VVS4_0(W(1,11),W(1,9),W(1,7),GC_285,AMP(53))
      CALL FFV5_3(W(1,3),W(1,4),GC_398,MDL_MZ,MDL_WZ,W(1,11))
C     Amplitude(s) for diagram number 54
      CALL VVS4_0(W(1,11),W(1,9),W(1,7),GC_285,AMP(54))
      CALL FFV3_3(W(1,3),W(1,4),GC_387,MDL_MZ,MDL_WZ,W(1,11))
C     Amplitude(s) for diagram number 55
      CALL VVS4_0(W(1,11),W(1,9),W(1,7),GC_285,AMP(55))
      CALL FFS2_3(W(1,1),W(1,2),GC_521,MDL_MH,MDL_WH,W(1,11))
C     Amplitude(s) for diagram number 56
      CALL VVS4_0(W(1,8),W(1,9),W(1,11),GC_285,AMP(56))
      CALL FFS2_3(W(1,1),W(1,2),GC_522,MDL_MH,MDL_WH,W(1,11))
C     Amplitude(s) for diagram number 57
      CALL VVS4_0(W(1,8),W(1,9),W(1,11),GC_285,AMP(57))
      CALL FFS2_3(W(1,1),W(1,2),GC_523,MDL_MH,MDL_WH,W(1,11))
C     Amplitude(s) for diagram number 58
      CALL VVS4_0(W(1,8),W(1,9),W(1,11),GC_285,AMP(58))
      CALL FFS2_3(W(1,1),W(1,2),GC_524,MDL_MH,MDL_WH,W(1,11))
C     Amplitude(s) for diagram number 59
      CALL VVS4_0(W(1,8),W(1,9),W(1,11),GC_285,AMP(59))
      CALL FFS2_3(W(1,1),W(1,2),GC_525,MDL_MH,MDL_WH,W(1,11))
C     Amplitude(s) for diagram number 60
      CALL VVS4_0(W(1,8),W(1,9),W(1,11),GC_285,AMP(60))
      JAMP(1)=+AMP(1)+AMP(2)+AMP(3)+AMP(4)+AMP(5)+AMP(6)+AMP(7)+AMP(8)
     $ +AMP(9)+AMP(10)+AMP(11)+AMP(12)+AMP(13)+AMP(14)+AMP(15)+AMP(16)
     $ +AMP(17)+AMP(18)+AMP(19)+AMP(20)+AMP(21)+AMP(22)+AMP(23)+AMP(24)
     $ +AMP(25)+AMP(26)+AMP(27)+AMP(28)+AMP(29)+AMP(30)+AMP(31)+AMP(32)
     $ +AMP(33)+AMP(34)+AMP(35)+AMP(36)+AMP(37)+AMP(38)+AMP(39)+AMP(40)
     $ +AMP(41)+AMP(42)+AMP(43)+AMP(44)+AMP(45)+AMP(46)+AMP(47)+AMP(48)
     $ +AMP(49)+AMP(50)+AMP(51)+AMP(52)+AMP(53)+AMP(54)+AMP(55)+AMP(56)
     $ +AMP(57)+AMP(58)+AMP(59)+AMP(60)

      M27_MATRIX = 0.D0
      DO I = 1, NCOLOR
        ZTEMP = (0.D0,0.D0)
        DO J = 1, NCOLOR
          ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
        ENDDO
        M27_MATRIX = M27_MATRIX+ZTEMP*DCONJG(JAMP(I))
      ENDDO

      END

      SUBROUTINE M27_GET_VALUE(P, ALPHAS, NHEL ,ANS)
      IMPLICIT NONE
C     
C     CONSTANT
C     
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=6)
C     
C     ARGUMENTS 
C     
      REAL*8 P(0:3,NEXTERNAL),ANS
      INTEGER NHEL
      DOUBLE PRECISION ALPHAS
      REAL*8 PI
CF2PY INTENT(OUT) :: ANS
CF2PY INTENT(IN) :: NHEL
CF2PY INTENT(IN) :: P(0:3,NEXTERNAL)
CF2PY INTENT(IN) :: ALPHAS
C     ROUTINE FOR F2PY to read the benchmark point.    
C     the include file with the values of the parameters and masses 
      INCLUDE 'coupl.inc'

      PI = 3.141592653589793D0
      G = 2* DSQRT(ALPHAS*PI)
      CALL UPDATE_AS_PARAM()
      IF (NHEL.NE.0) THEN
        CALL M27_SMATRIXHEL(P, NHEL, ANS)
      ELSE
        CALL M27_SMATRIX(P, ANS)
      ENDIF
      RETURN
      END

      SUBROUTINE M27_INITIALISEMODEL(PATH)
C     ROUTINE FOR F2PY to read the benchmark point.    
      IMPLICIT NONE
      CHARACTER*512 PATH
CF2PY INTENT(IN) :: PATH
      CALL SETPARA(PATH)  !first call to setup the paramaters    
      RETURN
      END

      LOGICAL FUNCTION M27_IS_BORN_HEL_SELECTED(HELID)
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=6)
      INTEGER    NCOMB
      PARAMETER (NCOMB=64)
C     
C     ARGUMENTS
C     
      INTEGER HELID
C     
C     LOCALS
C     
      INTEGER I,J
      LOGICAL FOUNDIT
C     
C     GLOBALS
C     
      INTEGER HELC(NEXTERNAL,NCOMB)
      COMMON/M27_PROCESS_NHEL/HELC

      INTEGER POLARIZATIONS(0:NEXTERNAL,0:5)
      COMMON/M27_BORN_BEAM_POL/POLARIZATIONS
C     ----------
C     BEGIN CODE
C     ----------

      M27_IS_BORN_HEL_SELECTED = .TRUE.
      IF (POLARIZATIONS(0,0).EQ.-1) THEN
        RETURN
      ENDIF

      DO I=1,NEXTERNAL
        IF (POLARIZATIONS(I,0).EQ.-1) THEN
          CYCLE
        ENDIF
        FOUNDIT = .FALSE.
        DO J=1,POLARIZATIONS(I,0)
          IF (HELC(I,HELID).EQ.POLARIZATIONS(I,J)) THEN
            FOUNDIT = .TRUE.
            EXIT
          ENDIF
        ENDDO
        IF(.NOT.FOUNDIT) THEN
          M27_IS_BORN_HEL_SELECTED = .FALSE.
          RETURN
        ENDIF
      ENDDO

      RETURN
      END
