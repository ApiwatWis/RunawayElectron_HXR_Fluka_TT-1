*
*=== source ===========================================================*
*
      SUBROUTINE SOURCE ( NOMORE )        

      INCLUDE 'dblprc.inc'
      INCLUDE 'dimpar.inc'
      INCLUDE 'iounit.inc'
*
*----------------------------------------------------------------------*
*                                                                      *
*     Copyright (C) 2003-2019:  CERN & INFN                            *
*     All Rights Reserved.                                             *
*                                                                      *
*     New source for FLUKA9x-FLUKA20xy:                                *
*                                                                      *
*     Created on 07 January 1990   by    Alfredo Ferrari & Paola Sala  *
*                                                   Infn - Milan       *
*                                                                      *
*  This is just an example of a possible user written source routine.  *
*  note that the beam card still has some meaning - in the scoring the *
*  maximum momentum used in deciding the binning is taken from the     *
*  beam momentum.  Other beam card parameters are obsolete.            *
*                                                                      *
*       Output variables:                                              *
*                                                                      *
*              Nomore = if > 0 the run will be terminated              *
*                                                                      *
*----------------------------------------------------------------------*
*
      INCLUDE 'beamcm.inc'
      INCLUDE 'fheavy.inc'
      INCLUDE 'flkstk.inc'
      INCLUDE 'ioiocm.inc'
      INCLUDE 'ltclcm.inc'
      INCLUDE 'paprop.inc'
      INCLUDE 'sourcm.inc'
      INCLUDE 'sumcou.inc'
*     
      LOGICAL LFIRST, LISNUT
*
      SAVE LFIRST
      DATA LFIRST / .TRUE. /

*  Assume the column is along the Z-axis from ZMIN to ZMAX
      INTEGER :: ii
      REAL(KIND=8) RNDR, RNDT, RNDP, RDUMMY  
      REAL(KIND=8) Rminor, theta, phi
      REAL(KIND=8) BTX, BTY, BTZ, BB, NREG, IDISC
      REAL(KIND=8) Bx, By, Bz, Btheta, Bphi
      REAL(KIND=8) chi, xi, valpha, vbeta, vgamma
      REAL(KIND=8) valpha2, vbeta2, vgamma2 
      REAL(KIND=8) valpha3, vbeta3, vgamma3, v3
      REAL(KIND=8) Pmax, a1, p_r, u_acc
      LOGICAL ACCEPTED

*  Statement function:
      LISNUT (IJ) = INDEX ( PRNAME (IJ), 'NEUTRI' ) .GT. 0
*======================================================================*
*                                                                      *
*                 BASIC VERSION                                        *
*                                                                      *
*======================================================================*
      NOMORE = 0
*  +-------------------------------------------------------------------*
*  |  First call initializations:
      IF ( LFIRST ) THEN
*  |  *** The following 3 cards are mandatory ***
         TKESUM = ZERZER
         LFIRST = .FALSE.
         LUSSRC = .TRUE.
*  |  *** User initialization ***
         WRITE(*, *) 'User-defined source.f'
         WRITE(*, *) 'Author: Apiwat Wisitsorasak'
      END IF
*  |
*  +-------------------------------------------------------------------*
*  Push one source particle to the stack. Note that you could as well
*  push many but this way we reserve a maximum amount of space in the
*  stack for the secondaries to be generated
*  Npflka is the stack counter: of course any time source is called it
*  must be =0
      NPFLKA = NPFLKA + 1
*  Wt is the weight of the particle
      WTFLK  (NPFLKA) = ONEONE
      WEIPRI = WEIPRI + WTFLK (NPFLKA)
*  Particle type (1=proton.....). Ijbeam is the type set by the BEAM
*  card
*  +-------------------------------------------------------------------*
*  |  (Radioactive) isotope:
      IF ( IJBEAM .EQ. -2 .AND. LRDBEA ) THEN
         IARES  = IPROA
         IZRES  = IPROZ
         IISRES = IPROM
         CALL STISBM ( IARES, IZRES, IISRES )
         IJHION = IPROM  * 100000 + MOD ( IPROZ, 100 ) * 1000 + IPROA
         IJHION = IJHION * 100    + KXHEAV
         IONID  = IJHION
         CALL DCDION ( IONID )
         CALL SETION ( IONID )
         LFRPHN (NPFLKA) = .FALSE.
*  |
*  +-------------------------------------------------------------------*
*  |  Heavy ion:
      ELSE IF ( IJBEAM .EQ. -2 ) THEN
         IJHION = IPROM  * 100000 + MOD ( IPROZ, 100 ) * 1000 + IPROA
         IJHION = IJHION * 100    + KXHEAV
         IONID  = IJHION
         CALL DCDION ( IONID )
         CALL SETION ( IONID )
         ILOFLK (NPFLKA) = IJHION
*  |  Flag this is prompt radiation
         LRADDC (NPFLKA) = .FALSE.
*  |  Group number for "low" energy neutrons, set to 0 anyway
         IGROUP (NPFLKA) = 0
*  |  Parent radioactive isotope:
         IRDAZM (NPFLKA) = 0
*  |  Kinetic energy of the particle (GeV)
         TKEFLK (NPFLKA) = SQRT ( PBEAM**2 + AM (IONID)**2 )
     &                   - AM (IONID)
*  |  Particle momentum
         PMOFLK (NPFLKA) = PBEAM
*        PMOFLK (NPFLKA) = SQRT ( TKEFLK (NPFLKA) * ( TKEFLK (NPFLKA)
*    &                          + TWOTWO * AM (IONID) ) )
         LFRPHN (NPFLKA) = .FALSE.
*  |
*  +-------------------------------------------------------------------*
*  |  Normal hadron:
      ELSE
         IONID = IJBEAM
         ILOFLK (NPFLKA) = IJBEAM
*  |  Flag this is prompt radiation
         LRADDC (NPFLKA) = .FALSE.
*  |  Group number for "low" energy neutrons, set to 0 anyway
         IGROUP (NPFLKA) = 0
*  |  Parent radioactive isotope:
         IRDAZM (NPFLKA) = 0
*  |  Kinetic energy of the particle (GeV)
         TKEFLK (NPFLKA) = SQRT ( PBEAM**2 + AM (IONID)**2 )
     &                   - AM (IONID)
*  |  Particle momentum
         PMOFLK (NPFLKA) = PBEAM
*        PMOFLK (NPFLKA) = SQRT ( TKEFLK (NPFLKA) * ( TKEFLK (NPFLKA)
*    &                          + TWOTWO * AM (IONID) ) )
*  |  +----------------------------------------------------------------*
*  |  |  Check if it is a neutrino, if so force the interaction
*  |  |  (unless the relevant flag has been disabled)
         IF ( LISNUT (IJBEAM) .AND. LNUFIN ) THEN
            LFRPHN (NPFLKA) = .TRUE.
*  |  |
*  |  +----------------------------------------------------------------*
*  |  |  Not a neutrino
         ELSE
            LFRPHN (NPFLKA) = .FALSE.
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
      END IF
*  |
*  +-------------------------------------------------------------------*
*  From this point .....
*  Particle generation (1 for primaries)
      LOFLK  (NPFLKA) = 1
*  User dependent flag:
      LOUSE  (NPFLKA) = 0
*  No channeling:
      KCHFLK (NPFLKA) = 0
      ECRFLK (NPFLKA) = ZERZER
*  Extra infos:
      INFSTK (NPFLKA) = 0
      LNFSTK (NPFLKA) = 0
      ANFSTK (NPFLKA) = ZERZER
*  Parent variables:
      IPRSTK (NPFLKA) = 0
      EKPSTK (NPFLKA) = ZERZER
*  User dependent spare variables:
      DO 100 ISPR = 1, MKBMX1
         SPAREK (ISPR,NPFLKA) = ZERZER
 100  CONTINUE
*  User dependent spare flags:
      DO 200 ISPR = 1, MKBMX2
         ISPARK (ISPR,NPFLKA) = 0
 200  CONTINUE
*  Save the track number of the stack particle:
      ISPARK (MKBMX2,NPFLKA) = NPFLKA
      NPARMA = NPARMA + 1
      NUMPAR (NPFLKA) = NPARMA
      NEVENT (NPFLKA) = 0
      DFNEAR (NPFLKA) = +ZERZER
*  ... to this point: don't change anything
      AKNSHR (NPFLKA) = -TWOTWO
*  Particle age (s)
      AGESTK (NPFLKA) = +ZERZER
*  Compute directions and birth positions for a toroidal source
*  xxx Begin: [Apiwat] Modify the source
*  Particle geometry: toroidal-particle source

*  Particle coordinates
*     RNDR = FLRNDM (RDUMMY)  ! minor radius ! use the rejection sampling instead
      RNDT = FLRNDM (RDUMMY)  ! theta
      RNDP = FLRNDM (RDUMMY)  ! phi

* Maximum value of the radial distribution function
      Pmax = 1.0  ! Since the maximum of 1 - (r/a)^2 is 1 at r = 0
      a1 = 20.0   ! minor radius of the source (limiter)

* Start rejection sampling
      ACCEPTED = .FALSE.

      DO WHILE (.NOT. ACCEPTED)
* Sample RNDR uniformly between 0 and a
          RNDR = a1 * FLRNDM(RDUMMY)  ! Uniformly distributed between 0 and a

* Compute the target distribution for this value of r
          p_r = 1.0 - (RNDR / a1)**2

* Sample a random number u from uniform distribution [0, 1]
          u_acc = FLRNDM(RDUMMY)

* If u is less than or equal to the ratio of the target function and Pmax, accept the sample
          IF (u_acc <= p_r / Pmax) THEN
              ACCEPTED = .TRUE.
          END IF
      END DO

      Rminor = RNDR
      theta = 2.0 * 3.14159 * RNDT
      phi = 2.0 * 3.14159 * RNDP

      XFLK (NPFLKA) = (65.0 + Rminor * COS(theta)) * COS(phi)        ! X position on the circle
      YFLK (NPFLKA) = (65.0 + Rminor * COS(theta)) * SIN(phi)        ! Y position on the circle
      ZFLK (NPFLKA) = Rminor * SIN(theta)               ! Z position on the torus

*  Compute the magnetic field
      CALL MAGFLD(XFLK(NPFLKA), YFLK(NPFLKA), ZFLK(NPFLKA), 
     &            BTX, BTY, BTZ, BB, NREG, IDISC )
*  Compute the components of the magnetic field
      Bx = BB * BTX
      By = BB * BTY
      Bz = BB * BTZ
*  Compute the angles Btheta and Bphi of the B vector
      Btheta = ATAN2(SQRT(Bx**2 + By**2), Bz)
      Bphi = ATAN2(By, Bx)
      
*  Compute the velocity vector
*  - First, fixed the pitch angle (chi) to some small angle
      chi = WHASOU(1) / 180.0 * 3.14159
*  - Second, fixed the azimuthal angle (xi) to some random angle
      xi = 2.0 * 3.14159 * FLRNDM(RDUMMY)
*  - Compute components of the velocity vector
      valpha = SIN(chi)*COS(xi)
      vbeta = SIN(chi)*SIN(xi)
      vgamma = COS(chi)

      valpha2 = valpha*COS(Btheta) + vgamma*SIN(Btheta)
      vbeta2 = vbeta
      vgamma2 = -valpha*SIN(Btheta) + vgamma*COS(Btheta)

      valpha3 = valpha2*COS(Bphi) - vbeta2*SIN(Bphi)
      vbeta3 = valpha2*SIN(Bphi) + vbeta2*COS(Bphi)
      vgamma3 = vgamma2

      v3 = SQRT(valpha3**2 + vbeta3**2 + vgamma3**2)

*  Cosines (tx,ty,tz)
      TXFLK  (NPFLKA) = valpha3 / v3
      TYFLK  (NPFLKA) = vbeta3 / v3
      TZFLK  (NPFLKA) = vgamma3 / v3

*  Polarization cosines:
      TXPOL  (NPFLKA) = -TWOTWO
      TYPOL  (NPFLKA) = +ZERZER
      TZPOL  (NPFLKA) = +ZERZER
*  Energy, direct determine from BEAM card
      TKEFLK (NPFLKA) = ABS(PBEAM)
*  Particle momentum      
      PMOFLK (NPFLKA) = SQRT((TKEFLK(NPFLKA) + AM(IONID))**2  - AM(IONID)**2 )

*  Print out the source particle
      ii = NPFLKA
      WRITE(*, 111) ILOFLK(ii), XFLK(ii), YFLK(ii), ZFLK(ii), TXFLK(ii), TYFLK(ii), 
     &              TZFLK(ii), TKEFLK(ii), PMOFLK(ii), chi, xi, BTX, BTY, BTZ, BB,
     &              Btheta, Bphi, Rminor, theta, phi
 111  FORMAT(I4, 19(1X, E12.4))

*  xxx End: [Apiwat] Modify the source

*  Calculate the total kinetic energy of the primaries: don't change
*  +-------------------------------------------------------------------*
*  |  (Radioactive) isotope:
      IF ( IJBEAM .EQ. -2 .AND. LRDBEA ) THEN
*  |
*  +-------------------------------------------------------------------*
*  |  Heavy ion:
      ELSE IF ( ILOFLK (NPFLKA) .EQ. -2 .OR.
     &          ILOFLK (NPFLKA) .GT. 100000 ) THEN
         TKESUM = TKESUM + TKEFLK (NPFLKA) * WTFLK (NPFLKA)
*  |
*  +-------------------------------------------------------------------*
*  |  Standard particle:
      ELSE IF ( ILOFLK (NPFLKA) .NE. 0 ) THEN
         TKESUM = TKESUM + ( TKEFLK (NPFLKA) + AMDISC (ILOFLK(NPFLKA)) )
     &          * WTFLK (NPFLKA)
*  |
*  +-------------------------------------------------------------------*
*  |
      ELSE
         TKESUM = TKESUM + TKEFLK (NPFLKA) * WTFLK (NPFLKA)
      END IF
*  |
*  +-------------------------------------------------------------------*
      RADDLY (NPFLKA) = ZERZER
*  Here we ask for the region number of the hitting point.
*     NREG (NPFLKA) = ...
*  The following line makes the starting region search much more
*  robust if particles are starting very close to a boundary:
      CALL GEOCRS ( TXFLK (NPFLKA), TYFLK (NPFLKA), TZFLK (NPFLKA) )
      CALL GEOREG ( XFLK  (NPFLKA), YFLK  (NPFLKA), ZFLK  (NPFLKA),
     &              NRGFLK(NPFLKA), IDISC )
*  Do not change these cards:
      CALL GEOHSM ( NHSPNT (NPFLKA), 1, -11, MLATTC )
      NLATTC (NPFLKA) = MLATTC
      CMPATH (NPFLKA) = ZERZER
      CALL SOEVSV
      RETURN
*=== End of subroutine Source =========================================*
      END
