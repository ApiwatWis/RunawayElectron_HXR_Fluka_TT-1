*
*===magfld=============================================================*
*
      SUBROUTINE MAGFLD ( X, Y, Z, BTX, BTY, BTZ, B, NREG, IDISC )

      INCLUDE 'dblprc.inc'
      INCLUDE 'dimpar.inc'
      INCLUDE 'iounit.inc'
*
*----------------------------------------------------------------------*
*                                                                      *
*     Copyright (C) 2003-2019:  CERN & INFN                            *
*     All Rights Reserved.                                             *
*                                                                      *
*     Created  in     1988         by    Alberto Fasso`                *
*                                                                      *
*     Input variables:                                                 *
*            x,y,z = current position                                  *
*            nreg  = current region                                    *
*     Output variables:                                                *
*            btx,bty,btz = cosines of the magn. field vector           *
*            B = magnetic field intensity (Tesla)                      *
*            idisc = set to 1 if the particle has to be discarded      *
*                                                                      *
*----------------------------------------------------------------------*
*
      INCLUDE 'cmemfl.inc'
      INCLUDE 'csmcry.inc'

      REAL(KIND=8) theta, phi, Bt, Bp, s_m, s, rho, q, rhon
      REAL(KIND=8) Bx, By, Bz
*
*  +-------------------------------------------------------------------*
*  |  Earth geomagnetic field:
      IF ( LGMFLD ) THEN
         CALL GEOFLD ( X, Y, Z, BTX, BTY, BTZ, B, NREG, IDISC )
         RETURN
      END IF
*  |
*  +-------------------------------------------------------------------*
      IDISC = 0
*     uniform magnetic field:
*      BTX   = UMGFLD
*      BTY   = VMGFLD
*      BTZ   = WMGFLD
*      B     = BIFUNI
*  Apiwat defined the following lines, for B in tokamak

*     Convert (x, y, z) to (r, phi, theta)
      s = SQRT(X**2 + Y**2)
      phi = ATAN2(Y, X)
      theta = ATAN2(Z, s - 65.0D0)

*     Compute Bt (From Hall's paper), Bp (assumed q profile)
      s_m = s / 100.0
      Bt = -8.93946648*s_m**3 + 20.64395579*s_m**2 - 17.42630629*s_m 
      Bt = Bt + 6.21272779
      rho = (65.0D0*COS(phi)-X)**2 + (65.0D0*SIN(phi)-Y)**2 + Z**2
      rho = SQRT(rho)
      rhon = rho / 20.0D0    
      q = 1.0 + 3.0*rhon**2
      Bp = rho * Bt / 65.0D0 / q

*     Compute Bx, By, Bz
      Bx = Bt*SIN(phi) - Bp*SIN(theta)*COS(phi)
      By = -Bt*COS(phi) - Bp*SIN(theta)*SIN(phi)
      Bz = Bp*COS(theta)

*     Finally, compute B and its cosine values
      B     = SQRT(Bt**2 + Bp**2)
      BTX   = Bx / B
      BTY   = By / B
      BTZ   = Bz / B
      
      RETURN
*=== End of subroutine Magfld =========================================*
      END

