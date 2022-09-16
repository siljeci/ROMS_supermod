#include "cppdefs.h"
      MODULE extract_obs_mod
#if (defined FOUR_DVAR || defined VERIFICATION) && defined OBSERVATIONS
!
!svn $Id: extract_obs.F 1037 2020-08-15 01:56:27Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2020 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine extracts model data at the requested observations      !
!  positions (Xobs,Yobs,Zobs).  The extraction is done via linear      !
!  interpolation. The (Xobs,Yobs) positions must be in fractional      !
!  grid coordinates.  Zobs can be in fractional  grid coordinates      !
!  (Zobs >= 0) or actual depths (Zobs < 0), if applicable.             !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     Imin       Global I-coordinate lower bound threshold for         !
!                  requested state field type.                         !
!     Imax       Global I-coordinate upper bound threshold for         !
!                  requested state field type.                         !
!     Jmin       Global J-coordinate lower bound threshold for         !
!                  requested state field type.                         !
!     Jmax       Global J-coordinate upper bound threshold for         !
!                  requested state field type.                         !
!     LBi        I-dimension Lower bound.                              !
!     UBi        I-dimension Upper bound.                              !
!     LBj        J-dimension Lower bound.                              !
!     UBj        J-dimension Upper bound.                              !
!     LBk        K-dimension Lower bound.                              !
!     UBk        K-dimension Upper bound.                              !
!     ifield     State field identification to process.                !
!     Mobs       Observation dimension in the calling program.         !
!     NobsSTR    Starting observation to process.                      !
!     NobsEND    Last observations to process.                         !
!     Xmin       Global minimum fractional I-coordinate to consider.   !
!     Xmax       Global maximum fractional I-coordinate to consider.   !
!     Ymin       Global minimum fractional J-coordinate to consider.   !
!     Ymax       Global maximum fractional J-coordinate to consider.   !
!     time       Current model time (seconds).                         !
!     dt         Model baroclinic time-step (seconds).                 !
!     ObsType    Observations type.                                    !
!     Tobs       Observations time (days).                             !
!     Xobs       Observations X-locations (grid coordinates).          !
!     Yobs       Observations Y-locations (grid coordinates).          !
!     Zobs       Observations Z-locations (grid coordinates or meters).!
!     A          Model array (2D or 3D) to process.                    !
!     Adepth     Depths (meter; negative).                             !
!     Amask      Land-sea masking.                                     !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     ObsVetting Observation screenning flag, 0: reject or 1: accept.  !
!     Aobs       Extracted model values at observation positions.      !
!     Zobs       Observations Z-locations (grid coordinates).          !
!                                                                      !
!  The interpolation weights matrix, Hmat(1:8), is as follows:         !
!                                                                      !
!                               8____________7                         !
!                               /.          /| (i2,j2,k2)              !
!                              / .         / |                         !
!                            5/___________/6 |                         !
!                             |  .        |  |                         !
!                             |  .        |  |         Grid Cell       !
!                             | 4.........|..|3                        !
!                             | .         |  /                         !
!                             |.          | /                          !
!                  (i1,j1,k1) |___________|/                           !
!                             1           2                            !
!                                                                      !
!  Notice that the indices i2 and j2 are reset when observations are   !
!  located exactly at the eastern and/or northern boundaries. This is  !
!  needed to avoid out-of-range array computations.                    !
!                                                                      !
!  All the observations are assumed to in fractional coordinates with  !
!  respect to RHO-points:                                              !
!                                                                      !
!                                                                      !
!  M      r..u..r..u..r..u..r..u..r..u..r..u..r..u..r..u..r..u..r      !
!         :                                                     :      !
!  Mm+.5  v  p++v++p++v++p++v++p++v++p++v++p++v++p++v++p++v++p  v      !
!         :  +     |     |     |     |     |     |     |     +  :      !
!  Mm     r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r      !
!         :  +     |     |     |     |     |     |     |     +  :      !
!  Mm-.5  v  p--v--p--v--p--v--p--v--p--v--p--v--p--v--p--v--p  v      !
!         :  +     |     |     |     |     |     |     |     +  :      !
!         r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r      !
!         :  +     |     |     |     |     |     |     |     +  :      !
!         v  p--v--p--v--p--v--p--v--p--v--p--v--p--v--p--v--p  v      !
!         :  +     |     |     |     |     |     |     |     +  :      !
!         r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r      !
!         :  +     |     |     |     |     |     |     |     +  :      !
!  2.5    v  p--v--p--v--p--v--p--v--p--v--p--v--p--v--p--v--p  v      !
!         :  +     |     |     |     |     |     |     |     +  :      !
!  2.0    r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r      !
!         :  +     |     |     |     |     |     |     |     +  :      !
!  1.5    v  p--v--p--v--p--v--p--v--p--v--p--v--p--v--p--v--p  v      !
!         :  +     |     |     |     |     |     |     |     +  :      !
!  1.0    r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r      !
!         :  +     |     |     |     |     |     |     |     +  :      !
!  0.5    v  p++v++p++v++p++v++p++v++p++v++p++v++p++v++p++v++p  v      !
!         :                                                     :      !
!  0.0    r..u..r..u..r..u..r..u..r..u..r..u..r..u..r..u..r..u..r      !
!                                                                      !
!           0.5   1.5   2.5                          Lm-.5 Lm+.5       !
!                                                                      !
!        0.0   1.0   2.0                                  Lm    L      !
!                                                                      !
!=======================================================================
!
      USE mod_kinds

      implicit none

      PUBLIC extract_obs2d
# ifdef SOLVE3D
      PUBLIC extract_obs3d
#  ifdef SUPERMOD
      PUBLIC extract_obs3d_supermod
#  endif
# endif

      CONTAINS
!
!***********************************************************************
      SUBROUTINE extract_obs2d (ng, Imin, Imax, Jmin, Jmax,             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ifield, Mobs, NobsSTR, NobsEND,         &
     &                          Xmin, Xmax, Ymin, Ymax,                 &
     &                          time, dt,                               &
     &                          ObsType, ObsVetting,                    &
     &                          Tobs, Xobs, Yobs,                       &
     &                          A,                                      &
# ifdef MASKING
     &                          Amask,                                  &
# endif
     &                          Aobs)
!***********************************************************************
!
      USE mod_ncparam,  ONLY : isUbar, isVbar
      USE mod_fourdvar, ONLY : ObsState2Type
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Imin, Imax, Jmin, Jmax
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: ifield, Mobs, NobsSTR, NobsEND
!
      real(r8), intent(in) :: Xmin, Xmax, Ymin, Ymax
      real(dp), intent(in) :: time, dt
!
# ifdef ASSUMED_SHAPE
      integer, intent(in) :: ObsType(:)

      real(dp), intent(in) :: Tobs(:)
      real(r8), intent(in) :: Xobs(:)
      real(r8), intent(in) :: Yobs(:)
      real(r8), intent(in) :: A(LBi:,LBj:)
#  ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:,LBj:)
#  endif

      real(r8), intent(inout) :: ObsVetting(:)
      real(r8), intent(inout) :: Aobs(:)
# else
      integer, intent(in) :: ObsType(Mobs)

      real(dp), intent(in) :: Tobs(Mobs)
      real(r8), intent(in) :: Xobs(Mobs)
      real(r8), intent(in) :: Yobs(Mobs)
      real(r8), intent(in) :: A(LBi:UBi,LBj:UBj)
#  ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(inout) :: ObsVetting(Mobs)
      real(r8), intent(inout) :: Aobs(Mobs)
# endif
!
!  Local variable declarations.
!
      integer :: ic, iobs, i1, i2, j1, j2

      real(dp) :: TimeLB, TimeUB

      real(r8) :: p1, p2, q1, q2, wsum

      real(r8), dimension(8) :: Hmat
!
!-----------------------------------------------------------------------
!  Interpolate from requested 2D state field when appropriate.
!-----------------------------------------------------------------------
!
      TimeLB=(time-0.5_dp*dt)/86400.0_dp
      TimeUB=(time+0.5_dp*dt)/86400.0_dp
!
      DO iobs=NobsSTR,NobsEND
        IF ((ObsType(iobs).eq.ifield).and.                              &
     &      ((TimeLB.le.Tobs(iobs)).and.(Tobs(iobs).lt.TimeUB)).and.    &
     &      ((Xmin.le.Xobs(iobs)).and.(Xobs(iobs).lt.Xmax)).and.        &
     &      ((Ymin.le.Yobs(iobs)).and.(Yobs(iobs).lt.Ymax))) THEN
          IF (ObsType(iobs).eq.ObsState2Type(isUbar)) THEN
            i1=INT(Xobs(iobs)+0.5_r8)        ! 2D U-grid type variable
            j1=INT(Yobs(iobs))
          ELSE IF (ObsType(iobs).eq.ObsState2Type(isVbar)) THEN
            i1=INT(Xobs(iobs))               ! 2D V-grid type variable
            j1=INT(Yobs(iobs)+0.5_r8)
          ELSE
            i1=INT(Xobs(iobs))               ! 2D RHO-grid type variable
            j1=INT(Yobs(iobs))
          END IF
          i2=i1+1
          j2=j1+1
          IF (i2.gt.Imax) THEN
            i2=i1                 ! Observation at the eastern boundary
          END IF
          IF (j2.gt.Jmax) THEN
            j2=j1                 ! Observation at the northern boundary
          END IF
          p2=REAL(i2-i1,r8)*(Xobs(iobs)-REAL(i1,r8))
          q2=REAL(j2-j1,r8)*(Yobs(iobs)-REAL(j1,r8))
          p1=1.0_r8-p2
          q1=1.0_r8-q2
          Hmat(1)=p1*q1
          Hmat(2)=p2*q1
          Hmat(3)=p2*q2
          Hmat(4)=p1*q2
# ifdef MASKING
          Hmat(1)=Hmat(1)*Amask(i1,j1)
          Hmat(2)=Hmat(2)*Amask(i2,j1)
          Hmat(3)=Hmat(3)*Amask(i2,j2)
          Hmat(4)=Hmat(4)*Amask(i1,j2)
          wsum=0.0_r8
          DO ic=1,4
            wsum=wsum+Hmat(ic)
          END DO
          IF (wsum.gt.0.0_r8) THEN
            wsum=1.0_r8/wsum
            DO ic=1,4
              Hmat(ic)=Hmat(ic)*wsum
            END DO
          END IF
# endif
          Aobs(iobs)=Hmat(1)*A(i1,j1)+                                  &
     &               Hmat(2)*A(i2,j1)+                                  &
     &               Hmat(3)*A(i2,j2)+                                  &
     &               Hmat(4)*A(i1,j2)
# ifdef MASKING
          IF (wsum.gt.0.0_r8) ObsVetting(iobs)=1.0_r8
# else
          ObsVetting(iobs)=1.0_r8
# endif
        END IF
      END DO

      RETURN
      END SUBROUTINE extract_obs2d

# ifdef SOLVE3D
!
!***********************************************************************
      SUBROUTINE extract_obs3d (ng, Imin, Imax, Jmin, Jmax,             &
     &                          LBi, UBi, LBj, UBj, LBk, UBk,           &
     &                          ifield, Mobs, NobsSTR, NobsEND,         &
     &                          Xmin, Xmax, Ymin, Ymax,                 &
     &                          time, dt,                               &
     &                          ObsType, ObsVetting,                    &
     &                          Tobs, Xobs, Yobs, Zobs,                 &
     &                          A, Adepth,                              &
#  ifdef MASKING
     &                          Amask,                                  &
#  endif
     &                          Aobs)
!***********************************************************************
!
      USE mod_param
!
      USE mod_ncparam,  ONLY : isUvel, isVvel
      USE mod_fourdvar, ONLY : ObsState2Type
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Imin, Imax, Jmin, Jmax
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
      integer, intent(in) :: ifield, Mobs, NobsSTR, NobsEND
!
      real(r8), intent(in) :: Xmin, Xmax, Ymin, Ymax
      real(dp), intent(in) :: time, dt
!
#  ifdef ASSUMED_SHAPE
      integer, intent(in) :: ObsType(:)

      real(dp), intent(in) :: Tobs(:)
      real(r8), intent(in) :: Xobs(:)
      real(r8), intent(in) :: Yobs(:)
      real(r8), intent(in) :: A(LBi:,LBj:,LBk:)
      real(r8), intent(in) :: Adepth(LBi:,LBj:,LBk:)
#   ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:,LBj:)
#   endif
      real(r8), intent(inout) :: ObsVetting(:)
      real(r8), intent(inout) :: Zobs(:)
      real(r8), intent(inout) :: Aobs(:)
#  else
      integer, intent(in) :: ObsType(Mobs)

      real(dp), intent(in) :: Tobs(Mobs)
      real(r8), intent(in) :: Xobs(Mobs)
      real(r8), intent(in) :: Yobs(Mobs)
      real(r8), intent(in) :: A(LBi:UBi,LBj:UBj,LBk:UBk)
      real(r8), intent(in) :: Adepth(LBi:UBi,LBj:UBj,LBk:UBk)
#   ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(inout) :: ObsVetting(Mobs)
      real(r8), intent(inout) :: Zobs(Mobs)
      real(r8), intent(inout) :: Aobs(Mobs)
#  endif
!
!  Local variable declarations.
!
      integer :: i, ic, iobs, i1, i2, j1, j2, k, k1, k2

      real(dp) :: TimeLB, TimeUB

      real(r8) ::  Zbot, Ztop, dz, p1, p2, q1, q2, r1, r2
      real(r8) :: w11, w12, w21, w22, wsum

      real(r8), dimension(8) :: Hmat
!
!-----------------------------------------------------------------------
!  Interpolate from requested 3D state field.
!-----------------------------------------------------------------------
!
      TimeLB=(time-0.5_dp*dt)/86400.0_dp
      TimeUB=(time+0.5_dp*dt)/86400.0_dp
!
      DO iobs=NobsSTR,NobsEND
        IF ((ObsType(iobs).eq.ifield).and.                              &
     &      ((TimeLB.le.Tobs(iobs)).and.(Tobs(iobs).lt.TimeUB)).and.    &
     &      ((Xmin.le.Xobs(iobs)).and.(Xobs(iobs).lt.Xmax)).and.        &
     &      ((Ymin.le.Yobs(iobs)).and.(Yobs(iobs).lt.Ymax))) THEN
          IF (ObsType(iobs).eq.ObsState2Type(isUvel)) THEN
            i1=INT(Xobs(iobs)+0.5_r8)        ! 3D U-grid type variable
            j1=INT(Yobs(iobs))
          ELSE IF (ObsType(iobs).eq.ObsState2Type(isVvel)) THEN
            i1=INT(Xobs(iobs))               ! 3D V-grid type variable
            j1=INT(Yobs(iobs)+0.5_r8)
          ELSE
            i1=INT(Xobs(iobs))               ! 3D RHO-grid type variable
            j1=INT(Yobs(iobs))
          END IF
          i2=i1+1
          j2=j1+1
          IF (i2.gt.Imax) THEN
            i2=i1                 ! Observation at the eastern boundary
          END IF
          IF (j2.gt.Jmax) THEN
            j2=j1                 ! Observation at the northern boundary
          END IF
          p2=REAL(i2-i1,r8)*(Xobs(iobs)-REAL(i1,r8))
          q2=REAL(j2-j1,r8)*(Yobs(iobs)-REAL(j1,r8))
          p1=1.0_r8-p2
          q1=1.0_r8-q2
          w11=p1*q1
          w21=p2*q1
          w22=p2*q2
          w12=p1*q2
          IF (Zobs(iobs).gt.0.0_r8) THEN
            k1=MAX(1,INT(Zobs(iobs)))          ! Positions in fractional
            k2=MIN(INT(Zobs(iobs))+1,N(ng))    ! levels
            r2=REAL(k2-k1,r8)*(Zobs(iobs)-REAL(k1,r8))
            r1=1.0_r8-r2
          ELSE
            Ztop=Adepth(i1,j1,N(ng))
            Zbot=Adepth(i1,j1,1    )
            IF (Zobs(iobs).ge.Ztop) THEN
              r1=0.0_r8                        ! If shallower, ignore.
              r2=0.0_r8
              ObsVetting(iobs)=0.0_r8
            ELSE IF (Zbot.ge.Zobs(iobs)) THEN
              r1=0.0_r8                        ! If deeper, ignore.
              r2=0.0_r8
              ObsVetting(iobs)=0.0_r8
            ELSE
              DO k=N(ng),2,-1                  ! Otherwise, interpolate
                Ztop=Adepth(i1,j1,k  )         ! to fractional level
                Zbot=Adepth(i1,j1,k-1)
                IF ((Ztop.gt.Zobs(iobs)).and.(Zobs(iobs).ge.Zbot)) THEN
                  k1=k-1
                  k2=k
                END IF
              END DO
              dz=Adepth(i1,j1,k2)-Adepth(i1,j1,k1)
              r2=(Zobs(iobs)-Adepth(i1,j1,k1))/dz
              r1=1.0_r8-r2
              Zobs(iobs)=REAL(k1,r8)+r2        ! overwrite
            END IF
          END IF
          IF ((r1+r2).gt.0.0_r8) THEN
            Hmat(1)=w11*r1
            Hmat(2)=w21*r1
            Hmat(3)=w22*r1
            Hmat(4)=w12*r1
            Hmat(5)=w11*r2
            Hmat(6)=w21*r2
            Hmat(7)=w22*r2
            Hmat(8)=w12*r2
#  ifdef MASKING
            Hmat(1)=Hmat(1)*Amask(i1,j1)
            Hmat(2)=Hmat(2)*Amask(i2,j1)
            Hmat(3)=Hmat(3)*Amask(i2,j2)
            Hmat(4)=Hmat(4)*Amask(i1,j2)
            Hmat(5)=Hmat(5)*Amask(i1,j1)
            Hmat(6)=Hmat(6)*Amask(i2,j1)
            Hmat(7)=Hmat(7)*Amask(i2,j2)
            Hmat(8)=Hmat(8)*Amask(i1,j2)
            wsum=0.0_r8
            DO ic=1,8
              wsum=wsum+Hmat(ic)
            END DO
            IF (wsum.gt.0.0_r8) THEN
              wsum=1.0_r8/wsum
              DO ic=1,8
                Hmat(ic)=Hmat(ic)*wsum
              END DO
            END IF
#  endif
            Aobs(iobs)=Hmat(1)*A(i1,j1,k1)+                             &
     &                 Hmat(2)*A(i2,j1,k1)+                             &
     &                 Hmat(3)*A(i2,j2,k1)+                             &
     &                 Hmat(4)*A(i1,j2,k1)+                             &
     &                 Hmat(5)*A(i1,j1,k2)+                             &
     &                 Hmat(6)*A(i2,j1,k2)+                             &
     &                 Hmat(7)*A(i2,j2,k2)+                             &
     &                 Hmat(8)*A(i1,j2,k2)
#  ifdef MASKING
            IF (wsum.gt.0.0_r8) ObsVetting(iobs)=1.0_r8
#  else
            ObsVetting(iobs)=1.0_r8
#  endif
#  ifndef ALLOW_BOTTOM_OBS
!
!  Reject observations that lie in the lower bottom grid cell (k=1) to
!  avoid clustering due shallowing of bathymetry during smoothing and
!  coarse level half-thickness (-h < Zobs < Adepth(:,:,1)) in deep
!  water.
!
            IF ((Zobs(iobs).gt.0.0_r8).and.(Zobs(iobs).le.1.0_r8)) THEN
              ObsVetting(iobs)=0.0_r8
            END IF
#  endif
          END IF
        END IF
      END DO
!
      RETURN

      END SUBROUTINE extract_obs3d
# ifdef SUPERMOD

!***********************************************************************
      SUBROUTINE extract_obs3d_supermod (ng, tile, model,               &
     &                          Imin, Imax, Jmin, Jmax,                 &
     &                          LBi, UBi, LBj, UBj, LBk, UBk,           &
     &                          ifield, Mobs, NobsSTR, NobsEND,         &
     &                          Xmin, Xmax, Ymin, Ymax,                 &
     &                          time, dt,                               &
     &                          ObsType, ObsVetting,                    &
     &                          Tobs, Xobs, Yobs, Zobs,                 &
     &                          A, Adepth,                              &
#  ifdef MASKING
     &                          Amask,                                  &
#  endif
     &                          Aobs)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars,      ONLY : EWperiodic, NSperiodic
      USE mod_ncparam,      ONLY : isUvel, isVvel
      USE mod_fourdvar,     ONLY : ObsState2Type,ObsMeta
      USE mp_exchange_mod,  ONLY : mp_exchange2d_supermod
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Imin, Imax, Jmin, Jmax, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
      integer, intent(in) :: ifield, Mobs, NobsSTR, NobsEND
!
      real(r8), intent(in) :: Xmin, Xmax, Ymin, Ymax
      real(dp), intent(in) :: time, dt
!
#  ifdef ASSUMED_SHAPE
      integer, intent(in) :: ObsType(:)

      real(dp), intent(in) :: Tobs(:)
      real(r8), intent(in) :: Xobs(:)
      real(r8), intent(in) :: Yobs(:)
      real(r8), intent(in) :: A(LBi:,LBj:,LBk:)
      real(r8), intent(in) :: Adepth(LBi:,LBj:,LBk:)
#   ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:,LBj:)
#   endif
      real(r8), intent(inout) :: ObsVetting(:)
      real(r8), intent(inout) :: Zobs(:)
      real(r8), intent(inout) :: Aobs(:)
#  else
      integer, intent(in) :: ObsType(Mobs)

      real(dp), intent(in) :: Tobs(Mobs)
      real(r8), intent(in) :: Xobs(Mobs)
      real(r8), intent(in) :: Yobs(Mobs)
      real(r8), intent(in) :: A(LBi:UBi,LBj:UBj,LBk:UBk)
      real(r8), intent(in) :: Adepth(LBi:UBi,LBj:UBj,LBk:UBk)
#   ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(inout) :: ObsVetting(Mobs)
      real(r8), intent(inout) :: Zobs(Mobs)
      real(r8), intent(inout) :: Aobs(Mobs)
#  endif
!
!  Local variable declarations.
!
      integer :: iobs, i, j, ic, jc, Itile, Jtile
      integer :: i1, j1, p1, p2, q1, q2
      integer :: LBi_ext, UBi_ext, LBj_ext, UBj_ext
      integer :: footprint_length, MaxFootprint, Nghost_supermod
      real(r8) :: wsum
      real(dp) :: TimeLB, TimeUB
      real(r8), allocatable :: A_extended(:,:) ! It will only cover the surface
      real(r8) :: Hmat(1:2,1:2)
      integer :: SWmat(1:2,1:2)
#   ifdef MASKING
      real(r8), allocatable :: Amask_extended(:,:)
#   endif
!     
!-----------------------------------------------------------------------
!  All the model values from grid points located in a square of size 
! 2* footprint_length around the observation are averaged.
! 
! NE_j   _ _ _ _ _ _ _ _ _ 
! 	|                 |		
! 	|                 |
!	|                 |
! 	|                 |
! 	|        +<------>|   L = footprint_length
! 	|             L   |
! 	|                 |
! 	|                 |
! SW_j	|_ _ _ _ _ _ _ _ _|    
!       SW_i    Xobs     NE_i
! 
! The length of the halo of each tile  for parallel computations in ROMS 
! is NghostPoints and is set in inp_par.F By default it is equal to 2,
!  and can be equal to 3 if periodic boundary conditions are used. 
! This might not be enough for us if the footprint length will probably
! exceed 2. Consequently, we need to temporarily work with a larger halo,
! using the arrays A_extended and Amask_extended.
! A_extended is initialized as:
!    _ _ _ _ _ _ _ _
!   |       0       | 
!   |     _ _ _     | 
!   |    |     |    | 
!   | 0  |  A  |  0 | 
!   |    |_ _ _|    | 
!   |               | 
!   |_ _ _ _0_ _ _ _| 
!    <-->  
!   L - NghostPoints  
! 
!   Then, mp_exchange is called to fill the 'border' outside of A.
!   Note that a special care should be taken when the tile is on the 
!   edge of the grid so that the array does not 'extend' outside of the
!   grid. 
!
! If observations are located between grid points, we need to combine to 
! interpolation and horizontal averaging. This method can be represented 
! more clearly in 1D. If our footprint is a segment with a length of three 
! points, then there are three points in the footprint if the observation
! is located exactly on a grid point, and two otherwise:
!
! grid:    +    +    +    +    +           +    +    +    +    +
! obs :              o                                 o   
! footprint:    |_________|                       |_________| 
!
! If we only average model values on two grid points in this case, then
! two observations would generate increments with significantly different
! amplitudes if one is located on a grid points and one between two. 
! Moreover, the 'effective' footrpint would not be centered exactly on the
! observaiton. To avoid that, an equivalent to hte third point can be 
! created by interpolating the two closest points:
! 
!                   x1       x2       x3       x4       x5
!                   +        +        +  o     +        +
!                               |_________________| 
!
!  x_average= ( 0.75*x2 + x3 + x4 +0.25*x5 ) /3
! 
! The same result can be obtained by computing the model averages for 
! footprints centered respectively on x3 and x4, and interpolate between
! the two:  
! 
!  x_average=  0.25 * (x2 + x3 + x4)/3  + +0.75*(x3 + x4 + x5)/3
!
!This second method is preferred in our case as it is the simplest
! of the two in two dimensions. The suffixes 1 and 2 play the same
! roles as in the interpolation in extract_obs3d:
!        
!                                1 point difference
!                                    <-> 
!     NE_j2  ___________________________
!           |   |                   |   |
!     NE_j1 |___|___________________|___|
!           |   |                   |   |
!           |   |                   |   |
!       j2  |   |                   |   |
!           |   |       +   +       |   |
!           |   |         o         |   |
!       j1  |   |       +   +       |   |
!           |   |                   |   |
!           |   |                   |   |
!     SW_j2 |___|___________________|___|
!           |   |                   |   |
!     SW_j1 |___|___________________|___|
! 
!        SW_i1 SW_i2   i1   i2    NE_i1 NE_i2
!           <---------->     <---------->
!        footprint_length   footprint_length  
! 
!---------------------------------------------------------------------
! 

!-----------------------------------------------------------------------
! First we need to determine the dimensions of the tiles 'extended' by
! the additional grid points. Get_bounds is the routine called in inp_par
!to definethe 'original' values of LBi,UBi,LBj,Ubj. It is here called 
! for a larger number of ghost points to get updated bounds 
!-----------------------------------------------------------------------


      ! The largest metadata value from all processed observations is computed.	
      MaxFootprint = 0.0_r8
      DO iobs=NobsSTR,NobsEND
        IF ((ObsType(iobs).eq.ifield).and.                              &
     &      ( MaxFootprint.lt.ObsMeta(iobs) )) THEN
           MaxFootprint = ObsMeta(iobs)
        END IF
      END DO      
      
      ! The '+1' is required as the limit for fractional observation locations
      ! goes beyond the 'integer locations' bounds i.e. Xmax=UBi+0.5
      Nghost_supermod = MAX ( MaxFootprint, Nghostpoints )+1

      ! We can use this value to find the 'extended' bounds of the tiles.
          CALL get_bounds (ng, tile, 0, Nghost_supermod, Itile, Jtile,      &
     &                     LBi_ext, UBi_ext, LBj_ext, UBj_ext)


      ! The extended arrays A_extended and Amask_extended are created
      allocate ( A_extended (LBi_ext:UBi_ext,LBj_ext:UBj_ext) )
      A_extended (:,:)=0.0_r8 !Could be any value, it will be overwritten
      A_extended (LBi:UBi,LBj:UBj) = A (LBi:UBi,LBj:UBj,N(ng))


      allocate ( Amask_extended (LBi_ext:UBi_ext,LBj_ext:UBj_ext))
      Amask_extended (:,:)=0.0_r8 !Could be any value, it will be overwritten
      Amask_extended (LBi:UBi,LBj:UBj) = Amask (LBi:UBi,LBj:UBj)


!     The values of the ghost points are fetched from the neighboring tiles
!     for both A_extended and Amask_extended.
      CALL mp_exchange2d_supermod (ng, tile, model, 2,                      &
      &                   LBi_ext, UBi_ext, LBj_ext, UBj_ext,               &
      &                          Nghost_supermod, EWperiodic(ng),           &
      &                          NSperiodic(ng), A_extended, Amask_extended)

!---------------------------------------------------------------------------

      TimeLB=(time-0.5_dp*dt)/86400.0_dp
      TimeUB=(time+0.5_dp*dt)/86400.0_dp

      DO iobs=NobsSTR,NobsEND
        IF ((ObsType(iobs).eq.ifield).and.                              &
     &      ((TimeLB.le.Tobs(iobs)).and.(Tobs(iobs).lt.TimeUB)).and.    &
     &      ((Xmin.le.Xobs(iobs)).and.(Xobs(iobs).lt.Xmax)).and.        &
     &      ((Ymin.le.Yobs(iobs)).and.(Yobs(iobs).lt.Ymax))) THEN
  
!---------------------------------------------------------------------------
! Corner and centers locations
!---------------------------------------------------------------------------
          footprint_length = ObsMeta(iobs) ! footprint length in grid points number

          i1 = INT(Xobs(iobs))
          j1 = INT(Yobs(iobs))

          SWmat (1,1) = i1 - footprint_length 
          SWmat (1,2) = j1 - footprint_length
          SWmat (2,1) = SWmat (1,1) + 1
          SWmat (2,2) = SWmat (1,2) + 1

!--------------------------------------------------------------------------
! Interpolation weights
!--------------------------------------------------------------------------
          p2 = Xobs(iobs)-REAL(i1,r8)
          q2 = Yobs(iobs)-REAL(j1,r8)
          p1 = 1.0_r8 - p2
          q1 = 1.0_r8 - q2
          Hmat(1,1) = p1 * q1
          Hmat(2,1) = p2 * q1
          Hmat(2,2) = p2 * q2
          Hmat(1,2) = p1 * q2
!-------------------------------------------------------------------------

! We need to make sure that the footprints do not reach outside the boundaries
        IF ( (Imin.le.SWmat(1,1)).and.(Jmin.le.SWmat(1,2)).and.          &
     &       (SWmat(2,1)+2*footprint_length.le.Imax).and.                &
     &       (SWmat(2,2)+2*footprint_length.le.Jmax) ) THEN

          Aobs(iobs) = 0

#   ifdef MASKING
          wsum = 0
          DO ic = 1, 2
             DO jc = 1, 2
                ! For each of the footprints, we loop on each of their 
                ! grid points.
                DO i = SWmat(ic,1), SWmat(ic,1) + 2 * footprint_length
                   DO j = SWmat(jc,2), SWmat(jc,2) + 2 * footprint_length
                      wsum = wsum + Amask_extended(i,j) * Hmat(ic,jc)
                      Aobs(iobs)=Aobs(iobs)+A_extended(i,j)*Hmat(ic,jc) &
                                        & * Amask_extended(i,j)
                   END DO
                END DO
             END DO
          END DO
          ! The sum of the weights should add up to the number of points in
          ! one footprint. Otherwise at least one point is masked, and the
          ! observation should be rejected
          IF (wsum.eq.( (1+2*footprint_length) ** 2 )) THEN
                Aobs(iobs) = Aobs(iobs) / wsum
                ObsVetting(iobs)=1.0_r8
          ELSE 
                ObsVetting(iobs)=0.0_r8
          END IF
#   else
          DO ic = 1, 2
             DO jc = 1, 2
                ! For each of the footprints, we loop on each of their 
                ! grid points.
                DO i = SWmat(ic,1), SWmat(ic,1) + 2 * footprint_length
                   DO j = SWmat(jc,2), SWmat(jc,2) + 2 * footprint_length
                      Aobs(iobs) = Aobs(iobs) + A_extended(i,j)* Hmat(ic,jc) 
                   END DO
                END DO
             END DO
          END DO
          wsum = (1+2*footprint_length) ** 2
          Aobs(iobs) = Aobs(iobs) / points_num
          ObsVetting(iobs)=1.0_r8
#   endif

        ELSE !i.e. if the observation footprint lie outside the grid
           ObsVetting(iobs) = 0.0_r8
        END IF
        END IF
      END DO
      RETURN
      deallocate ( A_extended )
      deallocate(Amask_extended)
      END SUBROUTINE extract_obs3d_supermod
# endif
# endif
#endif
      END MODULE extract_obs_mod
