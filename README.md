# ROMS_supermod
Modified ROMS code (svn Revision 1041) to add the supermod operator.

CPP flag to activate supermod operator:
```sh
#define SUPERMOD
```

Information about the extra observation type is specified in s4dvar.in. For passive microwave SST:
```sh
NextraObs = 1
ExtraIndex = 21
ExtraName = CoarseSST
```

The size of the supermodding footprint is (1+2L)dx, where (1+2L) is the number of grid cells and dx is the horizontal resolution of the model. L is specified in the observation netCDF file _meta_ variable. The _type_ variable in the observation netCDF file is also modified to match the index of the extra observation type.  

------------

Our modifications to the file Nonlinear/set_data.F do not have to be included in order to implement the supermod operator algorithm. In Nonlinear/set_data.F a piece of code has been added to set the prior shortwave radiation flux (srflx)
```sh
#  ifdef SHORTWAVE
!
!  Set prior shortwave radiation flux.  For consistency, we need to
!  process the same values using in the computation of the net heat
!  flux since there is a time interpolation from snapshots.
!
!
      CALL set_2dfld_tile (ng, tile, iNLM, idSrad,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%srflxG,                           &
     &                     FORCES(ng)%srflx,                            &
     &                     update)
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               __FILE__)) RETURN

#  endif
```
since srflx is used while setting the vertical diffusive tracer fluxes (FC) for each time step (in Nonlinear/pre_step3d.F if SOLAR_SOURCE is defined).

[![DOI](https://zenodo.org/badge/537367363.svg)](https://zenodo.org/badge/latestdoi/537367363)
