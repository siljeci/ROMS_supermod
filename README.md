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

[![DOI](https://zenodo.org/badge/537367363.svg)](https://zenodo.org/badge/latestdoi/537367363)
