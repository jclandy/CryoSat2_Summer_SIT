Landy et al. [2022] A year-round satellite sea ice thickness record from CryoSat-2
_______________________________________________________________________________________

MATLAB codes are provided to convert CryoSat-2 radar freeboards for the Arctic region to estimates of sea ice thickness.

Start with 'RFB_to_SIT.m' which converts radar freeboards from the 'Developers Product' netCDF files to sea ice thickness estimates provided in the 'Sea Ice Thickness' netCDF files.

Codes are then provided to exactly reproduce Figures 1-3 in the main paper.
_______________________________________________________________________________________

Auxiliary functions (inpaint_nans, nanmedfilt2, ncpolarm, polarstereo_fwd, polarstereo_inv) are used within the four other scripts.

The 'PIOMAS Grid' directory contains binary files for the native PIOMAS grid and is used in 'Figure2.m'.

The 'MASIE_region_mask.mat' file contains a mask for the NSIDC Multisensor Analyzed Sea Ice Extent - Northern Hemisphere (MASIE-NH) regions and is used in 'Figure3.m'.