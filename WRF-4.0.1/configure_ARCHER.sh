#!/bin/bash
module load cray-netcdf/4.3.2
export WRFIO_NCD_LARGE_FILE_SUPPORT=1

export NETCDF=/opt/cray/netcdf/4.3.2/intel/140
./configure

