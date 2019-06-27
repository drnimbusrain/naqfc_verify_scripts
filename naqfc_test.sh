#!/bin/bash -x
source /data/aqf2/barryb/anaconda2/bin/activate website

export OMP_NUM_THREADS=2
# get the current date (loop through current, or as far back as needed to fill)
date
yyyymmdd=$(date -d "-2 days" +%Y-%m-%d)
yyyymmdd48hr=$(date -d "-0 days" +%Y-%m-%d)
yyyymmddform=$(date -d "-2 days" +%Y%m%d)

# go to data directory
data_dir=/data/aqf3/patrickc/5xpm
mkdir $data_dir/$yyyymmddform
cd $data_dir
ln -sf /data/aqf2/barryb/5xpm/$yyyymmddform/* $yyyymmddform/.
 
# netcdf files
files_chem=$yyyymmddform/"aqm.t12z.aconc-sfc.ncf"
files_so2=$yyyymmddform/"aqm.t12z.aconc.ncf"
files_nox=$yyyymmddform/"aqm.t12z.nox.ncf"
files_met=$yyyymmddform/"aqm.t12z.metcro2d.ncf"
files_precip=$yyyymmddform/"aqm.t12z.precip.ncf"

#naqfc verify scripts 
pair=/data/aqf/patrickc/naqfc_verify_scripts/01.verify_pair.py
stats=/data/aqf/patrickc/naqfc_verify_scripts/02.verify_stats.py
taylor=/data/aqf/patrickc/naqfc_verify_scripts/03.verify_taylor_plots.py
bias=/data/aqf/patrickc/naqfc_verify_scripts/04.verify_spatial_bias.py
spatial=/data/aqf/patrickc/naqfc_verify_scripts/05.verify_spatial_plots.py
box=/data/aqf/patrickc/naqfc_verify_scripts/06.verify_box_plots.py

# pair the data
${pair} -f ${files_chem} -s 'OZONE' 'PM2.5' 

# spatial overlay plots
${spatial} -f ${files_chem} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair_chem.hdf -s 'OZONE' 'PM2.5'  -n ${yyyymmddform}

# spatial bias plots
${bias} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair_chem.hdf -s 'OZONE' 'PM2.5' -n ${yyyymmddform}

# create taylor
${taylor} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair_chem.hdf    -s 'OZONE' 'PM2.5' -n ${yyyymmddform}

