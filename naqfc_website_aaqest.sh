#!/bin/bash -x


export OMP_NUM_THREADS=2
# get the current date
date
yyyymmdd=$(date -d "-112 days" +%Y-%m-%d)
yyyymmddform=$(date -d "-112 days" +%Y%m%d)
yyyymmdd72hr=$(date -d "-109 days" +%Y-%m-%d)

# go to data directory
data_dir=/data/aqf/patrickc/CMAQ_DATA/cmaq5.0.2-72hr
#data_dir=/gpfs/hps2/ptmp/Patrick.C.Campbell/fengsha/gfs.$yyyymmdd/00
cd $data_dir

# netcdf files
files="aqm.${yyyymmddform}.t12z.aconc-pm25.ncf"

#naqfc verify scripts 
pair=/data/aqf/patrickc/naqfc_verify_scripts/01.verify_pair.py
stats=/data/aqf/patrickc/naqfc_verify_scripts/02.verify_stats.py
taylor=/data/aqf/patrickc/naqfc_verify_scripts/03.verify_taylor_plots.py
bias=/data/aqf/patrickc/naqfc_verify_scripts/04.verify_spatial_bias.py
spatial=/data/aqf/patrickc/naqfc_verify_scripts/05.verify_spatial_plots.py
box=/data/aqf/patrickc/naqfc_verify_scripts/06.verify_box_plots.py

# pair the data
${pair} -f ${files} 

# spatial overlay plots
${spatial} -f ${files} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd72hr}-12_pair.hdf -s {'OZONE','PM2.5'} -n ${yyyymmddform} 
for i in 'R1' 'R2' 'R3' 'R4' 'R5' 'R6' 'R7' 'R8' 'R9' 'R10'; do
  echo ${i}
  ${spatial} -f ${files} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd72hr}-12_pair.hdf -s {'OZONE','PM2.5'} -b True -e ${i} -n ${yyyymmddform}
done

# spatial bias plots
${bias} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd72hr}-12_pair.hdf -s {'OZONE','PM2.5'} -n ${yyyymmddform}
for i in 'R1' 'R2' 'R3' 'R4' 'R5' 'R6' 'R7' 'R8' 'R9' 'R10'; do
  echo ${i}
  ${bias} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd72hr}-12_pair.hdf -s {'OZONE','PM2.5'}  -b True -e ${i} -n ${yyyymmddform}
done 

# create taylor
${taylor} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd72hr}-12_pair.hdf -s {'OZONE','PM2.5'} -n ${yyyymmddform}
for i in 'R1' 'R2' 'R3' 'R4' 'R5' 'R6' 'R7' 'R8' 'R9' 'R10'; do
  echo ${i}
  ${taylor} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd72hr}-12_pair.hdf -s {'OZONE','PM2.5'} -b True -e ${i} -n ${yyyymmddform}
done

# make GIFS
##########################################################################################################
for i in '5X' 'R1' 'R2' 'R3' 'R4' 'R5' 'R6' 'R7' 'R8' 'R9' 'R10'; do
  for j in 'OZONE' 'PM2P5'; do
    echo "${i}.${j}"
    convert -delay 40 -loop 0 ${yyyymmddform}.${i}.${j}.spbias.*.jpg ${yyyymmddform}.${i}.${j}.spbias.gif
  done
done

for i in '5X' 'R1' 'R2' 'R3' 'R4' 'R5' 'R6' 'R7' 'R8' 'R9' 'R10'; do
  for j in 'OZONE' 'PM2P5'; do
    echo "${i}.${j}"
    convert -delay 40 -loop 0 ${yyyymmddform}.${i}.${j}.sp.*.jpg ${yyyymmddform}.${i}.${j}.sp.gif
  done
done

for i in '5X' 'R1' 'R2' 'R3' 'R4' 'R5' 'R6' 'R7' 'R8' 'R9' 'R10'; do
  for j in 'OZONE' 'PM2P5'; do
    echo "${i}.${j}"
    convert -delay 40 -loop 0 ${yyyymmddform}.${i}.${j}.tyr.*.jpg ${yyyymmddform}.${i}.${j}.tyr.gif
  done
done
##########################################################################################################

# Transfer data back to aaqest
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#ls -t *.jpg | xargs -I {} --max-procs 4 scp {} patrickc@aaqest.arl.noaa.gov:/data/aqf2/testbed-bak/www/testbed/TCHAI/hkim/${yyyymmddform}/
#ls -t *.gif | xargs -I {} --max-procs 4 scp {} patrickc@aaqest.arl.noaa.gov:/data/aqf2/testbed-bak/www/testbed/TCHAI/hkim/${yyyymmddform}/


