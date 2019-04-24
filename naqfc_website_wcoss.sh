#!/bin/bash -x

# load monet modules
module purge
module use -a /gpfs/hps3/emc/naqfc/noscrub/Barry.Baker/modulefiles
module load anaconda3/latest

export OMP_NUM_THREADS=2
# get the current date
date
yyyymmdd=$(date -d "+0 days" +%Y-%m-%d)
yyyymmddform=$(date -d "+0 days" +%Y%m%d)
yyyymmdd72hr=$(date -d "+3 days" +%Y-%m-%d)

# go to data directory
data_dir=/data/aqf/patrickc/CMAQ_DATA/cmaq5.0.2-72hr
#data_dir=/gpfs/hps2/ptmp/Patrick.C.Campbell/fengsha/gfs.$yyyymmdd/00
cd $data_dir

# netcdf files
files="aqm.${yyyymmddform}.t12z.aconc-pm25.ncf"

#naqfc verify scripts 
pair=/gpfs/hps3/emc/naqfc/noscrub/Patrick.C.Campbell/naqfc_verify_scripts/01.verify_pair.py
stats=/gpfs/hps3/emc/naqfc/noscrub/Patrick.C.Campbell/naqfc_verify_scripts/02.verify_stats.py
taylor=/gpfs/hps3/emc/naqfc/noscrub/Patrick.C.Campbell/naqfc_verify_scripts/03.verify_taylor_plots.py
bias=/gpfs/hps3/emc/naqfc/noscrub/Patrick.C.Campbell//naqfc_verify_scripts/04.verify_spatial_bias.py
spatial=/gpfs/hps3/emc/naqfc/noscrub/Patrick.C.Campbell/naqfc_verify_scripts/05.verify_spatial_plots.py
box=/gpfs/hps3/emc/naqfc/noscrub/Patrick.C.Campbell/naqfc_verify_scripts/06.verify_box_plots.py

# pair the data
${pair} -f ${files} 

# make spatial overlay plots for all regions
${spatial} -f ${files} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd72hr}-12_pair.hdf -s {'OZONE','PM2.5'} 
for i in 'R1' 'R2' 'R3' 'R4' 'R5' 'R6' 'R7' 'R8' 'R9' 'R10'; do
  echo ${i}
  ${spatial} -f ${files} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd72hr}-12_pair.hdf -s {'OZONE','PM2.5'} -b True -e ${i}
done

# create spatial bias plots
${bias} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd72hr}-12_pair.hdf -s {'OZONE','PM2.5'} 
for i in 'R1' 'R2' 'R3' 'R4' 'R5' 'R6' 'R7' 'R8' 'R9' 'R10'; do
  echo ${i}
  ${bias} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd72hr}-12_pair.hdf -s {'OZONE','PM2.5'}  -b True -e ${i}
done 

# create box plot
#${box} -p AIRNOW_CMAQ_2018-12-01-13_2019-01-31-12_pair.hdf -s {'OZONE','PM2.5'} -r True

# create fv3_aeronet_taylor
${taylor} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd72hr}-12_pair.hdf -s {'OZONE','PM2.5'}
for i in 'R1' 'R2' 'R3' 'R4' 'R5' 'R6' 'R7' 'R8' 'R9' 'R10'; do
  echo ${i}
  ${taylor} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd72hr}-12_pair.hdf -s {'OZONE','PM2.5'} -b True -e ${i}
done

# make GIFS
##########################################################################################################
for i in 'domain' 'R1' 'R2' 'R3' 'R4' 'R5' 'R6' 'R7' 'R8' 'R9' 'R10'; do
  for j in 'OZONE' 'PM2.5'; do
    echo "${i}.${j}"
    convert -delay 40 -loop 0 CMAQ_AIRNOW.${j}.${i}.sp.*.jpg CMAQ_AIRNOW.${j}.${i}.sp.gif
  done
done

for i in 'domain' 'R1' 'R2' 'R3' 'R4' 'R5' 'R6' 'R7' 'R8' 'R9' 'R10'; do
  for j in 'OZONE' 'PM2.5'; do
    echo "${i}.${j}"
    convert -delay 40 -loop 0 CMAQ_AIRNOW.${j}.${i}.sb.*.jpg CMAQ_AIRNOW.${j}.${i}.sb.gif
  done
done

for i in 'domain' 'R1' 'R2' 'R3' 'R4' 'R5' 'R6' 'R7' 'R8' 'R9' 'R10'; do
  for j in 'OZONE' 'PM2.5'; do
    echo "${i}.${j}"
    convert -delay 40 -loop 0 CMAQ_AIRNOW.${j}.${i}.tyr.*.jpg CMAQ_AIRNOW.${j}.${i}.tyr.gif
  done
done
##########################################################################################################

# Transfer data back to aaqest
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
ls -t *.jpg | xargs -I {} --max-procs 4 scp {} patrickc@aaqest.arl.noaa.gov:/data/aqf2/testbed-bak/www/testbed/TCHAI/hkim/${yyyymmddform}/
ls -t *.gif | xargs -I {} --max-procs 4 scp {} patrickc@aaqest.arl.noaa.gov:/data/aqf2/testbed-bak/www/testbed/TCHAI/hkim/${yyyymmddform}/


