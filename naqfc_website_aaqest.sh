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

#conda activate nco_stable for total NOX and PRECIP calculation
source /data/aqf2/barryb/anaconda2/bin/activate nco_stable
ncap2 -A -s "NOX=NO+NO2" $files_chem $files_nox
ncap2 -A -s "RT=RC+RN"   $files_met  $files_precip
source /data/aqf2/barryb/anaconda2/bin/deactivate

# pair the data
${pair} -f ${files_chem} -s 'OZONE' 'PM2.5' 'PM10' 'CO' 'NO' 'NO2' 'NO2Y' 
mv AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair.hdf AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair_chem.hdf
${pair} -f ${files_so2} -s 'SO2'
mv AIRNOW_CMAQ_${yyyymmdd}-12_${yyyymmdd48hr}-11_pair.hdf AIRNOW_CMAQ_${yyyymmdd}-12_${yyyymmdd48hr}-11_pair_so2.hdf
${pair} -f ${files_nox} -s 'NOX'
mv AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair.hdf AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair_nox.hdf
${pair} -f ${files_met}  -s 'BARPR' 'TEMP' 'WS' 'WD' 'SRAD' 'RHUM'
mv AIRNOW_CMAQ_${yyyymmdd}-12_${yyyymmdd48hr}-12_pair.hdf AIRNOW_CMAQ_${yyyymmdd}-12_${yyyymmdd48hr}-12_pair_met.hdf 
${pair} -f $files_precip -s 'PRECIP'
mv AIRNOW_CMAQ_${yyyymmdd}-12_${yyyymmdd48hr}-12_pair.hdf AIRNOW_CMAQ_${yyyymmdd}-12_${yyyymmdd48hr}-12_pair_precip.hdf

# spatial overlay plots
${spatial} -f ${files_chem}   -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair_chem.hdf    -s 'OZONE' 'PM2.5' 'PM10' 'CO' 'NO' 'NO2' 'NO2Y'       -n ${yyyymmddform}
${spatial} -f ${files_so2}    -p AIRNOW_CMAQ_${yyyymmdd}-12_${yyyymmdd48hr}-11_pair_so2.hdf     -s 'SO2'                                               -n ${yyyymmddform}
${spatial} -f ${files_nox}    -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair_nox.hdf     -s 'NOX'                                               -n ${yyyymmddform}
${spatial} -f ${files_met}    -p AIRNOW_CMAQ_${yyyymmdd}-12_${yyyymmdd48hr}-12_pair_met.hdf     -s 'BARPR' 'TEMP' 'WS' 'WD' 'SRAD'                     -n ${yyyymmddform}
${spatial} -f ${files_precip} -p AIRNOW_CMAQ_${yyyymmdd}-12_${yyyymmdd48hr}-12_pair_precip.hdf  -s 'PRECIP'                                            -n ${yyyymmddform}

for i in 'R1' 'R2' 'R3' 'R4' 'R5' 'R6' 'R7' 'R8' 'R9' 'R10'; do
  echo ${i}
  ${spatial} -f ${files_chem} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair_chem.hdf    -s  'OZONE' 'PM2.5' 'PM10' 'CO' 'NO' 'NO2' 'NO2Y'       -b True -e ${i} -n ${yyyymmddform}
  ${spatial} -f ${files_so2}  -p AIRNOW_CMAQ_${yyyymmdd}-12_${yyyymmdd48hr}-11_pair_so2.hdf     -s  'SO2'                                               -b True -e ${i} -n ${yyyymmddform}
  ${spatial} -f ${files_nox}  -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair_nox.hdf     -s  'NOX' 						-b True -e ${i} -n ${yyyymmddform}
  ${spatial} -f ${files_met}  -p AIRNOW_CMAQ_${yyyymmdd}-12_${yyyymmdd48hr}-12_pair_met.hdf     -s  'BARPR' 'TEMP' 'WS' 'WD' 'SRAD'                     -b True -e ${i} -n ${yyyymmddform}
  ${spatial} -f $files_precip -p AIRNOW_CMAQ_${yyyymmdd}-12_${yyyymmdd48hr}-12_pair_precip.hdf  -s  'PRECIP'                                            -b True -e ${i} -n ${yyyymmddform}
done

# spatial bias plots
${bias} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair_chem.hdf    -s 'OZONE' 'PM2.5' 'PM10' 'CO' 'NO' 'NO2' 'NO2Y'       -n ${yyyymmddform}
${bias} -p AIRNOW_CMAQ_${yyyymmdd}-12_${yyyymmdd48hr}-11_pair_so2.hdf     -s 'SO2'                                               -n ${yyyymmddform}
${bias} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair_nox.hdf     -s 'NOX'                                               -n ${yyyymmddform}
${bias} -p AIRNOW_CMAQ_${yyyymmdd}-12_${yyyymmdd48hr}-12_pair_met.hdf     -s 'BARPR' 'TEMP' 'WS' 'WD' 'SRAD' 'RHUM' 'PRECIP'     -n ${yyyymmddform}
${bias} -p AIRNOW_CMAQ_${yyyymmdd}-12_${yyyymmdd48hr}-12_pair_precip.hdf  -s 'PRECIP'                                            -n ${yyyymmddform}

for i in 'R1' 'R2' 'R3' 'R4' 'R5' 'R6' 'R7' 'R8' 'R9' 'R10'; do
  echo ${i}
  ${bias} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair_chem.hdf    -s 'OZONE' 'PM2.5' 'PM10' 'CO' 'NO' 'NO2' 'NO2Y'        -b True -e ${i} -n ${yyyymmddform}
  ${bias} -p AIRNOW_CMAQ_${yyyymmdd}-12_${yyyymmdd48hr}-11_pair_so2.hdf     -s 'SO2'                                                -b True -e ${i} -n ${yyyymmddform}
  ${bias} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair_nox.hdf     -s 'NOX'                                                -b True -e ${i} -n ${yyyymmddform}
  ${bias} -p AIRNOW_CMAQ_${yyyymmdd}-12_${yyyymmdd48hr}-12_pair_met.hdf     -s 'BARPR' 'TEMP' 'WS' 'WD' 'SRAD' 'RHUM' 'PRECIP'      -b True -e ${i} -n ${yyyymmddform}
  ${bias} -p AIRNOW_CMAQ_${yyyymmdd}-12_${yyyymmdd48hr}-12_pair_precip.hdf  -s 'PRECIP'                                             -b True -e ${i} -n ${yyyymmddform}
done

# create taylor
${taylor} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair_chem.hdf    -s 'OZONE' 'PM2.5' 'PM10' 'CO' 'NO' 'NO2' 'NO2Y'       -n ${yyyymmddform}
${taylor} -p AIRNOW_CMAQ_${yyyymmdd}-12_${yyyymmdd48hr}-11_pair_so2.hdf     -s 'SO2'                                               -n ${yyyymmddform}
${taylor} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair_nox.hdf     -s 'NOX'                                               -n ${yyyymmddform}
${taylor} -p AIRNOW_CMAQ_${yyyymmdd}-12_${yyyymmdd48hr}-12_pair_met.hdf     -s 'BARPR' 'TEMP' 'WS' 'WD' 'SRAD' 'RHUM' 'PRECIP'     -n ${yyyymmddform}
${taylor} -p AIRNOW_CMAQ_${yyyymmdd}-12_${yyyymmdd48hr}-12_pair_precip.hdf  -s 'PRECIP'                                            -n ${yyyymmddform}
for i in 'R1' 'R2' 'R3' 'R4' 'R5' 'R6' 'R7' 'R8' 'R9' 'R10'; do
  echo ${i}
  ${taylor} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair_chem.hdf   -s 'OZONE' 'PM2.5' 'PM10' 'CO' 'NO' 'NO2' 'NO2Y'       -b True -e ${i} -n ${yyyymmddform}
  ${taylor} -p AIRNOW_CMAQ_${yyyymmdd}-12_${yyyymmdd48hr}-11_pair_so2.hdf    -s 'SO2'                                               -b True -e ${i} -n ${yyyymmddform}
  ${taylor} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair_nox.hdf    -s 'NOX'                                               -b True -e ${i} -n ${yyyymmddform}
  ${taylor} -p AIRNOW_CMAQ_${yyyymmdd}-12_${yyyymmdd48hr}-12_pair_met.hdf    -s 'BARPR' 'TEMP' 'WS' 'WD' 'SRAD' 'RHUM' 'PRECIP'     -b True -e ${i} -n ${yyyymmddform}
  ${taylor} -p AIRNOW_CMAQ_${yyyymmdd}-12_${yyyymmdd48hr}-12_pair_precip.hdf -s 'PRECIP'                                            -b True -e ${i} -n ${yyyymmddform}
done

# make GIFS
##########################################################################################################
for i in '5X' 'R1' 'R2' 'R3' 'R4' 'R5' 'R6' 'R7' 'R8' 'R9' 'R10'; do
  for j in 'OZONE' 'PM2P5' 'PM10' 'CO' 'NO' 'NO2' 'NOX' 'NO2Y' 'BARPR' 'TEMP' 'WS' 'WD' 'SRAD' 'RHUM' 'PRECIP'; do
    echo "${i}.${j}"
    convert -delay 40 -loop 0 ${yyyymmddform}.${i}.${j}.spbias.*.jpg ${yyyymmddform}.${i}.${j}.spbias.gif
  done
done

for i in '5X' 'R1' 'R2' 'R3' 'R4' 'R5' 'R6' 'R7' 'R8' 'R9' 'R10'; do
  for j in 'OZONE' 'PM2P5' 'PM10' 'CO' 'NO' 'NO2' 'NOX' 'NO2Y' 'BARPR' 'TEMP' 'WS' 'WD' 'SRAD' 'PRECIP'; do
    echo "${i}.${j}"
    convert -delay 40 -loop 0 ${yyyymmddform}.${i}.${j}.sp.*.jpg ${yyyymmddform}.${i}.${j}.sp.gif
  done
done

for i in '5X' 'R1' 'R2' 'R3' 'R4' 'R5' 'R6' 'R7' 'R8' 'R9' 'R10'; do
  for j in 'OZONE' 'PM2P5' 'PM10' 'CO' 'NO' 'NO2' 'NOX' 'NO2Y' 'BARPR' 'TEMP' 'WS' 'WD' 'SRAD' 'RHUM' 'PRECIP'; do
    echo "${i}.${j}"
    convert -delay 40 -loop 0 ${yyyymmddform}.${i}.${j}.tyr.*.jpg ${yyyymmddform}.${i}.${j}.tyr.gif
  done
done
##########################################################################################################

mv *.jpg /data/aqf2/testbed-bak/www/testbed/TCHAI/hkim/${yyyymmddform}/
mv *.gif /data/aqf2/testbed-bak/www/testbed/TCHAI/hkim/${yyyymmddform}/

# Transfer data back to aaqest
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#ls -t *.jpg | xargs -I {} --max-procs 4 scp {} patrickc@aaqest.arl.noaa.gov:/data/aqf2/testbed-bak/www/testbed/TCHAI/hkim/${yyyymmddform}/
#ls -t *.gif | xargs -I {} --max-procs 4 scp {} patrickc@aaqest.arl.noaa.gov:/data/aqf2/testbed-bak/www/testbed/TCHAI/hkim/${yyyymmddform}/


