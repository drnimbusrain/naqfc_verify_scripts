#!/bin/bash -x


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
files=$yyyymmddform/"aqm.t12z.aconc-sfc.ncf"
files_met=$yyyymmddform/"aqm.t12z.metcro2d.ncf"

#naqfc verify scripts 
pair=/data/aqf/patrickc/naqfc_verify_scripts/01.verify_pair.py
stats=/data/aqf/patrickc/naqfc_verify_scripts/02.verify_stats.py
taylor=/data/aqf/patrickc/naqfc_verify_scripts/03.verify_taylor_plots.py
bias=/data/aqf/patrickc/naqfc_verify_scripts/04.verify_spatial_bias.py
spatial=/data/aqf/patrickc/naqfc_verify_scripts/05.verify_spatial_plots.py
box=/data/aqf/patrickc/naqfc_verify_scripts/06.verify_box_plots.py

# pair the data
${pair} -f ${files}     -s 'OZONE' 'PM2.5' 'PM10' 'CO' 'NO' 'NO2' 'NOX' 'NO2Y' 
mv AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair.hdf AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair_chem.hdf
${pair} -f ${files_met} -s 'BARPR' 'TEMP' 'WS' 'WD' 'SRAD' 'RHUM'
mv AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair.hdf AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair_met.hdf 

conda activate nco_stable
ncap2 -A -s "RT=RC+RN" $files_met $yyyymmddform/aqm.t12z.precip.ncf
deactivate

${pair} -f ${files_met} -s 'PRECIP'
mv AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair.hdf AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair_precip.hdf

# create taylor
${taylor} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair_chem.hdf    -s 'OZONE' 'PM2.5' 'PM10' 'CO' 'NO' 'NO2' 'NOX' 'NO2Y' -n ${yyyymmddform}
${taylor} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair_met.hdf     -s 'BARPR' 'TEMP' 'WS' 'WD' 'SRAD' 'RHUM' 'PRECIP'     -n ${yyyymmddform}
${taylor} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair_precip.hdf  -s 'PRECIP'                                            -n ${yyyymmddform}
for i in 'R1' 'R2' 'R3' 'R4' 'R5' 'R6' 'R7' 'R8' 'R9' 'R10'; do
  echo ${i}
  ${taylor} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair_chem.hdf -s 'OZONE' 'PM2.5' 'PM10' 'CO' 'NO' 'NO2' 'NOX' 'NO2Y' -b True -e ${i} -n ${yyyymmddform}
  ${taylor} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair_met.hdf  -s 'BARPR' 'TEMP' 'WS' 'WD' 'SRAD' 'RHUM' 'PRECIP'     -b True -e ${i} -n ${yyyymmddform}
done

# spatial bias plots
${bias} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair_chem.hdf -s 'OZONE' 'PM2.5' 'PM10' 'CO' 'NO' 'NO2' 'NOX' 'NO2Y' -n ${yyyymmddform}
${bias} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair_met.hdf  -s 'BARPR' 'TEMP' 'WS' 'WD' 'SRAD' 'RHUM' 'PRECIP'     -n ${yyyymmddform}
for i in 'R1' 'R2' 'R3' 'R4' 'R5' 'R6' 'R7' 'R8' 'R9' 'R10'; do
  echo ${i}
  ${bias} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair_chem.hdf -s 'OZONE' 'PM2.5' 'PM10' 'CO' 'NO' 'NO2' 'NOX' 'NO2Y'  -b True -e ${i} -n ${yyyymmddform}
  ${bias} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair_met.hdf  -s 'BARPR' 'TEMP' 'WS' 'WD' 'SRAD' 'RHUM' 'PRECIP'      -b True -e ${i} -n ${yyyymmddform}
done


# spatial overlay plots
${spatial} -f ${files} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair_chem.hdf -s 'OZONE' 'PM2.5' 'PM10' 'CO' 'NO' 'NO2' 'NOX' 'NO2Y' -n ${yyyymmddform}
${spatial} -f ${files} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair_met.hdf  -s 'BARPR' 'TEMP' 'WS' 'WD' 'SRAD' 'RHUM' 'PRECIP'     -n ${yyyymmddform}
for i in 'R1' 'R2' 'R3' 'R4' 'R5' 'R6' 'R7' 'R8' 'R9' 'R10'; do
  echo ${i}
  ${spatial} -f ${files} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair_chem.hdf -s  'OZONE' 'PM2.5' 'PM10' 'CO' 'NO' 'NO2' 'NOX' 'NO2Y' -b True -e ${i} -n ${yyyymmddform}
  ${spatial} -f ${files} -p AIRNOW_CMAQ_${yyyymmdd}-13_${yyyymmdd48hr}-12_pair_met.hdf  -s  'BARPR' 'TEMP' 'WS' 'WD' 'SRAD' 'RHUM' 'PRECIP'     -b True -e ${i} -n ${yyyymmddform}
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
  for j in 'OZONE' 'PM2P5' 'PM10' 'CO' 'NO' 'NO2' 'NOX' 'NO2Y' 'BARPR' 'TEMP' 'WS' 'WD' 'SRAD' 'RHUM' 'PRECIP'; do
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


