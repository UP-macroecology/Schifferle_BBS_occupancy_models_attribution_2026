#!/bin/bash

# this script is written for the HPC cluster of the MacroEco and PENC lab at the University of Potsdam
# running with 78 cores and 750 GB of memory on Debian GNU/Linux and the workload manager slurm,
# paths must be adapted before job submission

#SBATCH --job-name=attrici_input_prep
#SBATCH --mail-type=ALL
#SBATCH --mail-user=schifferle1@uni-potsdam.de
#SBATCH --cpus-per-task=2
#SBATCH --time=1-00:00:00
#SBATCH --mem=20gb
#SBATCH --nodelist=ecoc9z

# Script:   1_2c_attrici_input_preps.sh
# Purpose:  Prepare input for climate data detrending with ATTRICI tool
# Inputs:   data/Counterfactual_env_data/ISIMIP_GSWP3_W5E5/attrici_detrending/input_files/gswp3-w5e5_obsclim_tas_global_daily_1991_2000.nc
#           data/Counterfactual_env_data/ISIMIP_GSWP3_W5E5/attrici_detrending/input_files/gswp3-w5e5_obsclim_tas_global_daily_2001_2010.nc
#           data/Counterfactual_env_data/ISIMIP_GSWP3_W5E5/attrici_detrending/input_files/gswp3-w5e5_obsclim_tas_global_daily_2011_2019.nc
#           data/Env_data/ISIMIP_GSWP3_W5E5/gswp3-w5e5_obsclim_tas_lon-126.0to-66.0lat24.0to50.0_daily_1991_2000.nc
#           data/Env_data/ISIMIP_GSWP3_W5E5/gswp3-w5e5_obsclim_tas_lon-126.0to-66.0lat24.0to50.0_daily_2001_2010.nc
#           data/Env_data/ISIMIP_GSWP3_W5E5/gswp3-w5e5_obsclim_tas_lon-126.0to-66.0lat24.0to50.0_daily_2011_2019.nc
#           data/Env_data/ISIMIP_GSWP3_W5E5/gswp3-w5e5_obsclim_pr_lat24.0to50.0lon-126.0to-66.0_daily_1991_2000.nc
#           data/Env_data/ISIMIP_GSWP3_W5E5/gswp3-w5e5_obsclim_pr_lat24.0to50.0lon-126.0to-66.0_daily_2001_2010.nc
#           data/Env_data/ISIMIP_GSWP3_W5E5/gswp3-w5e5_obsclim_pr_lat24.0to50.0lon-126.0to-66.0_daily_2011_2019.nc
#           data/Env_data/ISIMIP_GSWP3_W5E5/gswp3-w5e5_obsclim_tasmin_lat24.0to50.0lon-126.0to-66.0_daily_1991_2000.nc
#           data/Env_data/ISIMIP_GSWP3_W5E5/gswp3-w5e5_obsclim_tasmin_lat24.0to50.0lon-126.0to-66.0_daily_2001_2010.nc
#           data/Env_data/ISIMIP_GSWP3_W5E5/gswp3-w5e5_obsclim_tasmin_lat24.0to50.0lon-126.0to-66.0_daily_2011_2019.nc
#           data/Env_data/ISIMIP_GSWP3_W5E5/gswp3-w5e5_obsclim_tasmax_lat24.0to50.0lon-126.0to-66.0_daily_1991_2000.nc
#           data/Env_data/ISIMIP_GSWP3_W5E5/gswp3-w5e5_obsclim_tasmax_lat24.0to50.0lon-126.0to-66.0_daily_2001_2010.nc
#           data/Env_data/ISIMIP_GSWP3_W5E5/gswp3-w5e5_obsclim_tasmax_lat24.0to50.0lon-126.0to-66.0_daily_2011_2019.nc
# Outputs:  data/Counterfactual_env_data/ISIMIP_GSWP3_W5E5/attrici_detrending/input_files/gswp3-w5e5_obsclim_tas_1994_2019.nc
#           data/Counterfactual_env_data/ISIMIP_GSWP3_W5E5/attrici_detrending/input_files/gswp3-w5e5_gmt_raw_1994_2019.nc  
#           data/Counterfactual_env_data/ISIMIP_GSWP3_W5E5/attrici_detrending/input_files/gswp3-w5e5_ssa_gmt_1994_2019.nc
#           data/Env_data/ISIMIP_GSWP3_W5E5/US_tas_1994_2019.nc
#           data/Env_data/ISIMIP_GSWP3_W5E5/US_pr_1994_2019.nc
#           data/Env_data/ISIMIP_GSWP3_W5E5/US_tasmin_1994_2019.nc
#           data/Env_data/ISIMIP_GSWP3_W5E5/US_tasmax_1994_2019.nc
#           data/Counterfactual_env_data/ISIMIP_GSWP3_W5E5/attrici_detrending/input_files/US_tasrange_1994_2019.nc
#           data/Counterfactual_env_data/ISIMIP_GSWP3_W5E5/attrici_detrending/input_files/US_tasskew_1994_2019.nc
# Runs on:  HPC (NAS Potsdam)
# Steps:    1. Generate smoothed global mean temperature time series
#           1.1. Merge global temperature files over time
#           1.2. Calculate global mean temperature (GMT)
#           1.3. Smooth GMT with ATTRICI
#           2. Merge tas, pr, tasmin, tasmax for the conterminous USA over time
#           3. Preprocess temperature: tas, tasmin, tasmax to tasrange and tasskew with ATTRICI

# 1. Generate smoothed global mean temperature (GMT) curve:

# merge global tas files:

# input and output files:

base_dir="/mnt/ibb_share/zurell_transfer/Schifferle_BBS_occupancy_models_2023/data/Counterfactual_env_data/ISIMIP_GSWP3_W5E5/attrici_detrending/input_files"

# merge files by time:

echo "Merging tas files."

cdo -seldate,1994-01-01,2019-12-31 -mergetime "${base_dir}/gswp3-w5e5_obsclim_tas_global_daily_1991_2000.nc" "${base_dir}/gswp3-w5e5_obsclim_tas_global_daily_2001_2010.nc" "${base_dir}/gswp3-w5e5_obsclim_tas_global_daily_2011_2019.nc" "${base_dir}/gswp3-w5e5_obsclim_tas_1994_2019.nc" 
  
echo "Merging tas files done."

# Calculate global mean temperature time series:

cdo -L fldmean "${base_dir}/gswp3-w5e5_obsclim_tas_1994_2019.nc" "${base_dir}/gswp3-w5e5_gmt_raw_1994_2019.nc"

echo "Creating GMT file done."

# Smooth global mean temperature time series with ATTRICI:

# switch to ATTRICI directory
cd attrici

source env/bin/activate

attrici ssa \
  "${base_dir}/gswp3-w5e5_gmt_raw_1994_2019.nc" \
  "${base_dir}/gswp3-w5e5_ssa_gmt_1994_2019.nc"

echo "Smoothing done."


# 2. Merge pr, tasmin and tasmax files by time:

echo "Merge tas, pr, tasmin, tasmax."

clim_dir="/mnt/ibb_share/zurell_transfer/Schifferle_BBS_occupancy_models_2023/data/Env_data/ISIMIP_GSWP3_W5E5"

tas1="${clim_dir}/gswp3-w5e5_obsclim_tas_lon-126.0to-66.0lat24.0to50.0_daily_1991_2000.nc"
tas2="${clim_dir}/gswp3-w5e5_obsclim_tas_lon-126.0to-66.0lat24.0to50.0_daily_2001_2010.nc"
tas3="${clim_dir}/gswp3-w5e5_obsclim_tas_lon-126.0to-66.0lat24.0to50.0_daily_2011_2019.nc"
tas_out="${clim_dir}/US_tas_1994_2019.nc"

cdo -seldate,1994-01-01,2019-12-31 -mergetime "$tas1" "$tas2" "$tas3" "$tas_out"

pr1="${clim_dir}/gswp3-w5e5_obsclim_pr_lat24.0to50.0lon-126.0to-66.0_daily_1991_2000.nc"
pr2="${clim_dir}/gswp3-w5e5_obsclim_pr_lat24.0to50.0lon-126.0to-66.0_daily_2001_2010.nc"
pr3="${clim_dir}/gswp3-w5e5_obsclim_pr_lat24.0to50.0lon-126.0to-66.0_daily_2011_2019.nc"
pr_out="${clim_dir}/US_pr_1994_2019.nc"
cdo -seldate,1994-01-01,2019-12-31 -mergetime "$pr1" "$pr2" "$pr3" "$pr_out"

tasmin1="${clim_dir}/gswp3-w5e5_obsclim_tasmin_lat24.0to50.0lon-126.0to-66.0_daily_1991_2000.nc"
tasmin2="${clim_dir}/gswp3-w5e5_obsclim_tasmin_lat24.0to50.0lon-126.0to-66.0_daily_2001_2010.nc"
tasmin3="${clim_dir}/gswp3-w5e5_obsclim_tasmin_lat24.0to50.0lon-126.0to-66.0_daily_2011_2019.nc"
tasmin_out="${clim_dir}/US_tasmin_1994_2019.nc"
cdo -seldate,1994-01-01,2019-12-31 -mergetime "$tasmin1" "$tasmin2" "$tasmin3" "$tasmin_out"

tasmax1="${clim_dir}/gswp3-w5e5_obsclim_tasmax_lat24.0to50.0lon-126.0to-66.0_daily_1991_2000.nc"
tasmax2="${clim_dir}/gswp3-w5e5_obsclim_tasmax_lat24.0to50.0lon-126.0to-66.0_daily_2001_2010.nc"
tasmax3="${clim_dir}/gswp3-w5e5_obsclim_tasmax_lat24.0to50.0lon-126.0to-66.0_daily_2011_2019.nc"
tasmax_out="${clim_dir}/US_tasmax_1994_2019.nc"
cdo -seldate,1994-01-01,2019-12-31 -mergetime "$tasmax1" "$tasmax2" "$tasmax3" "$tasmax_out"

echo "Merging done."


# 3. Preprocess temperature: tas, tasmin, tasmax to tasrange and tasskew:

attrici preprocess-tas "$tas_out" "$tasmin_out" "$tasmax_out" "${base_dir}/US_tasrange_1994_2019.nc" "${base_dir}/US_tasskew_1994_2019.nc" 

echo "Tas preprocessing done."