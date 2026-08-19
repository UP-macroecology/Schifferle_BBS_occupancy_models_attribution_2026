#!/usr/bin/env bash

# this script is written for the HPC cluster of the MacroEco and PENC lab at the University of Potsdam
# running with 78 cores and 750 GB of memory on Debian GNU/Linux and the workload manager slurm,
# paths must be adapted before job submission

#SBATCH --job-name="attrici_tasrange"
#SBATCH --mail-type=ALL
#SBATCH --mail-user=schifferle1@uni-potsdam.de
#SBATCH --cpus-per-task=2
#SBATCH --export=ALL,OMP_PROC_BIND=TRUE
#SBATCH --ntasks=10
#SBATCH --time=1-00:00:00
#SBATCH --nodelist=ecoc9
#SBATCH --mem=150gb

# Script:   1_2d_attrici_US_tasrange.sh
# Purpose:  Detrend tasrange from 1994 onward with ATTRICI tool (from tasrange and tasskew tasmin and tasmax are derived later)
# Inputs:   data/Counterfactual_env_data/ISIMIP_GSWP3_W5E5/attrici_detrending/input_files/gswp3-w5e5_ssa_gmt_1994_2019.nc
#           data/Env_data/ISIMIP_GSWP3_W5E5/US_tasrange_1994_2019.nc
#           data/Counterfactual_env_data/ISIMIP_GSWP3_W5E5/attrici_detrending/input_files/US_mask.nc
# Outputs:  data/Counterfactual_env_data/ISIMIP_GSWP3_W5E5/attrici_detrending/output/tasrange_detrended/US_tasrange_detrended_1994_2019.nc
# Runs on:  HPC (NAS Potsdam)

echo "Starting"

cd attrici

# load virtual environment if you don't queue with it activated:
source env/bin/activate

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

srun bash <<'EOF'

exec attrici \
     detrend \
     --gmt-file /import/ecoc9z/data-zurell/schifferle/BBS_occupancy_models_2023/data/Counterfactual_env_data/ISIMIP_GSWP3_W5E5/attrici_detrending/input_files/gswp3-w5e5_ssa_gmt_1994_2019.nc \
     --input-file /import/ecoc9z/data-zurell/schifferle/BBS_occupancy_models_2023/data/Counterfactual_env_data/ISIMIP_GSWP3_W5E5/attrici_detrending/input_files/US_tasrange_1994_2019.nc \
     --output-dir /import/ecoc9z/data-zurell/schifferle/BBS_occupancy_models_2023/data/Counterfactual_env_data/ISIMIP_GSWP3_W5E5/attrici_detrending/output/tasrange_detrended \
     --variable tasrange \
	 --mask-file /import/ecoc9z/data-zurell/schifferle/BBS_occupancy_models_2023/data/Counterfactual_env_data/ISIMIP_GSWP3_W5E5/attrici_detrending/input_files/US_mask.nc \
     --report-variables y cfact logp \
   	 --full-extrapolation \
     --task-id "$SLURM_PROCID" \
     --task-count "$SLURM_NTASKS" \
	   --seed 123
EOF

echo "Detrending done."

# Merge output:

echo "Merge output."

attrici merge-output /import/ecoc9z/data-zurell/schifferle/BBS_occupancy_models_2023/data/Counterfactual_env_data/ISIMIP_GSWP3_W5E5/attrici_detrending/output/tasrange_detrended/timeseries/tasrange  /import/ecoc9z/data-zurell/schifferle/BBS_occupancy_models_2023/data/Counterfactual_env_data/ISIMIP_GSWP3_W5E5/attrici_detrending/output/tasrange_detrended/US_tasrange_detrended_1994_2019.nc


echo "Output merging done."