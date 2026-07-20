#!/usr/bin/env bash

# this script is written for the HPC cluster of the MacroEco and PENC lab at the University of Potsdam
# running with 78 cores and 750 GB of memory on Debian GNU/Linux and the workload manager slurm,
# paths must be adapted before job submission

#SBATCH --job-name="attrici2"
#SBATCH --mail-type=ALL
#SBATCH --mail-user=schifferle1@uni-potsdam.de
#SBATCH --cpus-per-task=2
#SBATCH --export=ALL,OMP_PROC_BIND=TRUE
#SBATCH --ntasks=10
#SBATCH --time=23:00:00
#SBATCH --nodelist=ecoc9z
#SBATCH --mem=100gb


echo "Starting"

cd attrici

# load virtual environment if you don't queue with it activated:
source env/bin/activate

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

echo "Start detrending."

srun bash <<'EOF'

exec attrici \
     detrend \
     --gmt-file /import/ecoc9z/data-zurell/schifferle/BBS_occupancy_models_2023/data/Counterfactual_env_data/ISIMIP_GSWP3_W5E5/attrici_detrending/input_files/gswp3-w5e5_ssa_gmt_1995_2019.nc \
     --input-file /import/ecoc9z/data-zurell/schifferle/BBS_occupancy_models_2023/data/Counterfactual_env_data/ISIMIP_GSWP3_W5E5/attrici_detrending/input_files/US_pr_1901_2019.nc \
     --output-dir /import/ecoc9z/data-zurell/schifferle/BBS_occupancy_models_2023/data/Counterfactual_env_data/ISIMIP_GSWP3_W5E5/attrici_detrending/output/pr_detrended \
     --variable pr \
	   --start-date 1995-01-01 \
     --stop-date 2019-12-31 \
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


attrici merge-output /import/ecoc9z/data-zurell/schifferle/BBS_occupancy_models_2023/data/Counterfactual_env_data/ISIMIP_GSWP3_W5E5/attrici_detrending/output/pr_detrended/timeseries/pr  /import/ecoc9z/data-zurell/schifferle/BBS_occupancy_models_2023/data/Counterfactual_env_data/ISIMIP_GSWP3_W5E5/attrici_detrending/output/pr_detrended/US_pr_detrended_1901_2019.nc


echo "Output merging done."