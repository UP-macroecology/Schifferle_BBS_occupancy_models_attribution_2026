#!/usr/bin/env bash

# this script is written for the HPC cluster of the MacroEco and PENC lab at the University of Potsdam
# running with 78 cores and 750 GB of memory on Debian GNU/Linux and the workload manager slurm,
# paths must be adapted before job submission

#SBATCH --job-name="attrici_tasskew"
#SBATCH --mail-type=ALL
#SBATCH --mail-user=schifferle1@uni-potsdam.de
#SBATCH --cpus-per-task=5
#SBATCH --export=ALL,OMP_PROC_BIND=TRUE
#SBATCH --ntasks=1
#SBATCH --time=1-00:00:00
#SBATCH --nodelist=ecoc9z

# load necessary modules/packages here if you don't queue with them loaded
# e.g.: module purge; module load ...
#   or: spack load ...

# load virtual environment if you don't queue with it activated:
source env/bin/activate

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

srun bash <<'EOF'

exec attrici \
     detrend \
     --gmt-file /mnt/ibb_share/zurell_transfer/Schifferle_BBS_occupancy_models_2023/data/Counterfactual_env_data/ISIMIP_GSWP3_W5E5/attrici_detrending/input_files/20CRv3-ERA5_germany_ssa_gmt.nc \
     --input-file /mnt/ibb_share/zurell_transfer/Schifferle_BBS_occupancy_models_2023/data/Counterfactual_env_data/ISIMIP_GSWP3_W5E5/attrici_detrending/input_files/US_tasskew_1901_2019.nc \
     --output-dir /mnt/ibb_share/zurell_transfer/Schifferle_BBS_occupancy_models_2023/data/Counterfactual_env_data/ISIMIP_GSWP3_W5E5/attrici_detrending/output/tasskew_detrended \
     --variable tasskew \
	 --start-date 1995-01-01 \
     --stop-date 2019-12-31 \
	 --mask-file /mnt/ibb_share/zurell_transfer/Schifferle_BBS_occupancy_models_2023/data/Counterfactual_env_data/ISIMIP_GSWP3_W5E5/attrici_detrending/input_files/US_mask.nc \
     --report-variables y cfact logp \
     --overwrite \
     --task-id "$SLURM_PROCID" \
     --task-count "$SLURM_NTASKS"
EOF
