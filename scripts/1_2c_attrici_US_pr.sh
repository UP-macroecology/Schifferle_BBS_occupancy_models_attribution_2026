#!/usr/bin/env bash
#SBATCH --job-name="attrici1"
#SBATCH --mail-type=ALL
#SBATCH --mail-user=schifferle1@uni-potsdam.de
#SBATCH --cpus-per-task=2
#SBATCH --export=ALL,OMP_PROC_BIND=TRUE
#SBATCH --ntasks=4
#SBATCH --time=01:00:00
#SBATCH --nodelist=ecoc9z

# load necessary modules/packages here if you don't queue with them loaded
# e.g.: module purge; module load ...
#   or: spack load ...

# load virtual environment if you don't queue with it activated:
# e.e.: source venv/bin/activate

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

srun bash <<'EOF'

exec attrici \
     detrend \
     --gmt-file /mnt/ibb_share/zurell_transfer/Schifferle_BBS_occupancy_models_2023/data/Counterfactual_env_data/ISIMIP_GSWP3_W5E5/attrici_detrending/input_files/20CRv3-ERA5_germany_ssa_gmt.nc \
     --input-file /mnt/ibb_share/zurell_transfer/Schifferle_BBS_occupancy_models_2023/data/Counterfactual_env_data/ISIMIP_GSWP3_W5E5/attrici_detrending/input_files/US_pr_1901_2019.nc \
     --output-dir /mnt/ibb_share/zurell_transfer/Schifferle_BBS_occupancy_models_2023/data/Counterfactual_env_data/ISIMIP_GSWP3_W5E5/attrici_detrending/output/pr_detrended \
     --variable pr \
	 --start-date 1995-01-01 \
     --stop-date 2019-12-31 \
	 --mask-file /mnt/ibb_share/zurell_transfer/Schifferle_BBS_occupancy_models_2023/data/Counterfactual_env_data/ISIMIP_GSWP3_W5E5/attrici_detrending/input_files/US_mask3.nc \
     --report-variables y cfact logp \
     --overwrite \
     --task-id "$SLURM_PROCID" \
     --task-count "$SLURM_NTASKS"
EOF