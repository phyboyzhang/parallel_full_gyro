#!/bin/bash
##SBATCH -p public
##SBATCH -n 128
#SBATCH -N 1
#SBATCH --ntasks-per-core 1
#SBATCH --ntasks-per-node 16
#SBATCH -t 24:00:00
#SBATCH --output job_output.text
#SBATCH --error err_job.text
#SBATCH --job-name com_fulgy


RUN_DIR=/home/qmlu/zsx163/parallel_full_gyro/run
BUILD_DIR=/home/qmlu/zsx163/parallel_full_gyro/build/bin/
EXEC=$BUILD_DIR/test_orbit_comparison

## Here you have to load every module use for compilation :
##module load my_module_1
##module load my_module_2

cd $RUN_DIR
#mpirun --bind-to core $EXEC $RUN_DIR/dksim4d_polar_multi_mu_ng.nml
mpirun $EXEC 

