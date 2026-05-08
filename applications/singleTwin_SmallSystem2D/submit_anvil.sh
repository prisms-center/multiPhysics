#! /bin/bash

#SBATCH -J prisms_mp_singletwin
#SBATCH -o prisms_mp_%j_stdout.txt
#SBATCH -e prisms_mp_%j_stderr.txt
#SBATCH -A MSS160003
#SBATCH -p wholenode
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH -t 2:00:00

# NOTE: this job is currently expected to FAIL on wholenode
# due to lack of memory. A memory error during initialization
# is causing OOM. Running once more with a 1 hr runtime to
# reproduce the error with additional debugging/logging statements.
# The actual job should take much longer than an hour, and runs
# fine on highmem with 128 cores/1 node.

source /anvil/projects/x-mss160003/dii_w_petsc/dealii_candi/configuration/enable.sh

module list
pwd
date

fulltargetpath="${SCRATCH}/prisms_mp_singletwin_compare_matlab_seed16"

echo "Target path = ${fulltargetpath}"
mkdir -p "${fulltargetpath}"

cp parameters_cp.prm ${fulltargetpath}
cp parameters_pf.prm ${fulltargetpath}
cp LatentHardeningRatio.txt ${fulltargetpath}
cp orientations_twin_single.txt ${fulltargetpath}
cp slipDirections.txt ${fulltargetpath}
cp slipNormals.txt ${fulltargetpath}
cp twinDirections.txt ${fulltargetpath}
cp twinNormals.txt ${fulltargetpath}
cp grainID_single_32x32x32.txt ${fulltargetpath}
cp main ${fulltargetpath}

cd ${fulltargetpath}
mkdir -p results_cp
mkdir -p results_pf
mpirun -n $SLURM_NTASKS ./main parameters_cp.prm parameters_pf.prm

echo "Done, at:"
date

module purge

