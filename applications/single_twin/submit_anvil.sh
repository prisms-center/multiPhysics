#!/bin/bash
#SBATCH -J sangwon
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH -A MSS160003
#SBATCH -p highmem
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --cpus-per-task=1
#SBATCH --mem=500G
#SBATCH -t 4:00:00
#SBATCH --mail-user=abhiks@umich.edu
#SBATCH --mail-type=all

module purge
module load gcc/11.2.0
module load openmpi/4.0.6
module list

pwd
date

# Creating target directory on scratch
locdir=${PWD##*/}
targpath=$SCRATCH
fulltargpath=$targpath/$locdir

echo "Running in $fulltargpath"

mkdir -p $fulltargpath

cp parameters.prm $fulltargpath
cp main $fulltargpath
cp LatentHardeningRatio.txt $fulltargpath
cp slipDirections.txt $fulltargpath
cp slipNormals.txt $fulltargpath
cp twinDirections.txt $fulltargpath
cp BCinfo.txt $fulltargpath
cp twinNormals.txt $fulltargpath
cp grainID.txt $fulltargpath
cp orientations.txt $fulltargpath

cd $fulltargpath

ulimit -s unlimited

echo "Start simulation"

mpirun --mca opal_common_ucx_opal_mem_hooks 1 -n $SLURM_NTASKS ./main parameters.prm

echo "Finished"
date
