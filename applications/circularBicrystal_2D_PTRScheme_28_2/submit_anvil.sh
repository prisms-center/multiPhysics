#!/bin/bash

#SBATCH -J dmontiel_job
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH -A MSS160003
#SBATCH -p wholenode
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH -t 3:00:00
#SBATCH --mail-user=dmontiel@umich.edu
#SBATCH --mail-type=all

module list
pwd
date

#Creating target directory on scratch and copying files into it
locdir=${PWD##*/}
targpath=$SCRATCH
fulltargpath=$targpath/$locdir

echo $fulltargpath

mkdir $fulltargpath
cp parameters.prm $fulltargpath
cp main $fulltargpath
cp LatentHardeningRatio.txt $fulltargpath
cp slipDirections.txt $fulltargpath
cp slipNormals.txt $fulltargpath
cp twinDirections.txt $fulltargpath
cp twinNormals.txt $fulltargpath
cp grainID.txt $fulltargpath
cp orientations.txt $fulltargpath

cd $fulltargpath
srun -n $SLURM_NTASKS ./main parameters.prm


