#!/bin/bash

#SBATCH -J dmontiel_job
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH -A prisms_project1
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH -t 12:00:00
#SBATCH --mail-user=dmontiel@umich.edu
#SBATCH --mail-type=all

#~/env2.sh
module list
pwd
date

#Creating target directory on scratch and copying files into it
locdir=${PWD##*/}
targpath=/scratch/prisms_project_root/prisms_project1/dmontiel
fulltargpath=$targpath/$locdir

echo $fulltargpath

mkdir $fulltargpath
cp parameters_pf.prm $fulltargpath
cp prm_comp_3.prm $fulltargpath
cp main $fulltargpath
cp BCinfo.txt $fulltargpath
cp LatentHardeningRatio.txt $fulltargpath
cp slipDirections.txt $fulltargpath
cp slipNormals.txt $fulltargpath
cp twinDirections.txt $fulltargpath
cp twinNormals.txt $fulltargpath
cp grainID_single_64x64x64.txt $fulltargpath
cp orientations_twin_single.txt $fulltargpath

cd $fulltargpath
srun ./main prm_comp_3.prm
