#!/bin/bash
#SBATCH --job-name=protein_vac
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -p fgpu
#SBATCH --gres=gpu:a4000:1
#SBATCH --time=02:00:00
#SBATCH --output=vac.log

module load cuda12.2

echo "CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES"
nvidia-smi

./run_md.sh -g DMPVDPDNEAYEMPSEEGYQDYEPEA -s Gen