#!/bin/bash
#SBATCH --job-name={job-name}          # nom du job
##SBATCH --partition=gpu_p2          # de-commente pour la partition gpu_p2
#SBATCH --nodes={nodes}                    # on demande un noeud 
#SBATCH --ntasks-per-node=1          # avec une tache par noeud (= nombre de GPU ici)
#SBATCH --gres=gpu:1                 # nombre de GPU (1/4 des GPU)
#SBATCH --cpus-per-task=10           # nombre de coeurs CPU par tache (1/4 du noeud 4-GPU)
##SBATCH --cpus-per-task=3           # nombre de coeurs CPU par tache (pour gpu_p2 : 1/8 du noeud 8-GPU)
# /!\ Attention, "multithread" fait reference à l'hyperthreading dans la terminologie Slurm
#SBATCH --hint=nomultithread         # hyperthreading desactive
#SBATCH --time={walltime}:00:00              # temps maximum d'execution demande (HH:MM:SS)
#SBATCH --output=log.slurm        # nom du fichier de sortie, %x is the jobname, %j the jobid
#SBATCH --error=%x.o%j           # nom du fichier d'erreur (ici commun avec la sortie)
#SBATCH --account=drd@v100              # To specify cpu accounting: <account> = echo $IDRPROJ
##SBATCH --qos=qos_gpu-dev           # Uncomment for job requiring less than 2 hours

# nettoyage des modules charges en interactif et herites par defaut
module purge
 
# chargement des modules
module load pytorch-gpu/py3/{pytorch_version}
 
# Avoid limited disk space of tmp dir
export TMPDIR=$JOBSCRATCH

# echo des commandes lancees
set -x
 
# execution du code
{code.previous_job}
srun /gpfsdswork/projects/rech/drd/uqx98co/applications/{path_to_lammps} -i lmp_input.in
{code.next_job}
