#!/bin/bash
#SBATCH -N 1
#SBATCH -J NT_test
#SBATCH -n 1 
#SBATCH --mem=16G
#SBATCH --output=log/NeuTran-%j.out
#SBATCH --error=log/NeuTran-%j.err
#SBATCH -t 500:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=694076215@qq.com

##echo "begin time"; date
##echo "job run on :"; hostname
##echo "with job id:" $SLURM_JOBID

module load valgrind
#接收传递的输入文件路径
input_file=$1
##do your job
cd $SLURM_SUBMIT_DIR
mkdir $PWD/$SLURM_JOBID
mkdir $PWD/$SLURM_JOBID/data
mkdir $PWD/$SLURM_JOBID/fig
cp "$input_file" $PWD/$SLURM_JOBID
mpiexec.hydra -np 1 valgrind ./edit "$input_file" $PWD/$SLURM_JOBID/data/ --tool=memcheck --leak-check=full
#mv $PWD/Outputfile/* $PWD/$SLURM_JOBID/data
mv $PWD/log/NeuTran-$SLURM_JOBID* $PWD/$SLURM_JOBID
mv $PWD/$SLURM_JOBID Outputfile
#mail -s "Job Output" 694076215@qq.com < output.txt