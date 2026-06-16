#!/bin/bash
#SBATCH -M wice
#SBATCH -N 1
#SBATCH -J NT_test
#SBATCH -n 1
#SBATCH --mem=2G
#SBATCH --output=log/NeuTran-%j.out
#SBATCH --error=log/NeuTran-%j.err
#SBATCH -t 3-00:00:00
#SBATCH -A lp_opt_plasma_edge_conf     # 计费/项目账号（你的 VSC project）
#SBATCH --mail-type=ALL
#SBATCH --mail-user=694076215@qq.com

##echo "begin time"; date
##echo "job run on :"; hostname
##echo "with job id:" $SLURM_JOBID


#接收传递的输入文件路径
input_file=$1
##do your job
cd $SLURM_SUBMIT_DIR
export LD_LIBRARY_PATH=/data/leuven/379/vsc37950/nt/opt/hdf5/lib:/vsc-hard-mounts/leuven-apps/rocky9/skylake/2025b/software/GCCcore/14.3.0/lib64
mkdir $PWD/$SLURM_JOBID
mkdir $PWD/$SLURM_JOBID/data
mkdir $PWD/$SLURM_JOBID/fig
cp "$input_file" $PWD/$SLURM_JOBID
srun -n 1 ./bin/edit "$input_file" $PWD/$SLURM_JOBID/data/
#mv $PWD/Outputfile/* $PWD/$SLURM_JOBID/data
mv $PWD/log/NeuTran-$SLURM_JOBID* $PWD/$SLURM_JOBID
mv $PWD/$SLURM_JOBID Outputfile
#mail -s "Job Output" 694076215@qq.com < output.txt
