#!/bin/bash
#SBATCH -M mindwell
#SBATCH -p batch_graniterapids
#SBATCH -N 1
#SBATCH -J NT_test
#SBATCH -n 1
#SBATCH --mem=2G
#SBATCH --output=log/NeuTran-%j.out
#SBATCH --error=log/NeuTran-%j.err
#SBATCH -t 3-00:00:00
#SBATCH -A lp_opt_plasma_edge_conf
#SBATCH --mail-type=ALL
#SBATCH --mail-user=694076215@qq.com

# 接收传递的输入文件路径
input_file=$1

cd $SLURM_SUBMIT_DIR

mkdir -p log
mkdir -p Outputfile

export LD_LIBRARY_PATH=/data/leuven/379/vsc37950/nt/opt/hdf5/lib:/vsc-hard-mounts/leuven-apps/rocky9/skylake/2025b/software/GCCcore/14.3.0/lib64:$LD_LIBRARY_PATH

workdir=$PWD/$SLURM_JOBID

mkdir -p $workdir
mkdir -p $workdir/data
mkdir -p $workdir/fig

cp "$input_file" $workdir/

srun -n 1 ./bin/edit "$input_file" $workdir/data/
run_status=$?

mv $PWD/log/NeuTran-$SLURM_JOBID* $workdir/ 2>/dev/null
mv $workdir Outputfile/
exit $run_status
