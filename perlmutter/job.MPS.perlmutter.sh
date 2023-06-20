#!/bin/bash -l

# example for a 2 node, 8 GPU run


#SBATCH -A m4093_g
#SBATCH -C gpu
#SBATCH -q regular
#SBATCH -t 00:10:00
#SBATCH -N 2
#SBATCH -n 16
#SBATCH --ntasks-per-node=8
#SBATCH -c 16
#SBATCH -G 8
##SBATCH --gpus-per-node=4
#SBATCH -e MURaM.err
#SBATCH -o MURaM.out


module purge
module load PrgEnv-nvidia
module load craype-x86-milan
module load cray-fftw/3.3.8.13
module load cudatoolkit
module list

export HEFFTE_HOME=/pscratch/sd/c/cmille73/nvhpc_2111_dam/heffte/
export PATH=$HEFFTE_HOME/bin:$PATH
export LD_LIBRARY_PATH=$HEFFTE_HOME/lib:$LD_LIBRARY_PATH

export MURaM_FFTW_THREADS=1
export MPICH_GPU_SUPPORT_ENABLED=1
export CRAY_ACCEL_TARGET=nvidia80
### Source the test environment
echo "job starts"
echo $cwd
echo "PATH: "
echo $PATH
echo "LD_LIBRARY_PATH: "
echo $LD_LIBRARY_PATH

./clean

ulimit -s unlimited
nvidia-smi

export CUDA_LAUNCH_BLOCKING=0
#export UCX_TLS=rc,sm,cuda_copy,cuda_ipc,gdr_copy
#export OMPI_MCA_pml=ucx
#export OMPI_MCA_btl=self,vader,tcp,smcuda #openib
#export UCX_RNDV_SCHEME=get_zcopy
#export UCX_RNDV_THRESH=0
#export UCX_MAX_RNDV_RAILS=1
#export UCX_MEMTYPE_CACHE=n

srun ./perlMPSwrapper.sh &

srun ./mhd3d.x > MURaM.out

