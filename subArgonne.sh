#PBS -l nodes=1:ppn=8,mem=1gb,walltime=00:20:00
#PBS -N simu
#PBS -j oe
#PBS -t 1

cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=8

Rscript proposed.R $PBS_ARRAYID
