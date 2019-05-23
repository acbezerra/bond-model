##################################################
#Copy this file to your home directory and run it with qsub
#"qsub bondPricingBatch.sh" will start the job
###################################################

#!/bin/bash
#PBS -N jSVM
#PBS -l nodes=1:ppn=20,mem=20g
#PBS -j oe
#PBS -V
#PBS -t 1-1


cd $PBS_O_WORKDIR

#Create a job specific temp directory
mkdir -p ~/BondPricing/Julia/log_files/CVM/$PBS_JOBID
export JULIAWORKDIR=~/BondPricing/Julia/log_files/CVM/$PBS_JOBID

# Load Python and Julia Modules
module load python/intelpython3
export PATH=/home/artur/BondPricing/.julia/julia-1.1.0/bin:$PATH
export LD_LIBRARY_PATH=/home/artur/BondPricing/.julia/julia-1.1.0/bin:$LD_LIBRARY_PATH

# $PBS_NUM_PPN gives the number of processors to be used in each node.
# $PBS_ARRAYID gives the parameter combination # position of the coupon in the coupon grid -> Each batch gets one value of c!
echo $PBS_O_WORKDIR
echo $JULIAWORKDIR
echo $PBS_NODEFILE
echo $PBS_ARRAYID

# SYS Arguments:
# i. Number of Processors/Cores;
# ii. Parameter Combination (previously, $PBS_ARRAYID. Now: just enter combination number -> int);
# iii. Coupon Value (if > # of coupon values in coupon grid, run for all coupon values!);
# iv. Whether to skip calculations for a given coupon value if results already exist;
#     (if skip_sol = 0, run everything. Else, run only Finite Difference Calculations);
# v. Whether to skip 1st Part calculations if results already exist;
# vi. Whether to skip Equity Finite Differences if results already exist. Do not skip if Part I is run.
param_comb=$PBS_ARRAYID  # 1  # $PBS_ARRAYID

echo i. Number of processors/cores: $PBS_NUM_PPN
echo ii. Memory: $job_mem
echo iii. Maturity Parameter Combination: $param_comb

julia cvm_batch_script.jl $param_comb >> $JULIAWORKDIR/batch.log 2>&1

