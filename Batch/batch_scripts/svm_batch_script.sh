##################################################
#Copy this file to your home directory and run it with qsub
#"qsub bondPricingBatch.sh" will start the job
#This script is intended to reserve 12 processors for Matlab worker processes
#In your batch script you can run "parpool('local',12)" to start 12 workers on a node
###################################################

#!/bin/bash
#PBS -N jSVM
#PBS -l nodes=1:ppn=20,mem=50g
#PBS -j oe
#PBS -V
#PBS -t 1-1600
## 150


cd $PBS_O_WORKDIR

#Create a job specific temp directory
mkdir -p ~/BondPricing/Julia/log_files/SVM/$PBS_JOBID
export JULIAWORKDIR=~/BondPricing/Julia/log_files/SVM/$PBS_JOBID

# Load Python and Julia Modules
source ~/BondPricing/Julia/Batch/module_loader.sh

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
skip_julia=1
coupon_pos=999 # $PBS_ARRAYID  # 99
skip_c=1
skip_sol=1
skip_all_eqfd=1

echo i. Number of processors/cores: $PBS_NUM_PPN
echo ii. Memory: $job_mem
echo iii. Parameter Combination: $param_comb
echo iv. Skip Julia Calculations: $skip_julia
echo v. Coupon Position: $coupon_pos
echo vi. Rerun Finite Difference: $skip_c
echo vii. Rerun 1st Part Calculations: $skip_sol
echo viii. Rerun Equity Finite Differences Method: $skip_sol

julia svm_batch_script.jl $param_comb $skip_julia $coupon_pos $skip_c $skip_sol $skip_all_eqfd >> $JULIAWORKDIR/batch.log 2>&1
