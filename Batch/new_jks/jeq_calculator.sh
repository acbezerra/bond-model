##################################################
#Copy this file to your home directory and run it with qsub
#"qsub bondPricingBatch.sh" will start the job
#This script is intended to reserve 12 processors for Matlab worker processes
#In your batch script you can run "parpool('local',12)" to start 12 workers on a node
###################################################

#!/bin/bash
#PBS -N jSVM
#PBS -l nodes=1:ppn=12,mem=40g
#PBS -j oe
#PBS -V
#PBS -t 1-1 #70
## # 400   # 56


cd $PBS_O_WORKDIR

#Create a job specific temp directory
mkdir -p ~/BondPricing/bond-model/excluded/log_files/JEQ/$PBS_JOBID
export JULIAWORKDIR=~/BondPricing/bond-model/excluded/log_files/JEQ/$PBS_JOBID

# Load Python and Julia Modules
source ~/BondPricing/bond-model/Batch/module_loader.sh

# $PBS_NUM_PPN gives the number of processors to be used in each node.
# $PBS_ARRAYID gives the position in the vector of measures of safe firm
echo $PBS_O_WORKDIR
echo $JULIAWORKDIR
echo $PBS_NODEFILE
echo $PBS_ARRAYID

# SYS Arguments:
fixed_mu_s=0
pair_num=$PBS_ARRAYID #7
rerun_fi=1
rerun_misrep=1
run_pool=1
run_sep=1

echo i. Number of processors/cores: $PBS_NUM_PPN
#echo ii. Memory: $job_mem
echo ii. Position in the Safe Firm Measure array: $mu_s_pos
echo iii. Rerun Full Information Results: $rerun_fi
echo iv. Rerun Misrepresentation Results: $rerun_misrep
echo v. Run Pooling Equilibrium: $run_pool
echo vi. Run Separating Equilibrium: $run_sep

julia jeq_calculator.jl $fixed_mu_s $pair_num $rerun_fi $rerun_misrep $run_pool $run_sep >> $JULIAWORKDIR/batch.log 2>&1
