#!/bin/bash
#PBS -N Preindustrial
#PBS -q serial
#PBS -l walltime=2:30:00
#PBS -l mem=7500mb
#PBS -l nodes=1:ppn=24
#PBS -m abe


# Make folders as shortcuts

DATADIR="/home/obj30718/beegfs/Data"
CODE_FOLD="/home/obj30718/SensitivityAnalysisLPJ/"
RESULT_FOLD="/home/obj30718/beegfs/Results_sensi/"

# Export library Path

export LD_LIBRARY_PATH=/home/obj30718/R/x86_64-pc-linux-gnu-library/3.3

# -------------------------------------

# Copy stuff from your folder to the computing nodes

cp -r $CODE_FOLD $TMPDIR

cp -r $DATADIR $TMPDIR/SensitivityAnalysisLPJ/LPJrunTest/

cd $TMPDIR/SensitivityAnalysisLPJ/

# Let the script run ( info : must be executable: chmod +x filename)
# export OMP_NUM_THREADS=1

Rscript ./R_scripts/_script_to_execute_.R

cp -r $TMPDIR/SensitivityAnalysisLPJ/LPJrunTest/Results/_script_to_execute_.rds $RESULT_FOLD

#cp -r $TMPDIR/SensitivityAnalysisLPJ/LPJrunTest/runDirectory1/guess.log $RESULT_FOLD
#cp -r $TMPDIR/SensitivityAnalysisLPJ/LPJrunTest/runDirectory1/main_new.ins $RESULT_FOLD
