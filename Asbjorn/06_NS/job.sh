#!/bin/sh 
### General options 
### -- specify queue -- 
#BSUB -q hpc
### -- Set the job name 
#BSUB -J Vortex_t8Re100_20x4_dt5
### -- ask for number of cores (default: 1) -- 
#BSUB -n 12 
### -- specify that the cores must be on the same host -- 
#BSUB -R "span[hosts=1]"
### -- specify that we need 2GB of memory per core/slot -- 
#BSUB -R "rusage[mem=2GB]"
### -- specify that we want the job to get killed if it exceeds 20 GB per core/slot -- 
##BSUB -M 20GB
### -- set walltime limit: hh:mm -- 
#BSUB -W 23:00 
#BSUB -N                # send notification at completion
#BSUB -u jesgan@dtu.dk
### -- set the email address -- 
### -- Specify the output and error file. %J is the job-id -- 
### -- -o and -e mean append, -oo and -eo mean overwrite -- 
#BSUB -o Output_%J_%I.out 
#BSUB -e Error_%J_%I.err 

echo "Number of procs per thread"
echo $LSB_DJOB_NUMPROC

### Value to import 
#SOMEVALUE=2.5

######### program to run
echo $PATH
pwd
echo "Start time:"
date

# NOW EXPORT THE PARAMETER YOU CAN READ THROUGH getenv()
#export VALUE=$SOMEVALUE; echo $VALUE

#sleep "$((LSB_JOBINDEX * 5))"
#matlab -nodisplay -r 'matlab_test_script;exit' #FRA FEM VIB
matlab -batch main_vortex_shedding;exit #Fra HPC center

echo "End time:"
date 
