#!/bin/bash

# necessary but nonetheless demonstrates good practice.
export NAMD="/Projects/dhardy/namd_builds/NAMD_3.0alpha11_Linux-x86_64-multicore-CUDA/namd3"
#export NAMD=$(which namd2)
############################## Python and WESTPA ###############################
# Next inform WESTPA what python it should use.  
#export WEST_PYTHON=$(which python2.7)
#export WEST_PYTHON=$(/Projects/arango/anaconda3/envs/westpa-2.0/bin/python333)
#export WEST_PYTHON="/Projects/arango/anaconda3/envs/py27/bin/python2.7"
# Check to make sure that the environment variable WEST_ROOT is set. 
# Here, the code '[ -z "$WEST_ROOT"]' will return TRUE if WEST_ROOT is not set,
# causing us to enter the if-statement, print an error message, and exit.


# Set up environment for westpa
#export WEST_PYTHON=$(which python)
# Actviate a conda environment containing westpa, openmm, mdtraj;
# You may need to create this first (see install instructions)
#source /Projects/arango/anaconda3/envs/westpa-2020.02/westpa-2020.03/westpa.sh 
#source activate westpa-2020.02 
#source activate westpa-openmm-mdtraj


#if [ -z "$WEST_ROOT" ]; then
#  echo "The environment variable WEST_ROOT is not set."
#  echo "Try running 'source westpa.sh' from the WESTPA installation directory"
#  exit 1
#fi

#temporary fix
#export WEST_SIM_ROOT=/Scr/arango/membrane_AA/WESTPA/AA_memb_WESTPA
export WEST_SIM_ROOT=$PWD
#if [ -z "$WEST_SIM_ROOT" ]; then
#  export WEST_SIM_ROOT="$PWD"
#fi

export SIM_NAME=$(basename $WEST_SIM_ROOT)


