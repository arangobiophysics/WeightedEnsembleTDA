#!/bin/bash
#
# get_pcoord.sh
#
# This script is run when calculating initial progress coordinates for new 
# initial states (istates).  This script is NOT run for calculating the progress
# coordinates of most trajectory segments; that is instead the job of runseg.sh.

# If we are debugging, output a lot of extra information.
if [ -n "$SEG_DEBUG" ] ; then
  set -x
  env | sort
fi

# Make sure we are in the correct directory
cd $WEST_SIM_ROOT
source env.sh
cd $WEST_STRUCT_DATA_REF


# Make a temporary file in which to store output from the python script
DIST=$(mktemp)


# Symlink a file needed for analysis
ln -sfn $WEST_SIM_ROOT/namd_config/nacl.psf structure.psf
ln -sfn $WEST_SIM_ROOT/namd_config/nacl.pdb structure.pdb
ln -sfn $WEST_SIM_ROOT/westpa_scripts/test_persistent_homology.py test_persistent_homology.py 

# Use a custom python script to calculate the distance between the Na+ and Cl-
# ions. This script looks for files named 'nacl.psf' and 'seg.dcd'.
/Scr/hyunpark/anaconda3/envs/deeplearning/bin/python test_persistent_homology.py --data_dir $WEST_SIM_ROOT/westpa_scripts --psf structure.psf --pdb structure.pdb --trajs seg.dcd --selections "name CA" > $DIST

# Pipe the relevant part of the output file (the distance) to $WEST_PCOORD_RETURN
cat $DIST | tail -n 1 > $WEST_PCOORD_RETURN

#debugging
echo "$WEST_PCOORD_RETURN" >> pcoord2.dat
echo "$DIST" >> pcoord.dat
# Remove the temporary file to clean up
rm $DIST

if [ -n "$SEG_DEBUG" ] ; then
  head -v $WEST_PCOORD_RETURN
fi