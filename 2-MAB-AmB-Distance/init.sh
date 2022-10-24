#!/bin/bash

# Set up simulation environment
source env.sh

# Clean up from previous/ failed runs

read -p "Clear old simulation files? (y/n)? `echo $'\n> '`" CONT
if [ "$CONT" = "y" ]; then
	echo "Purging old WE data!"
	rm -rf traj_segs seg_logs istates west.h5 
	mkdir   seg_logs traj_segs istates

	mkdir -p bstates/unbound
	mkdir -p job_logs
	rm -f job_logs/*
	cp prep/1.eqn.coor bstates/unbound/seg.coor
	cp prep/1.eqn.dcd  bstates/unbound/seg.dcd
	cp prep/1.eqn.vel  bstates/unbound/seg.vel
	cp prep/1.eqn.xsc  bstates/unbound/seg.xsc

	cp prep/1.psf   namd_config/nacl.psf
	cp prep/1.pdb   namd_config/nacl.pdb
else
	echo "Old data preserved!"
fi


# Set pointer to bstate and tstate
BSTATE_ARGS="--bstate-file bstates/bstates.txt"
TSTATE_ARGS="--tstate bound,4.9"

echo ${BSTATE_ARGS}
# Run w_init
w_init \
  $BSTATE_ARGS \
  $TSTATE_ARGS \
  --segs-per-state 5 \
  --work-manager=threads "$@"
  
