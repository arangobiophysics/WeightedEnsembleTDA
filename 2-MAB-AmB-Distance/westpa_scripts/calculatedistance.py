import persim
import ripser
import MDAnalysis as mda
import argparse
from typing import *
import functools
import itertools 

parent_universe = mda.Universe('structure.psf', 'seg.dcd')
#--------------------
#parent_universe.trajectory[-1]
#REFERENCE FINAL STRUCTURE
reference_universe = mda.Universe('/Scr/arango/Sobolev-Hyun/test_simulations/adk_equilibrium/adk4AKE.psf', '/Scr/arango/Sobolev-Hyun/test_simulations/adk_equilibrium/ADK_last.pdb')
parent_sel = reference_universe.select_atoms('name CA')
#Doing up to H1
parent_homology = ripser.ripser(parent_sel.atoms.positions, maxdim=1)['dgms']

#--------------------
#Comparing the dcd files to the desired snapshot reference
x=[]
for ts in parent_universe.trajectory:
    print(ts.frame)
    Pa1 =  parent_universe.select_atoms("name CA")
    R1 = ripser.ripser(Pa1.atoms.positions, maxdim=1)['dgms']
    perwas = persim.wasserstein(parent_homology[1], R1[1])
    print("{:.03f}".format(perwas))
