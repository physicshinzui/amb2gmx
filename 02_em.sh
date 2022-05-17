#!/bin/bash
set -eu

# Energy Minimisation
GMX=gmx
$GMX grompp -f templates/em.mdp \
            -c mol_solv_ions.gro \
            -o em \
            -p system.top \
            -r mol_solv_ions.gro \
            -maxwarn 1

$GMX mdrun -v -deffnm em # -ntmpi 1 -ntomp 40 #>& em.job
