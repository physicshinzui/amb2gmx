#!/bin/bash
#PJM -L "rscunit=ito-a"         # System name: ito-a (CPU nodes) / ito-b (GPU node)
#PJM -L "rscgrp=ito-a-oc190166" # Gropu name: ito-a-oc190166 / ito-b-oc190166
#PJM -L "vnode=1"               # The maximum number of nodes to be used
#PJM -L "vnode-core=36"         # The maximum number of cores to be used
#PJM -L "elapse=96:00:00"          # The maximum time for computation
module load gromacs/2020.6-cpu
set -u

GMX=gmx_mpi #/home/app/gromacs/2020.6/cpu/bin/gmx_mpi

FORCES="500"
Nruns=10

for FORCE in $FORCES; do 
    mkdir -p f$FORCE
    cd f$FORCE
    for id in $(seq 1 $Nruns); do
        mkdir -p $id
        cd $id
        
        # Comment: if outputs does not exist, execute the run.
        if [ ! -e "nvt_eq_${id}.gro" ]; then 
            echo "NVT equilibration runs are running..."
            cat ../../templates/nvt_eq.mdp | sed -e "s!#{RAND}!${RANDOM}!g" > nvt_eq_${id}.mdp
            ${GMX} grompp -f nvt_eq_${id}.mdp -c ../../em.gro -r ../../em.gro -p ../../system.top -o nvt_eq_${id}.tpr
            ${GMX} mdrun -deffnm nvt_eq_${id} 
        fi

        if [ ! -e "npt_eq_${id}.gro" ]; then 
            echo "NPT equilibration runs are running..."
            cp ../../templates/npt_eq.mdp npt_eq_${id}.mdp
            ${GMX} grompp -f npt_eq_${id}.mdp -c nvt_eq_${id}.gro -r nvt_eq_${id}.gro -p ../../system.top -o npt_eq_${id}.tpr
            ${GMX} mdrun -deffnm npt_eq_${id}
        fi

        if [ ! -e "smd_${id}.xtc" ]; then
            echo "Steered MD are runnig..."
            sed -e "s/FORCE/${FORCE}/g" ../../templates/smd_const_f.mdp > smd_${id}.mdp
            ${GMX} grompp -f smd_${id}.mdp -c npt_eq_${id}.gro -t npt_eq_${id}.cpt -r npt_eq_${id}.gro -p ../../system.top -n ../../index.ndx -o smd_${id}.tpr
            ${GMX} mdrun -deffnm smd_${id}
        fi

        cd ..
    done
    cd ..
done  
