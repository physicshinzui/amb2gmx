#!/bin/bash
#PJM -L "rscunit=ito-a"         # System name: ito-a (CPU nodes) / ito-b (GPU node)
#PJM -L "rscgrp=ito-a-oc190166" # Gropu name: ito-a-oc190166 / ito-b-oc190166
#PJM -L "vnode=1"               # The maximum number of nodes to be used
#PJM -L "vnode-core=36"         # The maximum number of cores to be used
#PJM -L "elapse=96:00:00"          # The maximum time for computation
module load gromacs/2020.6-cpu
set -e

GMX=gmx_mpi

nwarnings=7
id=3

mkdir -p $id
cd $id

# Comment: if outputs does not exist, execute the run.
if [ ! -e "nvt_eq_${id}.gro" ]; then 
    echo "NVT equilibration runs are running..."
    cat ../templates/nvt_eq.mdp | sed -e "s!#{RAND}!${RANDOM}!g" > nvt_eq_${id}.mdp
    ${GMX} grompp -f nvt_eq_${id}.mdp -c ../em.gro -r ../em.gro -p ../system.top -o nvt_eq_${id}.tpr -maxwarn $nwarnings
    ${GMX} mdrun -deffnm nvt_eq_${id} 
fi

if [ ! -e "npt_eq_${id}.gro" ]; then 
    echo "NPT equilibration runs are running..."
    cp ../templates/npt_eq.mdp npt_eq_${id}.mdp
    ${GMX} grompp -f npt_eq_${id}.mdp -c nvt_eq_${id}.gro -r nvt_eq_${id}.gro -p ../system.top -o npt_eq_${id}.tpr -maxwarn $nwarnings
    ${GMX} mdrun -deffnm npt_eq_${id}
fi

if [ ! -e "npt_prod_${id}.gro" ]; then 
    echo "NPT production runs are running..."
    cp ../templates/npt_prod.mdp npt_prod_${id}.mdp
    ${GMX} grompp -f npt_prod_${id}.mdp -c npt_eq_${id}.gro -p ../system.top -o npt_prod_${id}.tpr -maxwarn $nwarnings
    ${GMX} mdrun -deffnm npt_prod_${id}
fi

cd ..
