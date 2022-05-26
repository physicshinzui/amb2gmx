#!/bin/bash
set -ue

GMX=gmx 

nwarnings=7
nruns=1
for id in $(seq 1 $nruns); do
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

    cd ..
done
