#!/bin/bash
set -eu

PDB=$1
FF="ff14SB"
GMX=gmx

# === (0) Preprocessing: 
${GMX} editconf -center 0 0 0 -f ${PDB} -o $(basename ${PDB})

# === (1) pdb4amber: 
pdb4amber -i $(basename $PDB) -o pdb4amber.pdb

mkdir -p log.pdb4amber ; mv pdb4amber_* log.pdb4amber

# === (2) tleap: topology file generation
cat templates/tleap.in | \
    sed -e "s!#{INPUT}!${PDB}!g" \
        -e "s!#{FF}!${FF}!g" > tleap.in

# NOTE: tleap generates system.prmtop and system.inpcrd
tleap -f tleap.in 
mkdir -p log.tleap ; mv tleap* log.tleap

# === (3) acpype: Convertion of Amber to Gromacs topology file

## ===(3.1) For protein
acpype -p system.prmtop -x system.inpcrd -b system
# NOTE: This generates `system.amb2gmx` directory

## ==(3.2) For ligand or cosolvent
#ligdir=lig
#acpype -p $ligdir/ligand.parm7 -x $ligdir/ligand.rst7 -b ligand

# === (4) Keep Protein and Ions only
cp -p system.amb2gmx/system_GMX.gro system.gro
cp -p system.amb2gmx/system_GMX.top system.top

echo '"Protein" | "Ion"' > inp
echo quit >> inp
${GMX} make_ndx -f system.gro -o index_ProeinIon.ndx < inp
echo Protein_Ion | ${GMX} trjconv -f system.gro -s system.gro -o system.gro -n index_ProeinIon.ndx

cat system.top   | \
        sed '$d' | \
        sed  -e "s/WAT/SOL/g" -e "" > tmp

mv tmp system.top
rm inp index_ProeinIon.ndx

# === (5) Generate position restraint for main chain
posres_itp=posre_mainchain.itp
echo MainChain | ${GMX} genrestr -f system.gro -o ${posres_itp}

# === (6) Insert ifdef for position restraints.
iloc=$(nl -b a system.top | grep "moleculetype" | awk 'NR==2 {print $1}' | bc)
gsed -i "${iloc} i #ifdef POSRES\n#include \"${posres_itp}\"\n#endif\n" system.top 
# NOTE: 
# This is GNU sed installed by Homebrew on mac. If you are in a linux, replace gsed with sed

# === (7) Build a system by gmx.
echo Protein | ${GMX} editconf -f system.gro -o box.gro -princ -bt cubic -d 1.5

## === (7.1) insert molecule for cosolvent MD: 
lig_crd=lig/em.gro
nmol=10
${GMX} insert-molecules -f box.gro -ci $lig_crd -o box_inserted.gro -nmol $nmol -scale 3

## === (7.2) Solvate
${GMX} solvate -cp box_inserted.gro -cs spc216.gro -o mol_solv.gro -p system.top

## === (7.3) Merge topologies: protein and ligand
# TODO: This part is hard coded and mannually done. I want to automate here.  
lig_top=lig/hs_GMX.top
vim -s src/get_atomtype.vim $lig_top # -> ligand.itp generated
vim -s src/top2itp.vim      $lig_top # -> atomtype.txt generated
echo "=============================================================================================="
echo "| Manually Insertion of ligand topology files (ligand.itp and atomtype.txt)... Go on? [Enter]|"
echo "=============================================================================================="
read 
vim system.top

## == (7.3) Rename gromacs default atom names of water molecule to those of Amber's
python src/change_water_atmname.py system.top
cp -p system.top system.top.bak
mv modified.top system.top

## ==(7.4) Add ions
read -p "How many warnings are allowd?" nwarns
${GMX} grompp -f templates/ions.mdp -c mol_solv.gro -p system.top -o ions.tpr -maxwarn $nwarns
echo SOL | ${GMX} genion -s ions.tpr -o mol_solv_ions.gro -p system.top -pname NA+ -nname CL- -neutral

bash trash.sh

cat <<EOF
=================================
Note: 
You may see a few warnings, this is because redundant atom types for the ligand given.
There are two ways: 
    i. set -maxwarn n 
    ii. remove the redandant atom types from system.top.
==================================
EOF

