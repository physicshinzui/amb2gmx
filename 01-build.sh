#!/bin/bash
set -eu

PDB=$1
FF="ff14SB"

# === (0) Preprocessing: 
gmx editconf -center 0 0 0 -f ${PDB} -o $(basename ${PDB})

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
acpype -p system.prmtop -x system.inpcrd -b system
# NOTE: This generates `system.amb2gmx` directory

# === (4) Keep Protein and Ions only
cp -p system.amb2gmx/system_GMX.gro system.gro
cp -p system.amb2gmx/system_GMX.top system.top

echo '"Protein" | "Ion"' > inp
echo quit >> inp
gmx make_ndx -f system.gro -o index_ProeinIon.ndx < inp
echo Protein_Ion | gmx trjconv -f system.gro -s system.gro -o system.gro -n index_ProeinIon.ndx

cat system.top   | \
        sed '$d' | \
        sed  -e "s/WAT/SOL/g" -e "" > tmp

mv tmp system.top
rm inp index_ProeinIon.ndx

# === (5) Generate position restraint for main chain
posres_itp=posre_mainchain.itp
echo MainChain | gmx genrestr -f system.gro -o ${posres_itp}

# === (6) Insert ifdef for position restraints.
iloc=$(nl -b a system.top | grep "moleculetype" | awk 'NR==2 {print $1}' | bc)
gsed -i "${iloc} i #ifdef POSRES\n#include \"${posres_itp}\"\n#endif\n" system.top 
# NOTE: 
# This is GNU sed installed by Homebrew on mac. If you are in a linux, replace gsed with sed

# === (7) Build a system by gmx.
echo Protein | gmx editconf -f system.gro -o box.gro -princ -bt dodecahedron -d 1.0

## === (7.1) insert molecule for cosolvent MD: 
#gmx editconf -f pdb4amber.pdb -o pdb4amber.pdb -d 1.5 -bt cubic
#nmol=10
#gmx insert-molecules -f pdb4amber.pdb -ci cosolv.pdb -o pdb4amber.pdb -nmol $nmol -scale 3

gmx solvate -cp box.gro -cs spc216.gro -o mol_solv.gro -p system.top

python src/change_water_atmname.py system.top
cp -p system.top system.top.bak
mv modified.top system.top

gmx grompp -f templates/ions.mdp -c mol_solv.gro -p system.top -o ions.tpr 
echo SOL | gmx genion -s ions.tpr -o mol_solv_ions.gro -p system.top -pname NA+ -nname CL- -neutral

bash trash.sh
