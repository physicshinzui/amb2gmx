#!/bin/bash 

read -p "ok?"
rm \#* *.gro pdb4amber* *itp tleap.in ions.tpr leap.log mdout.mdp em.*
rm -rf system.*
rm -rf pdb4amber.log
rm selection.txt 
rm index.ndx
rm -rf log.*
