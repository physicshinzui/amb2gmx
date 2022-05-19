#!/bin/bash 
src=../../src

python $src/change_water_atmname.py system.top

echo "~~~Output test~~~"
diff modified.top system.top > out 
diff diff.log out
