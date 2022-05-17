#!/usr/bin/env python
import parmed as pmd

from parmed import gromacs
gmx_top = gromacs.GromacsTopologyFile('system.amb2gmx/system_GMX.top')
gmx_top.write('00-inputs/unminimized-gas.amb2gmx/unminimized-gas_GMX.top',[[0,1]])

## for amber
#parm1 = pmd.load_file('system.prmtop', 'system.inpcrd')
#parm2 = pmd.load_file('00-inputs/unminimized-gas.parm7', '00-inputs/unminimized-gas.rst7')
#joined = parm1 + parm2
#joined.save('pro_cosolv.parm7')
#joined.save('pro_cosolv.rst7')
