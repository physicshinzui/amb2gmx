import parmed as pmd
amber = pmd.load_file('pro_cosolv.parm7', 'pro_cosolv.rst7')

# Save a GROMACS topology and GRO file
amber.save('pro_cosolv.top')
amber.save('pro_cosolv.gro')

