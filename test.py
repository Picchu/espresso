# example of dict fields used to populate an instance of QEInput in the database
# this format allows the user to specifies all the keywords of a quantum espresso
# input file. The information on the structure will be taken by the appropriate
# Structure object linked to the Calc as input
prefix = 'Li_C19_coarse_50Ry'
work_dir = "/Users/Picchu/Documents/Bosch/Project/qe_scripts/wf_test"
pseudo_dir = "/Users/Picchu/Documents/PhD/QEspresso/espresso-5.1.2/" +\
"pslibrary.1.0.0/pbesol/PSEUDOPOTENTIALS_US_HIGH"

control = {'title': prefix, 'calculation': 'vc-relax', 'restart_mode': 'from_scratch',
'outdir': work_dir, 'wfcdir': work_dir, 'pseudo_dir': pseudo_dir, 'prefix': prefix, 'verbosity': 'default',
'wf_collect': '.true.', 'etot_conv_thr': '1.0d-4', 'forc_conv_thr': '1.0d-3', 'tstress': '.true.', 'tprnfor': '.true.'}

system = {'ibrav': '0', 'occupations': 'smearing', 'smearing': 'mp', 'degauss': '0.03',
'ecutwfc': '50', 'ecutrho': '400'}

electrons = {'electron_maxstep': '100', 'conv_thr': '1.0d-8', 'mixing_mode': 'plain', 'mixing_beta': '0.2',
'mixing_ndim': '8', 'diagonalization': 'david', 'diago_david_ndim': '4'}

ions = {'ion_dynamics': 'bfgs'}

cell = {'cell_dynamics': 'bfgs'}
