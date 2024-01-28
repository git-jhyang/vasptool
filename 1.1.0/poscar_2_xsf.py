#!/usr/bin/env python3

import pyutil.lp34vasp as lp
import sys

if len(sys.argv) != 2: 
    lp.call_error('input error - '+sys.argv[0], '1 arg needed')

lp.assert_file(sys.argv[1], sys.argv[0])

path_poscar = sys.argv[1]
path        = lp.get_path(sys.argv[1])
path_outcar = path+'OUTCAR'
pos         = lp.get_poscar(path_poscar)
pos.vector  = None

if lp.assert_file(path_outcar, sys.argv[0]):
    out = lp.get_outcar(path_outcar)
    if (not out.is_end): 
        lp.call_error('output error','calculation not finished')
    pos.vector = out.get_force()
    pos.des    = 'total energy = %.8f eV'%(out.energy['E_sigma0'])

lp.write_xsf(pos)

