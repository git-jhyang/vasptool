#!/usr/bin/env python3

import lp34vasp as lp
import sys

if len(sys.argv) != 2: lp.call_error('input error - '+sys.argv[0], '1 arg needed')

lp.assert_file(sys.argv[1], sys.argv[0], isopen=False)

path = ''
path_split = sys.argv[1].split('/')
for str in path_split[:-1]:
   path += str + '/'

path_poscar = sys.argv[1]
pos         = lp.get_poscar(path_poscar)

path_outcar = path+'OUTCAR'

if lp.assert_file(path_outcar,isopen=False):
   out = lp.get_outcar(path_outcar)
   if (not out.is_end): lp.call_error('output error','calculation not finished')
   pos.force = out.get_force()
   pos.des   = 'total energy = %.8f eV'%(out.get_energy()['E_sigma0'])

lp.write_xsf(pos.filename+'.xsf', pos)

