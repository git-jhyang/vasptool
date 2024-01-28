#!/usr/bin/env python3

import lp34vasp as lp
import numpy as np
import sys, os

def call_help(str1='', str2=''):
   if (str1 != ''): lp.call_error(str1, str2, isexit=False)
   print(' usage   : '+sys.argv[0]+' [DOSCAR]')
   print(' default : [DOSCAR] = ./DOSCAR')
   print('            - file should include \'DOSCAR\'')
   print('            - directory should include \'OUTCAR\'')
   print(' options : [-h] : print this help')
   exit(1)

file_doscar = './DOSCAR'

if len(sys.argv) > 1:
   for arg in sys.argv[1:]:
      if arg[0] == '-':
         for str in arg:
            if str == 'h': call_help()
      elif 'DOSCAR' in arg:
         file_doscar = arg

file_outcar = file_doscar.replace('DOSCAR','OUTCAR')

out = lp.get_outcar(file_outcar)
inp = out.input
dos = lp.get_doscar(file_doscar, inp)

dos.write_dos('dos_total.dat', dos.tdos, dos.header_T)

if dos.lorbit != 0:
   try: 
      os.stat('pdos')
   except:
      os.mkdir('pdos')
   atoms = inp['atoms']
   nspec = inp['nspec']
   pdos_sum = np.zeros([len(atoms),dos.nedos,dos.ncol],dtype=np.float_)
   pdos = dos.get_pdos()
   idx = 0
   for iatom in range(dos.natoms):
      if (iatom >= sum(nspec[0:idx+1])): idx += 1
      pdos_sum[idx] = pdos_sum[idx] + pdos[iatom]
      atom = atoms[idx]
      file_output = 'pdos/dos_%3.3i_%s.dat'%(iatom+1, atom)
      dos.write_dos(file_output, pdos[iatom], dos.header_P)
   for idx in range(len(atoms)):
      file_output = 'pdos/dos_total_%s.dat'%(atoms[idx])
      for idos in range(dos.nedos): pdos_sum[idx][idos][0] = pdos[idx][idos][0]
      dos.write_dos(file_output, pdos_sum[idx], dos.header_P)
