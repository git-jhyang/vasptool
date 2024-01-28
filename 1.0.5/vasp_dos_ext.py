#!/usr/bin/env python3

import pyutil.lp34vasp as lp
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
    if not os.path.isdir('pdos'):
        os.mkdir('pdos')
    Atoms = inp['atoms']
    nSpec = inp['nspec']
    nType = len(Atoms)
    for iType,nAtom in enumerate(nSpec):
        Atom     = Atoms[iType]
        pdos_sum = np.zeros([dos.nedos,dos.ncol],dtype=np.float_)
        for iAtom in range(nAtom):
            idx         = iAtom + sum(nSpec[0:iType])
            pdos_sum   += dos.pdos[idx]
            file_output = 'pdos/dos_%3.3i_%s.dat'%(idx+1, Atom)
            dos.write_dos(file_output, dos.pdos[idx], dos.header_P)
        file_output = 'pdos/dos_total_%s.dat'%(Atom)
        pdos_sum[:,0] = dos.pdos[iType,:,0]
        dos.write_dos(file_output, pdos_sum, dos.header_P)
