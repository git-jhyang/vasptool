#!/usr/bin/env python3

import pyutil.lp34vasp as lp
import numpy as np
import sys

cnt = 0
dim = [1, 1, 1]

def call_help(str1='', str2=''):
    if (str1 != ''): lp.call_error(str1, str2, isexit=False)
    print(' usage   : '+sys.argv[0]+' [option] [filename] [nx] [ny] [nz]')
    print(' default : [filename] = None')
    print('           [option]   = None')
    print('           [n]        = 1')
    print(' options : [-h] : print this help')
    exit(1)

if len(sys.argv) == 1: 
    call_help('input arguments','need more than 1 arg')

for arg in sys.argv[1:]:
    if arg[0] == '-':
        for str in arg:
            if str == 'h': call_help()
    if lp.assert_num(arg):
        dim[cnt] = int(arg)
        cnt += 1
        if cnt == 3: 
            break
    else:
        isfile = lp.assert_file(arg, sys.argv[0], callerr=False)
        if isfile: path = arg

pos    = lp.get_poscar(path)
pos_re = lp.get_poscar(path)

for i in range(3): 
    pos_re.lvec[i] = pos.lvec[i]*dim[i]
for iType,nAtom in enumerate(pos.nSpec):
    pos_re.nSpec[iType] = nAtom*dim[0]*dim[1]*dim[2]
pos_re.nAtom = pos.nAtom*dim[0]*dim[1]*dim[2]
pos_re.coord = np.zeros([pos_re.nAtom,3], dtype = np.float64)
pos_re.vector = None
if pos.isselc: 
    pos_re.selc = [['T','T','T']]*pos_re.nAtom
if pos.iscart: 
    lp.cart_to_dir(pos.lvec, pos.coord)

cnt = 0
coo_dis = np.zeros(3,dtype=np.float64)
for i in range(pos.nAtom):
    coo_tmp = np.float64(pos.coord[i]/dim)
    for ix in range(dim[0]):
        coo_dis[0] = np.float64(ix/dim[0])
        for iy in range(dim[1]):
            coo_dis[1] = np.float64(iy/dim[1])
            for iz in range(dim[2]):
                coo_dis[2] = np.float64(iz/dim[2])
                pos_re.coord[cnt] = coo_tmp + coo_dis
                if pos.isselc: 
                    pos_re.selc[cnt] = pos.selc[i]
                cnt += 1

pos_re.f_name = 'POSCAR_sc_%ix%ix%i'%(dim[0],dim[1],dim[2])
lp.write_poscar(pos_re)
