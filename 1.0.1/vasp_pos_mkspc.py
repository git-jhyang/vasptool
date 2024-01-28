#!/usr/bin/env python

import lp24vasp as lp
import numpy as np
import sys

cnt = 0
dim = [1, 1, 1]

def call_help(str1='', str2=''):
   if (str1 != ''): lp.call_error(str1, str2, isexit=False)
   print ' usage   : '+sys.argv[0]+' [option] [filename] [nx] [ny] [nz]'
   print ' default : [filename] = None'
   print '           [option]   = None'
   print '           [n]        = 1'
   print ' options : [-h] : print this help'
   exit(1)

if len(sys.argv) == 1: call_help('input arguments','need more than 1 arg')

for arg in sys.argv[1:]:
   if arg[0] == '-':
      for str in arg:
         if str == 'h': call_help()
   elif not lp.assert_num(arg):
      lp.assert_file(arg, sys.argv[0], isopen = False)
      path = arg
   else:
      dim[cnt] = int(arg)
      cnt += 1
      if cnt == 3: break

pos    = lp.get_poscar(path)
pos_re = lp.get_poscar(path)

for i in range(3): pos_re.lvec[i][:] = pos.lvec[i][:]*dim[i]
for i in range(len(pos.nspec)): pos_re.nspec[i] = pos.nspec[i]*dim[0]*dim[1]*dim[2]
pos_re.natoms = pos.natoms*dim[0]*dim[1]*dim[2]
pos_re.coord = np.zeros([pos_re.natoms,3], dtype = np.float_)
if pos.isselc: pos_re.selc = [['T','T','T']]*pos_re.natoms
if pos.iscart: lp.cart_to_dir(pos.lvec, pos.coord)

cnt = 0
coo_tmp= np.zeros(3,dtype=np.float_)
coo_dis = np.zeros(3,dtype=np.float_)
for i in range(pos.natoms):
   for j in range(3): coo_tmp[j] = np.float_(pos.coord[i][j]/dim[j])
   for ix in range(dim[0]):
      coo_dis[0] = np.float(ix)/dim[0]
      for iy in range(dim[1]):
         coo_dis[1] = np.float(iy)/dim[1]
         for iz in range(dim[2]):
            coo_dis[2] = np.float(iz)/dim[2]
            pos_re.coord[cnt][:] = coo_tmp + coo_dis
            if pos.isselc: pos_re.selc[cnt] = pos.selc[i]
            cnt += 1

lp.write_poscar('POSCAR_sc_%ix%ix%i'%(dim[0],dim[1],dim[2]), pos_re)
