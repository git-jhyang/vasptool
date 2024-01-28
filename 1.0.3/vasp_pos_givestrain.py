#!/usr/bin/env python3

import lp34vasp as lp
import numpy as np
import sys

def call_help(str1='',str2=''):
   if (str1 != ''): lp.call_error(str1, str2, isexit=False)
   print(' usage   : '+sys.argv[0]+' [filename] [axis] [str_max] -n [images]')
   print(' default : [filename] = None')
   print('           [axis]     = -z   (strain direction)')
   print('           [str_max]  = 0.10 (maximum strain)')
   print('           [images]   = 4    (number of images)')
   print()
   print(' options : [-h]           : print this help')
   print('           [-i] [str_max] : isotropic strain')
   print('           [-x] [str_max] : strain along x direction')
   print('           [-y] [str_max] : strain along y direction')
   print('           [-z] [str_max] : strain along z direction')
   print('           [-n] [images]  : number of images')
   exit(1)

if (len(sys.argv) == 1): call_help('input error','need more than 1 arguments')
if (sys.argv[1] == '-h'): call_help()

lp.assert_file(sys.argv[1], sys.argv[0], isopen = False)
path = sys.argv[1]

axis   = [False, False, False, False]
strain = np.zeros(4, dtype=np.float32)
images = 4

for i, arg in enumerate(sys.argv):
   if (i < 2): continue
   if arg[0] != '-': continue
   if lp.assert_num(arg): continue
   for str in arg[1:]:
      if str == 'h': call_help()
      elif str == 'i':
         if (len(sys.argv) < i): call_help('input arguments','strain is missing')
         if (not lp.assert_num(sys.argv[i+1])): call_help('input arguments -i','not a number '+sys.argv[i+1])
         axis[0] = True
         strain[0] = np.float(sys.argv[i+1])
      elif str == 'x':
         if (len(sys.argv) < i): call_help('input arguments','strain is missing')
         if (not lp.assert_num(sys.argv[i+1])): call_help('input arguments -x','not a number '+sys.argv[i+1])
         axis[1] = True
         strain[1] = np.float(sys.argv[i+1])
      elif str == 'y':
         if (len(sys.argv) < i): call_help('input arguments','strain is missing')
         if (not lp.assert_num(sys.argv[i+1])): call_help('input arguments -y','not a number '+sys.argv[i+1])
         axis[2] = True
         strain[2] = np.float(sys.argv[i+1])
      elif str == 'z':
         if (len(sys.argv) < i): call_help('input arguments','strain is missing')
         if (not lp.assert_num(sys.argv[i+1])): call_help('input arguments -z','not a number '+sys.argv[i+1])
         axis[3] = True
         strain[3] = np.float(sys.argv[i+1])
      elif str == 'n':
         if (len(sys.argv) < i): call_help('input arguments','variable is missing')
         if (not lp.assert_num(sys.argv[i+1])): call_help('input arguments -n','not a number '+sys.argv[i+1])
         images = int(sys.argv[i+1])
      else:
         call_help('input arguments','not supported argument -'+str)

if (not (axis[0] or axis[1] or axis[2])): 
   axis[3] = True
   strain[3] = 0.1

if (axis[0]):
   for i in range(1,4):
      axis[i] = True
      strain[i] = strain[0]
       
pos    = lp.get_poscar(path)
pos_re = lp.get_poscar(path)

if pos.iscart: lp.cart_to_dir(pos.lvec, pos_re.coord)
if pos.iscart: pos_re.iscart = False

xyz = ['x','y','z']
for img in range(1,images+1):
   for i in range(3):
      if axis[i+1]: 
         ratio = 1e0 + strain[i+1]*img/images
         for j in range(3): pos_re.lvec[j][i] = pos.lvec[j][i]*ratio
   if axis[0]:
      str = 'strain_iso=%.3f'%(1e0+strain[0]*img/images)
   else:
      str = 'strain'
      for i in range(3):
         if axis[i+1]: str += '_%s_%.3f'%(xyz[i], 1e0+strain[i+1]*img/images)
   pos_re.des = str+'\n'
   lp.write_poscar('POSCAR_%s'%(str[7:]), pos_re)

