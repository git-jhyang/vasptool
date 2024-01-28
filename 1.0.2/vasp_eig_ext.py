#!/usr/bin/env python3

import lp34vasp as lp
import numpy as np
import sys

def call_help(str1='', str2=''):
   if (str1 != ''): lp.call_error(str1, str2, isexit=False)
   print(' usage   : '+sys.argv[0]+' [option] [EIGENVAL] [OUTCAR]')
   print(' default : [EIGENVAL] = ./EIGENVAL')
   print('           [OUTCAR]   = ../OUTCAR')
   print('             files should include \'EIGENVAL\' and \'OUTCAR\', respectively')
   print(' options : [-h] : print this help')
#   print('         : [-v] : vertically arranged outfile, default')
#   print('         : [-h] : horizentally arranged outfile')
   exit(1)

out_v = True
out_h = True
file_outcar = '../OUTCAR'
file_eigenv = './EIGENVAL'

if len(sys.argv) > 1:
   for arg in sys.argv[1:]:
      if arg[0] == '-':
         for str in arg:
            if str == 'h': call_help()
#            if str == 'v': out_v = True
#            if str == 'h': out_h = True
      elif 'OUTCAR' in arg:
         file_outcar = arg
      elif 'EIGENVAL' in arg:
         file_eigenv = arg

out = lp.get_outcar(file_outcar)
eig = lp.get_eigenval(file_eigenv)

if not lp.assert_dat(int(out.input['nelect']), eig.nelect):
   lp.call_error('input error - NELECT','two files from different calculation')
if not lp.assert_dat(out.input['ispin'], eig.ispin): 
   lp.call_error('input error - ISPIN','two files from different calculation')

lvec_rec = out.get_lattice()[1]
ef       = out.get_energy()['E_fermi']

if out.input['ncollin'] == 'T':
   icbm = eig.nelect 
else:
   icbm = int(eig.nelect/2)

ivbm = icbm - 1
evbm = eig.eigens[0][ivbm][0] - ef
kvbm = 0
svbm = 0

ecbm = eig.eigens[0][icbm][0] - ef
kcbm = 0
scbm = 0

ismetal = False
lenkpt  = np.zeros([eig.nkpts],dtype=np.float_)
vec     = np.zeros([3],dtype=np.float_)

for iKpt in range(eig.nkpts):
   eig.kpts[iKpt]   = np.matmul(eig.kpts[iKpt],lvec_rec)
   eig.eigens[iKpt] = eig.eigens[iKpt] - ef
   if iKpt != 0:
      vec = eig.kpts[iKpt] - eig.kpts[iKpt-1]
      lenkpt[iKpt] = lenkpt[iKpt-1] + np.sqrt(np.dot(vec,vec))
   for iSpin in range(eig.ispin):
      if ismetal: break
      if eig.eigens[iKpt][ivbm][iSpin] > evbm:
         evbm = eig.eigens[iKpt][ivbm][iSpin]
         kvbm = iKpt
         svbm = iSpin
         if evbm > 0: ismetal = True
      if eig.eigens[iKpt][icbm][iSpin] < ecbm: 
         ecbm = eig.eigens[iKpt][icbm][iSpin]
         kcbm = iKpt
         scbm = iSpin
         if ecbm < 0: ismetal = True

if not ismetal: egap = ecbm - evbm

if out_h:
   f = open('eigenval.h.dat','w')
   f.write('# Type: H\n')
   if ismetal:
      f.write('# %10s %10s %10s\n'%('ISPIN','NBANDS','NKPT'))
      f.write('# %10i %10i %10i\n'%(eig.ispin, eig.nbands, eig.nkpts))
   else:
      f.write('# %10s %10s %10s %10s %10s'%('ISPIN','NBANDS','NKPT','EF','E_GAP'))
      f.write(' %10s %10s %10s %10s %10s'%('VBM_N','VBM_K','VBM_S','VBM_E','VBM_E-EF'))
      f.write(' %10s %10s %10s %10s %10s\n'%('CBM_N','CBM_K','CBM_S','CBM_E','CBM_E-EF'))
      f.write('# %10i %10i %10i %10.6f %10.6f'%(eig.ispin, eig.nbands, eig.nkpts, ef, egap))
      f.write(' %10i %10i %10i %10.6f %10.6f'%(ivbm+1, kvbm+1, svbm+1, evbm+ef, evbm))
      f.write(' %10i %10i %10i %10.6f %10.6f\n'%(icbm+1, kcbm+1, scbm+1, ecbm+ef, ecbm))

   f.write('# %10s %10s %10s %10s '%('k_len','kx','ky','kz'))
   f.write('%12s %12s %12s\n'%('Band_1','Band_2','...'))

   for iSpin in range(eig.ispin):
      f.write('# spin component : %i\n'%(iSpin+1))
      for iKpt in range(eig.nkpts):
         f.write('  %10.6f'%(lenkpt[iKpt]))
         for i in range(3): f.write(' %10.6f'%(eig.kpts[iKpt][i]))
         for iBand in range(eig.nbands):
            f.write(' %12.6f'%(eig.eigens[iKpt][iBand][iSpin]))
         f.write('\n')
      f.write('#\n')
   f.close()

if out_v:
   f = open('eigenval.v.dat','w')
   f.write('# Type: V\n')
   if ismetal:
      f.write('# %10s %10s %10s\n'%('ISPIN','NBANDS','NKPT'))
      f.write('# %10i %10i %10i\n'%(eig.ispin, eig.nbands, eig.nkpts))
   else:
      f.write('# %10s %10s %10s %10s %10s'%('ISPIN','NBANDS','NKPT','EF','E_GAP'))
      f.write(' %10s %10s %10s %10s %10s'%('VBM_N','VBM_K','VBM_S','VBM_E','VBM_E-EF'))
      f.write(' %10s %10s %10s %10s %10s\n'%('CBM_N','CBM_K','CBM_S','CBM_E','CBM_E-EF'))
      f.write('# %10i %10i %10i %10.6f %10.6f'%(eig.ispin, eig.nbands, eig.nkpts, ef, egap))
      f.write(' %10i %10i %10i %10.6f %10.6f'%(ivbm+1, kvbm+1, svbm+1, evbm+ef, evbm))
      f.write(' %10i %10i %10i %10.6f %10.6f\n'%(icbm+1, kcbm+1, scbm+1, ecbm+ef, ecbm))

   for iBand in range(eig.nbands):
      f.write('# Eigenstate index : %6i\n'%(iBand+1))
      f.write('# %10s %10s %10s %10s %12s '%('k_len','kx','ky','kz','Eigen_1'))
      if eig.ispin == 2:
         f.write('%12s'%('Eigen_2'))
      f.write('\n')
      for iKpt in range(eig.nkpts):
         f.write('  %10.6f'%(lenkpt[iKpt]))
         for i in range(3): f.write(' %10.6f'%(eig.kpts[iKpt][i]))
         for iSpin in range(eig.ispin):
            f.write(' %12.6f'%(eig.eigens[iKpt][iBand][iSpin]))
         f.write('\n')
      f.write('#\n')
   f.close()

