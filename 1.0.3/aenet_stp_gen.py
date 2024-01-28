#!/usr/bin/env python3

import sys, os
import copy as cp
import lp34vasp as lp

def get_int(string, ref=None):
   while True:
      try:
         get = int(input(string))
      except ValueError:
         continue
      else:
         if ref is None: break
         if get-1 in range(ref): break
         continue
   return get

def get_floats(string):
   while True:
      gets = input(string).split()
      try:
         for i in gets: float(i)
      except ValueError:
         print(i, gets)
         continue
      else:
         break
   out = []
   for i in gets: out.append(float(i))
   return out

def conv_floats(string):
   buf = string.split()
   for idx,val in enumerate(buf):
      buf[idx] = float(val)
   return buf

readfile = False
for idx, arg in enumerate(sys.argv):
   if arg == '-r':
      if len(sys.argv) < idx+2: call_help('input error','Need more arguments')
      inp = lp.assert_file(sys.argv[idx+1], sys.argv[0])
      readfile = True

if readfile:
   NameSpec = inp.readline().split()
else:
   print()
   NameSpec = input(' Elements      : ').split()

nSpec    = len(NameSpec)
two_body = []

for i, Atom1 in enumerate(NameSpec):
   tmp_1 = [Atom1]
   for j in range(i,nSpec):
      Atom2 = NameSpec[j]
      tmp_2 = cp.copy(tmp_1)
      tmp_2.append(Atom2)
      two_body.append(tmp_2)

if readfile:
   func_bond   = int(inp.readline().split()[1])
   val_bond_Rc = conv_floats(inp.readline())
   nFuncB = len(val_bond_Rc)
   if func_bond == 2:
      val_bond_eta = conv_floats(inp.readline())
      val_bond_Rs  = conv_floats(inp.readline())
      nFuncB = nFuncB*len(val_bond_eta)*len(val_bond_Rs)
   if func_bond == 3:
      val_bond_kpa = conv_floats(inp.readline()) 
      nFuncB = nFuncB*len(val_bond_kpa)
   func_angle  = int(inp.readline().split()[1])
   val_angle_Rc   = conv_floats(inp.readline())
   val_angle_eta  = conv_floats(inp.readline())
   val_angle_lamb = conv_floats(inp.readline())
   val_angle_zeta = conv_floats(inp.readline())
   bond_remove    = inp.readline()
   angle_remove   = inp.readline()
else:
   str = ''
   cnt = 0
   for idx, pair in enumerate(two_body):
      cnt += 1
      if (idx < 10):  str += '  [%i]  %2s-%2s   '%(idx,pair[0],pair[1])
      if (idx >= 10): str += '  [%2i] %2s-%2s   '%(idx,pair[0],pair[1])
      if cnt % 5 == 0:
         print(str)
         str = ''
      elif idx+1 == len(two_body):
         print(str)
   idxs = input(' removable bonds : ')
   bond_remove = ''
   for idx in idxs.split():
      idx = int(idx)
      bond_remove += '%2s-%2s/%2s-%2s/'%(two_body[idx][0],two_body[idx][1],two_body[idx][1],two_body[idx][0])

   out = open('stp.dat','w')
   for atom in NameSpec: out.write('%3s'%(atom))
   print()
   print(' ===== Two-body potential parameters ===== ')
   print()
   print('    [1] G1 - simple radial function')
   print('    [2] G2 - G1 multiplied exponential function')
   print('    [3] G3 - G1 multiplied cos function')

   func_bond   = get_int('   - Function type : ',3)
   out.write('\nG %i function'%(func_bond))

   val_bond_Rc = get_floats('   - Rc     : ')
   out.write('\n')
   for val in val_bond_Rc: out.write(' %.4f'%(val))
   nFuncB = len(val_bond_Rc)
   if func_bond == 2:
      val_bond_eta = get_floats('   - eta    : ')
      out.write('\n')
      for val in val_bond_eta: out.write(' %.7f'%(val))
      val_bond_Rs  = get_floats('   - Rs     : ')
      out.write('\n')
      for val in val_bond_Rs: out.write(' %.4f'%(val))
      nFuncB = nFuncB*len(val_bond_eta)*len(val_bond_Rs)
   if func_bond == 3:
      val_bond_kpa = get_floats('   - kappa  : ')
      out.write('\n')
      for val in val_bond_kpa: out.write(' %.7f'%(val))
      nFuncB = nFuncB*len(val_bond_kpa)

   print()
   print()

   angle_remove = ''
   for atom in NameSpec:
      str = ''
      cnt = 0
      for idx, pair in enumerate(two_body):
         cnt += 1
         if (idx < 10):  str += '  [%i]  %2s-%2s-%2s   '%(idx,pair[0],atom,pair[1])
         if (idx >= 10): str += '  [%2i] %2s-%2s-%2s   '%(idx,pair[0],atom,pair[1])
         if cnt % 5 == 0:
            print(str)
            str = ''
         elif idx+1 == len(two_body):
            print(str)
      idxs = input(' removable angles : ')
      for idx in idxs.split():
         idx = int(idx)
         angle_remove += '%2s-%2s-%2s/'%(two_body[idx][0],atom,two_body[idx][1])
         angle_remove += '%2s-%2s-%2s/'%(two_body[idx][1],atom,two_body[idx][0])

   print(' ===== Three-body potential parameters ===== ')
   print()
   print('    [1] G4 - simple radial function')
   print('    [2] G5 - G1 multiplied exponential function')

   func_angle     = get_int('   - Function type : ',2)+3
   out.write('\nG %i function'%(func_angle))
   val_angle_Rc   = get_floats('   - Rc     : ')
   out.write('\n')
   for val in val_angle_Rc: out.write(' %.4f'%(val))

   val_angle_eta  = get_floats('   - eta    : ')
   out.write('\n')
   for val in val_angle_eta: out.write(' %.7f'%(val))

   val_angle_lamb = get_floats('   - lambda : ')
   out.write('\n')
   for val in val_angle_lamb: out.write(' %.2f'%(val))

   val_angle_zeta = get_floats('   - zeta   : ')
   out.write('\n')
   for val in val_angle_zeta: out.write(' %.2f'%(val))
   out.write('\n'+bond_remove+'\n'+angle_remove)
   out.close()

nFuncA = len(val_angle_Rc)*len(val_angle_eta)*len(val_angle_lamb)*len(val_angle_zeta)
for Atom1 in NameSpec:
   exFunc = 0
   out = open('%s.stp'%(Atom1),'w')
   out.write('DESCR\n\nEND DESCR\n\n')
   out.write('ATOM %s\n\n'%(Atom1))
   out.write('ENV %i\n'%(len(NameSpec)))
   for Atom2 in NameSpec: out.write('%s\n'%(Atom2))
   out.write('\nRMIN 0.75d0\n\nSYMMFUNC type=Behler2011\n')
   out.write('nFunc\n')
   for Atom2 in NameSpec:
      check='%2s-%2s'%(Atom1,Atom2)
      if check in bond_remove:
         exFunc += nFuncB
         continue
      for Rc in val_bond_Rc:
         if func_bond == 1: out.write('G=1 type2=%s Rc=%.4f\n'%(Atom2, Rc))
         if func_bond == 2:
            for Rs in val_bond_Rs:
               for eta in val_bond_eta: out.write('G=2 type2=%s eta=%.7f Rs=%.4f Rc=%.4f\n'%(Atom2, eta, Rs, Rc))
         if func_bond == 3:
            for kpa in val_bond_kpa: out.write('G=3 type2=%s kappa=%.7f Rc=%.4f\n'%(Atom2, kpa, Rc))
   for pair in two_body:
      check1 = '%2s-%2s'%(Atom1,pair[0])
      check2 = '%2s-%2s'%(Atom1,pair[1])
      if check1 in bond_remove or check2 in bond_remove:
         exFunc += nFuncA
         continue
#      if check1 in bond_remove or check2 in bond_remove: out.write('# ')
      check1 = '%2s-%2s-%2s'%(pair[0],Atom1,pair[1])
      check2 = '%2s-%2s-%2s'%(pair[1],Atom1,pair[0])
      if check1 in angle_remove or check2 in angle_remove:
         exFunc += nFuncA
         continue
      for Rc in val_angle_Rc:
         for eta in val_angle_eta:
            for lamb in val_angle_lamb:
               for zeta in val_angle_zeta:
                  out.write('G=%i type2=%s type3=%s eta=%.7f '%(func_angle, pair[0], pair[1], eta))
                  out.write('lambda=%.2f zeta=%.2f Rc=%.4f\n'%(lamb, zeta, Rc))
   out.close()
   nFunc = nFuncB*nSpec + nFuncA*len(two_body) - exFunc
   os.system("sed -i 's/nFunc/%i/g' %s.stp"%(nFunc,Atom1))

