#!/usr/bin/env python3

import sys, os
import lp34vasp as lp
import random as rnd

def call_help(str1='',str2=''):
   if (str1 != ''): lp.call_error(str1, str2, isexit=False)
   print(' usage   : '+sys.argv[0]+' [DIRECTORY] [-v RATE] [-t RATE]')
   print(' default : [DIRECTORY] = None')
   print('           [-v RATE]   = 0.0 ')
   print('           [-t RATE]   = 0.1 ')
   print(' options : [-h] : print this help')
   exit(1)

path = './'
vrate = 0.0
trate = 0.1

if (len(sys.argv) > 1):
   for idx,arg in enumerate(sys.argv[1:]):
      if arg == '-h': call_help()
      if os.path.isdir(arg): path = arg
      if arg == '-v': vrate = float(sys.argv[idx+2])
      if arg == '-t': trate = float(sys.argv[idx+2])
else:
   call_help('input error','need more than 1 arg')

if path[-1] != '/': path+='/'

fTrain = [f for f in os.listdir(path) if os.path.isfile(path+f)]
fVal   = []
fTest  = []

nFile  = len(fTrain)
nVal   = int(nFile*vrate)
nTest  = int(nFile*trate)

rnd.seed()
while (True):
   if len(fVal) < nVal:
      idx = int(rnd.random()*len(fTrain))
      fVal.append(fTrain.pop(idx))
   elif len(fTest) < nTest: 
      idx = int(rnd.random()*len(fTrain))
      fTest.append(fTrain.pop(idx))
   else:
      break

nTrain = nFile - nTest - nVal

print(' Directory "'+path+'" includes %i files '%(nFile))
print('   Train set      : %i'%(nTrain))
if nVal != 0: print('   Validation set : %i (nTrain + nVal = %i)'%(nVal, nTrain+nVal))
print('   Test set       : %i'%(nTest))

out = open('list.train','w')
for f in fTrain: out.write(path+f+'\n')
out.close()

if nVal != 0:
   out = open('list.validation','w')
   for f in fVal: out.write(path+f+'\n')
   out.close()

out = open('list.test','w')
for f in fTest: out.write(path+f+'\n')
out.close()
