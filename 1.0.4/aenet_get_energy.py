#!/usr/bin/env python3

import sys

if len(sys.argv) < 2:
   print('   usage : ',sys.argv[0],' [PREDICT_OUTPUT] ')
   exit(1)

file_pred = sys.argv[1]
inp       = open(file_pred)
isRead    = False
isGetEne1 = False
isGetEne2 = False

for line in inp:
   if 'Energy evaluation' in line: isRead = True
   if not isRead: continue
   if (not isGetEne1) and ('File name' in line):
      inp_xsf = open(line.split()[3])
      ene_inp = inp_xsf.readline().split()[4]
      inp_xsf.close()
      isGetEne1 = True
   if (not isGetEne2) and ('Total energy' in line):
      ene_out = line.split()[3]
      isGetEne2 = True
   if isGetEne1 and isGetEne2:
      print(ene_inp, ene_out)
      isGetEne1 = False
      isGetEne2 = False

inp.close()
