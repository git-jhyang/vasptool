#!/usr/local/bin/python3.6

import os, socket

def select_option(n,chr):
   print()
   chr_all = chr+' [1-%i] : '%(n)
   while True:
      try:
         idx = int(input(chr_all))
      except ValueError:
         continue
      if idx > 0 and idx < n+1:
         return idx-1

path_all = ['/APP/enhpc/CODES/VASP/bin',
'/APP/enhpc/CODES/VASP/V_5.4/vasp.5.4.1/bin',
'/APP/enhpc/CODES/LAMMPS/bin']

path_mpi  = '/APP/enhpc/mpi/openmpi-intel/bin'

rigel = ['all.q','short.q']
rigel_info = ['all.q(12)*11','short.q(12)*3']
deneb = ['all.q','short.q','E5_2680']
deneb_info = ['all.q(20)*3','short.q(20)*2','E5_2680(28)*4']

# get que name
print()
hostname = socket.gethostname()

k = 0
if hostname == 'nodemaster01':
   for str in rigel_info:
      k += 1
      print('  [%i] %s'%(k,str))
   que = rigel[select_option(len(rigel),' Select QUE name')]
elif hostname == 'nodemaster02':
   for str in deneb_info:
      k += 1
      print('  [%i] %s'%(k,str))
   que = deneb[select_option(len(deneb),' Select QUE name')]
else:
   print(" !!! ERROR : Wrong HOSTNAME '%s' !!!"%(hostname))
   exit(1)

# get list of programs
print()

k=0
prg_all = []
for i in range(0,len(path_all)):
   prgs = os.popen('ls %s'%(path_all[i])).read().split()
   print (' Executions in path : %s'%(path_all[i]))
   for j in range(0,len(prgs)):
      k+=1
      print('  [%i] %s'%(k,prgs[j]))
      prg_all.append(prgs[j])

k = select_option(len(prg_all),' Select execution')
prg = prg_all[k]

for i in range(0,len(path_all)):
   prgs = os.popen('ls %s'%(path_all[i])).read().split()
   k -= len(prgs)
   if k < 0:
      prg_path = path_all[i]
      break

isvasp = False
for str in prg_path.split('/'):
   if str == 'VASP':
      isvasp = True
      break

# get job name
print()
job_name = input(' Your job name : ')

# get # of cores 
print()
while True:
   try:
      num_core = int(input(' Number of cores : '))
   except ValueError:
      continue
   else:
      break

# write que file
O = open('q'+job_name,'w')

O.write("#!/bin/sh\n")
O.write("#$ -S /bin/sh\n")
O.write("#$ -cwd\n")
O.write("#$ -V\n")
O.write("#$ -v MPIHOME=%s\n"%(path_mpi))
O.write("#$ -v PRGHOME=%s\n"%(prg_path))
O.write("#$ -j y\n")
O.write("#$ -q %s\n\n"%(que))

O.write("#$ -N %s\n"%(job_name))
O.write("#$ -pe mpich_ef %i\n\n"%(num_core))
O.write("#$ -o log.out\n\n")

O.write("$MPIHOME/mpirun -np $NSLOTS -machinefile $TMPDIR/machines $PRGHOME/%s "%(prg))

if isvasp:
   O.write("> vasp.out\n")
else:
   O.write("< lammps.in > lammps.out\n")

O.close()

# ask submit
str = input(' Want you submit your job? [y/n] (default : y) : ')
if len(str) == 0:
   os.system('qsub q%s'%(job_name))
elif str[0] == 'N' or str[0] == 'n':
   exit(1)
else:
   os.system('qsub q%s'%(job_name))
