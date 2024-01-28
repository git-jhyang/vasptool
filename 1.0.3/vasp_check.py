#!/usr/bin/env python3

import os, time
import lp34vasp as lp

if not os.path.isfile('INCAR'): lp.call_error('bad location','run this code in VASP calculation directory')
if not os.path.isfile('OUTCAR'): lp.call_error('OUTCAR is missing','Not on running')

loop    = int(os.popen('grep \'energy  without\' OUTCAR | wc -l').read())
eneline = os.popen('grep \'energy  without\' OUTCAR | tail -1').read().split()

if len(eneline) == 0:
   enes = [0.0, 0.0]
else:
   enes = [float(eneline[3]),float(eneline[6])]

endline = os.popen('grep Elapsed OUTCAR').read().split()
isend   = False
if len(endline) != 0: isend = True

if isend:
   dt_tot = float(os.popen('grep Elapsed OUTCAR').read().split()[3])/3600
elif os.path.isfile('IBZKPT'):
   ts = time.ctime(os.path.getmtime('IBZKPT')).replace(':',' ').split()
   te = time.ctime(os.path.getmtime('OUTCAR')).replace(':',' ').split()
   month = {'Jan':31, 'Fab':28, 'Mar':31, 'Apl':30, 'May':31, 'Jun':30}
   month = {'Jul':31, 'Aug':31, 'Sep':30, 'Oct':31, 'Nov':30, 'Dec':31}
   dt_tot = (float(te[2]) - float(ts[2]))*24
   if dt_tot < 0: dt_tot += month[ts[1]]*24
   dt_tot += float(te[3]) - float(ts[3]) + (float(te[4]) - float(ts[4]))/60
else:
   dt_tot = 0.0

stat = 'Runneing'
if isend: stat = 'Converged'

print(' %10s // %4d Loop // Ew.o.S = %15.8f / Esig.0 = %15.8f // %7.2f hour'%(stat, loop, enes[0], enes[1], dt_tot))

