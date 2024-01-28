#!/usr/bin/env python3

import os, time
import pyutil.util as pu

if not os.path.isfile('INCAR'): 
    pu.call_error('bad location','run this code in VASP calculation directory')
if not os.path.isfile('OUTCAR'): 
    pu.call_error('OUTCAR is missing','Not on running')

loop    = int(os.popen('grep \'energy  without\' OUTCAR | wc -l').read())
eneline = pu.get_line('energy  without','OUTCAR')

if len(eneline) == 0:
    enes = [0.0, 0.0]
else:
    enes = [float(eneline[3]),float(eneline[6])]

endline = pu.get_line('Elapsed','OUTCAR')
isend   = False
if len(endline) != 0: 
    isend = True

stat = 'Runneing'

if isend:
    stat = 'Converged'
    dt_tot = float(endline[3])/3600
elif os.path.isfile('IBZKPT'):
    ts = time.ctime(os.path.getmtime('IBZKPT')).replace(':',' ').split()
    te = time.ctime(os.path.getmtime('OUTCAR')).replace(':',' ').split()
    month = {'Jan':31, 'Fab':28, 'Mar':31, 'Apl':30, 'May':31, 'Jun':30}
    month = {'Jul':31, 'Aug':31, 'Sep':30, 'Oct':31, 'Nov':30, 'Dec':31}
    dt_tot = (float(te[2]) - float(ts[2]))*24
    if dt_tot < 0: 
        dt_tot += month[ts[1]]*24
    dt_tot += float(te[3]) - float(ts[3]) + (float(te[4]) - float(ts[4]))/60
else:
    dt_tot = 0.0


print(' %10s // %4d Loop // Ew.o.S = %15.8f /'%(stat, loop, enes[0]), end='')
print(' Esig.0 = %15.8f // %7.2f hour'%(enes[1], dt_tot))

