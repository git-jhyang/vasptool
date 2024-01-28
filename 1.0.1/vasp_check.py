#!/usr/bin/env python2

import os, time
import lp24vasp as lp

if not os.path.isfile('INCAR'): lp.call_error('bad location','run this code in VASP calculation directory')

out = lp.get_outcar('OUTCAR')
energy = out.get_energy()['E_sigma0']

ts = time.ctime(os.path.getmtime('IBZKPT')).replace(':',' ').split()
te = time.ctime(os.path.getmtime('OUTCAR')).replace(':',' ').split()

month = {'Jan':31, 'Fab':28, 'Mar':31, 'Apl':30, 'May':31, 'Jun':30, 'Jul':31, 'Aug':31, 'Sep':30, 'Oct':31, 'Nov':30, 'Dec':31}

dt_tot = (te[2] - ts[2])*24
if dt_dat < 0: dt_dat += month[ts[1]]*24
dt_tot += te[3] - ts[3] + float(te[4] - ts[4])/60

stat = 'Runneing'
if out.is_end: stat = 'Converged'

print ' %10s   //   Loop: %4d   //   E: %16.6f eV   //   Time: %7.2f hour'%(stat, out.loop, energy, dt_tot)
