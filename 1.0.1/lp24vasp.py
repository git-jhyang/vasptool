"""
   lp24vasp
   ------------------------------------------------------------------
   python2 module for vasp in/out files written by Yang Jin-Hoon
   ------------------------------------------------------------------
   usage : 

   import lp24vasp as lp
   lp.FUNCTION

   more information can be seen by
   print lp.FUNCINON.__doc__
   ------------------------------------------------------------------
   - functions -
   .assert_num 	.assert_file   .get_str    .get_path
   .chk_len       .call_error	   .angle      .cart_to_dir
   .dir_to_cart   .lvec_rebuild

   - classes -
   .get_outcar    .get_poscar    .get_eigenval     .write_xsf
   .write_poscar
"""

import os
import numpy as np

def call_error(str1, str2, isexit=True):
   """
	call_error(str1, str2, isexit=True)
	------------------------------------------------------------------
	print error message and exit program. if isexit=False, do not exit
	------------------------------------------------------------------
	' !!! ERROR : 'str1' : 'str2'
   """
   print ' !!! ERROR : '+str1+' : '+str2
   if isexit: exit(1)

def assert_num(dat):
   """
	assert_num(dat):
	------------------------------------------------------------------
	check whether given data 'dat' is a number or not.
	if a number, return 'True', else, return 'False'.
   """
   try:
      np.float_(dat)
   except ValueError:
      return False
   else:
      return True

def assert_file(f, str, isopen=True):
   """
	assert_file(filename, str, isopen = True)
	------------------------------------------------------------------
	check whether file 'filename' exists or not.
	if there is no file, call error with 'str' for error check.
	if there is a file and 'isopen=True', return 'open(filename)'.
   """
   if os.path.isfile(f):
      if isopen: return open(f)
   else:
      call_error('file does not exists - '+str, f)

def assert_dat(dat1, dat2):
   """
   assert_dat(dat1, dat2)
   ------------------------------------------------------------------
   check whether 'dat1' and 'dat2' are equal or not.
   if datas are equal, return True, else, return False
   """
   if dat1 == dat1:
      return True
   else:
      return False

def get_str(str, file, n=-1):
   """
	get_str(str, filename, n = -1)
	------------------------------------------------------------------
	open file 'filename' and take a line that contains string 'str'
	and return. if there is any line contins 'str', return 'None'.
	------------------------------------------------------------------
	- options for 'n' -
	n < 0 (default) : return last one line.
	n >= 0          : return 'n'-th line (start from 0) and break loop.
			  if 'n' is larger than number of lines, return 
			  last one line.
	n = 'all'	: return all lines.
   """
   f = assert_file(file, 'get_str', isopen = True)
   buf = []
   for line in f:
      if str in line:
         buf.append(line.split())
         if (type(n) == int) and (np.shape(buf)[0]-1 == n): break
   f.close()
   if np.shape(buf)[0] == 0:
      return None
   elif type(n) == int:
      if n > np.shape(buf)[0]-1 or n < 0: return buf[-1]
      return buf[n]
   elif n == 'all' or n[0] == 'a':
      return buf
   else:
      return buf

def get_path(path, n=0):
   """
	get_path(path, n = 0)
	------------------------------------------------------------------
	if input 'path' include filename at the end, return 'path' without 
	filename. if integer abs('n') is larger than 0, substract 'n' 
	directories from the end in 'path' and return.
   """
   isfile = False
   if '~' in path: path = os.path.expanduser(path)
   if type(n) is not int: call_error('input error - get_path','non-integer \'n\'')
   if os.path.isfile(path): isfile = True
   if isfile:
      path_buf = path.split('/')
      path = ''
      for i in range(len(path_buf)-1): path += path_buf[i]+'/'
   if n == 0: return path
   if (len(path.split('/')) - abs(n) - 1 < 0): return path
   path_buf = path.replace('/',' ').split()
   path = ''
   for i in range(len(path_buf)-abs(n)):
      path += path_buf[i]+'/'
   return path

def chk_len(list, n, fin=True, bool=False):
   """ 
	chk_len(list, n, fin=True, bool=False)
	------------------------------------------------------------------
	check whether 'm' (the first dimension of input array 'list') 
	is larger than input integer 'n'.
	------------------------------------------------------------------
	if 'm' is larger than 'n', return list[n]
        if 'm' is smaller than 'n' and fin = True, exit program
	if 'm' is smaller than 'n' and fin = False, return 'm' 
        if bool is True, fin become False and check length of list
        if length of list is 'n', return True, else return False
   """
   if type(n) is not int: call_error('input error - chk_len','non-integer \'n\'')
   if bool:
      if np.shape(list)[0] == n: return True
      if np.shape(list)[0] != n: return False
   else:
      if np.shape(list)[0] <= n:
         if fin: call_error('input error - chk_len','%2i < %2i'%(np.shape(list)[0]-1,n))
         if not fin: return np.shape(list)[0]
      else:
         return list[n]

def angle(vec1, vec2):
   """
	angle(vec1, vec2)
	------------------------------------------------------------------
	return angle between two vectors, 'vec1' and 'vec2' using arccos.
   """
   if np.shape(vec1) != (3,):
      call_error('input error - angle','dim - %s'%(np.shape(vec1),))
   if np.shape(vec2) != (3,):
      call_error('input error - angle','dim - %s'%(np.shape(vec2),))
   len1 = np.float_(np.dot(vec1, vec1))
   len2 = np.float_(np.dot(vec2, vec2))

   if len1 != 0: vec1 = vec1/len1
   if len2 != 0: vec2 = vec2/len2
   if len1 == 0 or len2 == 0: return np.float_(0)
   return np.float_(np.arccos(np.dot(vec1, vec2)))
       
def dir_to_cart(lvec, coord):
   """
	dir_to_cart(lvec, coord)
	------------------------------------------------------------------
	convert fractional coordinate to cartesian and return.
	------------------------------------------------------------------
	- input -
	lvec  : lattice vector (dim: 3*3)
	coord : fractional coordinate in xyz (dim: n*3)
   """
   if np.shape(lvec) != (3,3):
      call_error('input error - dir_to_cart','dim - %s'%(np.shape(lvec),))
   if np.shape(coord)[1] != 3:
      call_error('input error - dir_to_cart','dim - %s'%(np.shape(coord),))

   latconts = [np.sqrt(np.dot(v, v)) for v in lvec]
   alp  = angle(lvec[1][:], lvec[2][:])
   bet  = angle(lvec[2][:], lvec[0][:])
   gam  = angle(lvec[0][:], lvec[1][:])
   cosa = np.cos(alp)
   cosb = np.cos(bet)
   cosg = np.cos(gam)
   sing = np.sin(gam)
   vol  = np.prod(latconts)*np.sqrt(np.float_(1e0) - cosa**2 - cosb**2 - cosg**2 + \
          np.float_(2)*cosa*cosb*cosg) 
   natoms = np.shape(coord)[0]

   convm = np.zeros([3,3],dtype=np.float_)

   convm[0][:] = [latconts[0], latconts[1]*cosg, latconts[2]*cosb]
   convm[1][:] = [0, latconts[1]*sing, latconts[2]*(cosa - cosb*cosg)/sing]
   convm[2][:] = [0, 0, vol/(latconts[0]*latconts[1]*sing)]

   for i in range(natoms):
      for j in range(3):
         while coord[i][j] >= 1: coord[i][j] -= np.float_(1.)
         while coord[i][j] <  0: coord[i][j] += np.float_(1.)
      coord[i][:] = np.matmul(convm, coord[i][:])
   return coord

def cart_to_dir(lvec, coord):
   """
	cart_to_dir(lvec, coord)
	------------------------------------------------------------------
	convert cartesian coordinate to fractional and return.
	------------------------------------------------------------------
	- input -
	lvec  : lattice vector (dim: 3*3)
	coord : cartesian coordinate in xyz (dim: n*3)
   """
   if np.shape(lvec) != (3,3):
      call_error('input error - dir_to_cart','dim - %s'%(np.shape(lvec),))
   if np.shape(coord)[1] != 3: 
      call_error('input error - dir_to_cart','dim - %s'%(np.shape(coord),))

   latconts = [np.sqrt(np.dot(v, v)) for v in lvec]
   alp  = angle(lvec[1][:], lvec[2][:])
   bet  = angle(lvec[2][:], lvec[0][:])
   gam  = angle(lvec[0][:], lvec[1][:])
   cosa = np.cos(alp)
   cosb = np.cos(bet)
   cosg = np.cos(gam)
   sing = np.sin(gam)
   vol  = np.prod(latconts)*np.sqrt(np.float_(1e0) - cosa**2 - cosb**2 - cosg**2 + \
          np.float_(2)*cosa*cosb*cosg) 
   natoms = np.shape(coord)[0]

   convm = np.zeros([3,3],dtype=np.float_)

   convm[0][:] = [1./latconts[0], -cosg/(latconts[0]*sing), \
                  latconts[1]*latconts[2]*(cosg*cosa - cosb)/(sing*vol)]
   convm[1][:] = [0, 1./(latconts[1]*sing), \
                  latconts[0]*latconts[2]*(cosb*cosg - cosa)/(sing*vol)]
   convm[2][:] = [0, 0, latconts[0]*latconts[1]*sing/vol]

   for i in range(natoms):
      coord[i][:] = np.matmul(convm, coord[i][:])
      for j in range(3):
         while coord[i][j] >= 1: coord[i][j] -= np.float_(1)
         while coord[i][j] <  0: coord[i][j] += np.float_(1)
   return coord

def lvec_rebuild(lvec):
   """
	lvec_rebuild(lvec)
	------------------------------------------------------------------
	rebuild lattice vector 'lvec' and return.
   """
   if np.shape(lvec) != (3,3):
      call_error('input error - lvec_rebuild','dim - %s'%(np.shape(lvec),))
   latconts = [np.sqrt(np.dot(v, v)) for v in lvec]
   alp  = angle(lvec[1][:], lvec[2][:])
   bet  = angle(lvec[2][:], lvec[0][:])
   gam  = angle(lvec[0][:], lvec[1][:])
   cosa = np.cos(alp)
   cosb = np.cos(bet)
   cosg = np.cos(gam)
   
   mat = np.zeros([3,3], dtype=np.float_)

   mat[0][:] = [latconts[0], 0, 0]
   mat[1][:] = [latconts[1]*cosg, np.sqrt(latconts[1]**2 - mat[1][0]**2), 0]
   mat[2][:] = [latconts[2]*cosb, \
                (cosa*latconts[1]*latconts[2] - mat[1][0]*mat[2][0])/mat[1][1], \
                np.sqrt(latconts[2]**2 - mat[2][0]**2 - mat[2][1]**2)]
   return mat

class get_outcar(object):
   """ 
   get_outcar(filename)
   ------------------------------------------------------------------
   read VASP OUTCAR and return various information.
   ------------------------------------------------------------------
   - output -
   .filename : filename of OUTCAR
   .input    : return INCAR variables in python dictionary
      key(integer) :	nsw, nions, ispin, nkpts, nbands
      key(string)  :	lorbit, lwave, lcharg, lvtot, lsorbit
   .atoms    : return atom names 
   .is_band  : 'True' if band structure calculation
   .is_stat  : 'True' if calculation is static self-consistant
   .is_end   : 'True' if calculation finished

   .get_energy()   : return energies in python dictionary
      key(float) : E_fermi, E_free, E_w/oS, E_sigma0
   .get_magnetize() : return n*3 magnetization vector if LORBIT >= 10
   .get_total_magnetize() : return total magnetization
   .get_force() : return n*3 force vector
   """
   def __init__(self, filename):
      self.filename = filename
      f = assert_file(self.filename,'get_outcar', isopen=True)
      f.seek(-1,2)
      self.eof = f.tell()
      f.close()
      self.input = self.get_inputs()
      if self.input['nsw'] >= 1: self.energy = self.get_energy()
      self.atoms = self.get_atoms()
      self.is_band = False
      if get_str('interpolating k-points', self.filename, n=1) is not None:
         self.is_band = True
      self.is_stat = False
      if not self.is_band and self.input['nsw'] == 0:
         self.is_stat = True
      self.loop   = 0 
      self.is_end = self.get_ended()

   def get_inputs(self):
      input = {}
      input['nsw']     = int(chk_len(get_str('  NSW', self.filename, n=0), 2))
      input['nedos']   = int(chk_len(get_str('  NEDOS', self.filename, n=0), 5))
      input['nions']   = int(chk_len(get_str('  NIONS', self.filename, n=0), 11))
      input['ispin']   = int(chk_len(get_str('  ISPIN', self.filename, n=0), 2))
      input['lorbit']  = int(chk_len(get_str('  LORBIT', self.filename, n=0), 2))
      input['lwave']   = chk_len(get_str('  LWAVE', self.filename, n=0), 2)
      input['lcharg']  = chk_len(get_str('  LCHARG', self.filename, n=0), 2)
      input['lvtot']   = chk_len(get_str('  LVTOT', self.filename, n=0), 2)
      input['lsorbit'] = chk_len(get_str('  LSORBIT', self.filename, n=0), 2)
      input['ncollin'] = chk_len(get_str('  LNONCOLLINEAR', self.filename, n=0), 2)
      input['nkpts']   = int(chk_len(get_str('  NKPTS', self.filename, n=0), 3))
      input['nbands']  = int(chk_len(get_str('  NBANDS=', self.filename, n=0), 14))
      input['nelect']  = float(chk_len(get_str('  NELECT =', self.filename, n=0), 2))
      input['atoms']   = self.get_atoms()
      input['nspec']   = get_str('ions per type =', self.filename, n=0)[4:]
      for idx,num in enumerate(input['nspec']): input['nspec'][idx] = int(num)
      return input

   def get_energy(self):
      enes = {}
      enes['E_fermi']  = np.float_(chk_len(get_str('E-fermi', self.filename), 2))
      enes['E_free']   = np.float_(chk_len(get_str('free  energy', self.filename), 4))
      enes['E_w/oS']   = np.float_(chk_len(get_str('energy  without', self.filename), 3))
      enes['E_sigma0'] = np.float_(chk_len(get_str('energy  without', self.filename), 6))
      return enes

   def get_ended(self):
      isend = False
      lines = get_str('EDIFF is reached', self.filename, n='all')
      self.loop = np.shape(lines)[0]
      if self.is_stat:
         if self.loop == 1: isend = True
      else:
         line = get_str('reached required accuracy', self.filename)
         if line is not None: isend = True
      return isend

   def get_atoms(self):
      f = open(self.filename)
      buf = []
      atoms = []
      while True:
         line = f.readline()
         if 'VRHFIN' in line: buf.append(line)
         if f.tell() > self.eof*0.1: break
      f.close()     
      for line in buf:
         line = line.replace(':',' ').replace('=',' ').split()
         atoms.append(line[1])
      return atoms

   def get_lattice(self):
      f = open(self.filename)
      f.seek(0)
      loc_lvec = None
      while (f.tell() <= self.eof):
         line = f.readline()
         if 'direct lattice vectors' in line: loc_lvec = f.tell()
      if loc_lvec is None: call_error('reading error - get_lattice','lattice is missing')
      f.seek(loc_lvec)
      lvec = np.zeros([2,3,3],dtype=np.float_)
      for i in range(3):
         line = f.readline().split()
         for j in range(3):
            lvec[0][i][j] = np.float_(line[j])
            lvec[1][i][j] = np.float_(line[j+3])
      return lvec

   def get_magnetize(self):
      iscollin = True
      mag = np.zeros([self.input['nions'], 3],dtype=np.float_)
      if self.input['ncollin'] == 'T': iscollin = False
      if not (self.input['ispin'] == 2 or (not iscollin)): 
         call_error('input error - get_magnetize','spin polarization required')
      if self.input['lorbit'] < 10: 
         call_error('input error - get_magnetize','LORBIT >= 10 required')
      f = open(self.filename)
      f.seek(0)
      magx = None
      while (f.tell() <= self.eof):
         line = f.readline()
         if 'magnetization (x)' in line: magx = f.tell()
         if 'magnetization (y)' in line: magy = f.tell()
         if 'magnetization (z)' in line: magz = f.tell()
      if magx is None: call_error('reading error - get_magnetize','polarization is missing')
      f.seek(magx)
      for i in range(3): f.readline() # empty
      for i in range(self.input['nions']):
         if self.input['ispin'] == 2: mag[i][2] = np.float_(f.readline().split()[4])
         if iscollin: mag[i][0] = np.float_(f.readline().split()[4])
      if self.input['ispin'] == 2: return mag
      f.seek(magy)
      for i in range(3): buf = f.readline() # empty
      for i in range(self.input['nions']):
         mag[i][1] = np.float_(f.readline().split()[4])
      f.seek(magz)
      for i in range(3): buf = f.readline() # empty
      for i in range(self.input['nions']):
         mag[i][2] = np.float_(f.readline().split()[4])
      return mag

   def get_total_magnetize(self):
      iscollin = True
      if self.input['ncollin'] == 'T': iscollin = False
      if not (self.input['ispin'] == 2 or (not iscollin)):
         call_error('input error - get_total_magnetize','spin polarization required')
      if not iscollin: mag = np.zeros([3], dtype=np.float_)
      f = open(self.filename)
      f.seek(0)
      line_data = None
      while (f.tell() <= self.eof):
         line = f.readline()
         if ' number of electron ' in line: line_data = line
      if line_data is None: call_error('reading error - get_total_magnetize','magnetization is missing')
      if iscollin: mag = np.float_(line_data.split()[5])
      if not iscollin:
         for i in range(3): mag[i] = np.float_(line_data.split()[i+5])
      return mag

   def get_force(self):
      force = np.zeros([self.input['nions'], 3], dtype=np.float_)
      f = open(self.filename)
      f.seek(0)
      loc_force = None
      while (f.tell() <= self.eof):
         line = f.readline()
         if 'TOTAL-FORCE (eV/Angst)' in line: loc_force = f.tell()
      if loc_force is None: call_error('reading error - get_force','force is missing')
      f.seek(loc_force)
      f.readline() # buffer line
      for i in range(self.input['nions']):
         line = f.readline().split()
         for j in range(3): force[i][j] = np.float_(line[j+3])
      return force

class get_doscar(object):
   """
   get_doscar(filename, input)
   ------------------------------------------------------------------
   read VASP DOSCAR and return dataset.
   input : INCAR tags get from lp.get_outcar().input
   ------------------------------------------------------------------
   - output -
   .tdos()           : total DOS with E-EF
   .get_pdos()       : lm and site projected DOS (PDOS) if LORBIT is not 0
   .write_dos(n,d,h) : write dos datafile. requires 'filename' f,
                       'dos data' d, and 'header' h
   .nedos    : number of grid
   .natoms   : number of atoms, required by PDOS
   .ef       : fermi level
   .header_T : header of total DOS
   .header_P : header of PDOS
   """
   def __init__(self, filename, inp):
      self.filename = filename
      self.f = assert_file(self.filename, 'get_doscar', isopen=True)
      self.natoms = int(chk_len(self.f.readline().split(), 0))
      for i in range(4): self.f.readline() # read buffer lines
      line = self.f.readline().split() # data line
      self.nedos = int(chk_len(line, 2)) - 1
      self.ef    = np.float_(chk_len(line, 3))
      if (not assert_dat(self.nedos, inp['nedos'])):
         call_error('input error','mismatch btw DOSCAR and OUTCAR')
      self.ispin    = inp['ispin']
      self.iscollin = True
      if inp['ncollin'] == 'T': self.iscollin = False
      self.lorbit   = inp['lorbit']
      self.tdos     = self.__get_tdos()
      self.header_T = ['E-EF','Tot']
      if self.ispin == 2:
         self.header_T.append(self.header_T[1]+'_D')
         self.header_T[1] += '_U'
      if self.lorbit == 0:
         return
      elif self.lorbit < 11:
         header_tmp = ['s','p','d']
      else:
         header_tmp = ['s','px','py','pz','dxy','dyz','dz2','dxz','dx2y2']
      self.header_P = ['E-EF']
      if self.ispin == 2:
         for str in header_tmp:
            self.header_P.append(str+'_U')
            self.header_P.append(str+'_D')
      elif not self.iscollin:
         for str in header_tmp:
            self.header_P.append(str+'-Tot')
            self.header_P.append(str+'-mx')
            self.header_P.append(str+'-my')
            self.header_P.append(str+'-mz')
   def __get_tdos(self):
      tdos = np.zeros([self.nedos,self.ispin+1],dtype=np.float_)
      self.f.readline() # first line, remove
      for i in range(self.nedos):
         line = self.f.readline().split()
         tdos[i][0] = np.float_(chk_len(line,0)) - self.ef
         for j in range(self.ispin): tdos[i][1+j] = np.float_(chk_len(line,1+j))
      return tdos
   def get_pdos(self):
      """
      .get_pdos() : return PDOS with dimension of NATOMS*NEDOS*DATA
      """
      if self.lorbit == 0: call_err('LORBIT is 0','No partial dos')
      ncol = len(self.header_P)
      pdos = np.zeros([self.natoms,self.nedos,ncol],dtype=np.float_)
      for iAtom in range(self.natoms):
         self.f.readline() # data line
         self.f.readline() # first line, remove
         for idos in range(self.nedos):
            line = self.f.readline().split()
            for idx,dat in enumerate(line): pdos[iAtom][idos][idx] = np.float_(dat)
            pdos[iAtom][idos][0] -= self.ef
      return pdos
   def write_dos(self, fname_out, dos, header):
      """
      .write_dos(f,d,h) : write dos datafile with filename 'f' using data 'd'
                          and data description header 'h'
      """
      f = open(fname_out,'w')
      f.write('#%9s'%(header[0]))
      for str in header[1:]: f.write(' %10s'%(str))
      f.write('\n')
      for idos in range(self.nedos):
         for dat in dos[idos]:
            f.write('%10.4f '%(dat))
         f.write('\n')

class get_poscar(object):
   """
	get_poscar(filename)
	------------------------------------------------------------------
	read VASP POCSCAR and return dataset.
	------------------------------------------------------------------
	- output -
	.filename : filename of POSCAR
	.des      : description of POSCAR at the first line
	.lvec     : lattice vector
	.style    : POSCAR style, 4 or 5
	.natoms   : total number of atoms
	.atoms    : atom species (empty array returned if .style = 4)
	.nspec    : number of each atom species
	.isselc   : True if selective dynamics used
	.iscart   : True if coordinate is cartesian
	.coord    : atomic coordination in xyz (dim : .natoms*3)
	.selc     : if isselc is True, return selective dynamics flags
   """
   def __init__(self, filename):
      self.filename = filename
      f = assert_file(self.filename, 'get_poscar', isopen = True)
      self.des = f.readline()
# get lattice constant
      buf = f.readline().split()
      if assert_num(buf[0]): latcont = np.float_(buf[0])
      if not assert_num(buf[0]): call_error('reading error - get_poscar','lcont - '+buf[0])

# get lattice vector
      self.lvec = np.zeros((3,3), dtype=np.float_)
      for ii in range(3):
         buf = f.readline().split()
         for ij, num in enumerate(buf):
            if ij == 3: break
            if assert_num(num): self.lvec[ii][ij] = np.float_(num)*latcont
            if not assert_num(num): call_error('reading error - get_poscar','lvec - '+num)

# get atom types and numbers
      buf = f.readline().split()
      self.style = 5
      self.natoms = 0
      if assert_num(buf[0]): self.style = 4
      if self.style == 5:
         self.atoms = buf
         self.nspec = f.readline().split()
      else:
         self.atoms = []
         self.nspec = buf
      for i,j in enumerate(self.nspec):
         if not assert_num(j): call_error('reading error - get_poscar','nspec - '+j)
         self.nspec[i] = int(j)
         self.natoms += int(j)
      self.isselc = False
      self.iscart = False

# get POSCAR type (selective dynamics // fractional // cartesian)
      buf = f.readline()
      if buf[0] == 'S':
         self.isselc = True
         buf = f.readline()
      if buf[0] == 'C' or buf[0] == 'c':
         self.iscart = True
      elif buf[0] == 'D' or buf[0] == 'd':
         self.iscart = False
      else:
         call_error('reading error - get_poscar','type - '+buf)

# get coordinate
      self.coord = np.zeros((self.natoms,3), dtype=np.float_)
      if self.isselc: self.selc = []
      for i in range(self.natoms):
         buf = f.readline().split()
         for j in range(3):
            if assert_num(buf[j]):
               bufdat = np.float_(buf[j])
               self.coord[i][j] = bufdat
            if not assert_num(buf[j]):
               call_error('reading error - get_poscar','coo - %3i'%(i)+buf[j])
         if self.isselc:
            if len(buf) != 6: call_error('reading error - get_poscar','selc - %3i'%(i))
            self.selc.append(buf[3:6])

class get_eigenval(object):
   """
	get_eigenval(filename)
	------------------------------------------------------------------
	read VASP EIGENVAL and return dataset.
	------------------------------------------------------------------
	- output -
	.filename : filename of EIGENVAL
	.ispin    : ISPIN
	.nelect   : number of electrons
	.nkpts    : number of kpoints
	.nbands   : number of bands
	.kpts     : kpoint coordination in xyz (dim: .nkpts*3)
	.weight   : weight of kpoint (dim: .nkpts)
	.eigens   : eigenvalues (dim: .nkpts*.nbands*.ispin)
   """
   def __init__(self, filename):
      self.filename = filename
      f = assert_file(self.filename, 'get_eigenval', isopen=True)
#      self.__isfile_eigenval()

# get information
      line = f.readline().split()
      self.ispin = int(chk_len(line, 3))
      for i in range(0,4): f.readline() # buffer line
      line = f.readline().split()
      self.nelect = int(chk_len(line,0))
      self.nkpts  = int(chk_len(line,1))
      self.nbands = int(chk_len(line,2))
      self.eigens = np.zeros([self.nkpts, self.nbands, self.ispin], dtype=np.float_)
      self.kpts   = np.zeros([self.nkpts,3], dtype=np.float_)
      self.weight = np.zeros([self.nkpts], dtype=np.float_)

# get eigenvalues
      for j in range(self.nkpts):
         f.readline() # empty line
         line = f.readline().split()
         self.weight[j] = np.float_(chk_len(line,3))
         for i in range(3): self.kpts[j][i] = np.float_(line[i])
         for i in range(self.nbands):
            line = f.readline().split()
            if len(line) < self.ispin+1:
               call_error('reading error - get_eigenval',self.filename)
            for k in range(self.ispin): self.eigens[j][i][k] = np.float_(line[k+1])
      f.close()
#   def __isfile_eigenval(self):
      
class write_xsf(object):
   """
	write_xsf(filename, dat, ispbc=True, vecs=None)
	------------------------------------------------------------------
	get dataset 'dat' and vector 'vecs = vec' then write XCrysDen 
	coordination file (.xsf) named 'filename'.
	------------------------------------------------------------------
	- requirements (output <- input) -
	.des    <- dat.des    : description of file at first line
	.lvec   <- dat.lvec   : lattice vector
	.natoms <- dat.natoms : total number of atoms
	.atoms  <- dat.atoms  : atom species
	.nspec  <- dat.nspec  : number of each atom species
	.coord  <- dat.coord  : atomic coordination in xyz (dim: .natoms*3)
	.iscart <- dat.iscart : if dat.coord is fractional, False
	.force  <- vecs = vec : (optional) vector, drawn in VESTA (dim: .natoms*3)
   """
   def __init__(self, filename, dat, ispbc=True, vecs=None):
      if not '.xsf' in filename: filename = filename+'.xsf'
      f           = open(filename,'w')
      self.des    = '# '+dat.des
      self.ispbc  = ispbc
      self.lvec   = dat.lvec
      self.natoms = dat.natoms
      self.atom   = dat.atoms
      self.nspec  = dat.nspec
      atoms_l     = [at for ii,at in enumerate(self.atom) for i in range(self.nspec[ii])]
      self.iscart = dat.iscart
      self.coord  = dat.coord
      if not self.iscart: self.coord = dir_to_cart(self.lvec, self.coord)
      self.force  = np.zeros(np.shape(self.coord))
      if vecs is not None: self.force = vecs
      if np.shape(self.force) != np.shape(self.coord):
         call_error('input error - write_xsf', \
                    'dim - %s %s'%(np.shape(self.force),np.shape(coord),))
      f.write(self.des+'\n')
      if (self.ispbc):
         f.write('CRYSTAL\n')
         f.write('PRIMVEC\n')
         for i in self.lvec: f.write(''.join(['%19.15f '%(j) for j in i])+'\n')
         f.write('PRIMCOORD\n')
         f.write('%i %i\n'%(self.natoms, 1))
      else:
         f.write('ATOMS\n')
      for i in range(self.natoms):
         f.write('%2s'%(atoms_l[i]))
         f.write(''.join(['%20.15f '%(j) for j in self.coord[i][:]]))
         if vecs is not None: f.write(''.join(['%14.9f '%(j) for j in self.force[i][:]]))
         f.write('\n')
      f.close()

class write_poscar(object):
   """
	write_poscar(filename, dat)
	------------------------------------------------------------------
	get dataset 'dat' then write VAPS coordination POSCAR named 'filname'.
	------------------------------------------------------------------
	- requirements (output <- input) -
	.des    <- dat.des    : description of file at first line
	.lvec   <- dat.lvec   : lattice vector
	.natoms <- dat.natoms : total number of atoms
	.style  <- dat.style  : POSCAR style (4 or 5)
	.atoms  <- dat.atoms  : atom species
	.nspec  <- dat.nspec  : number of each atom species
	.isselc <- dat.isselc : whether write selective dynamics or not
	.iscart <- dat.iscart : if dat.coord is fractional, False
	.coord  <- dat.coord  : atomic coordination in xyz (dim: .natoms*3)
	.selc   <- dat.selc   : selective dynamics flags
   """
   def __init__(self, filename, dat):
      if not (('POSCAR' in filename) or ('CONTCAR' in filename)):
         filename = 'CONTCAR'+filename
      f           = open(filename,'w')
      self.des    = dat.des
      self.lvec   = dat.lvec
      self.natoms = dat.natoms
      self.style  = dat.style
      self.atom   = dat.atoms
      self.nspec  = dat.nspec
      self.isselc = dat.isselc
      self.iscart = dat.iscart
      self.coord  = dat.coord
      if self.isselc: self.selc   = dat.selc
      if self.iscart: self.coord = cart_to_dir(self.lvec, self.coord)
      if self.isselc and (np.shape(self.selc) != np.shape(self.coord)):
         call_error('input error - write_poscar','dim - %s %s'%(\
			np.shape(self.selc), np.shape(self.coord),))
      f.write(self.des)
      f.write('%20.15f\n'%(np.float(1)))
      for i in self.lvec: f.write(''.join(['%20.15f '%(j) for j in i])+'\n')
      if self.style == 5:
         if np.shape(self.atom) != np.shape(self.nspec):
            call_error('input error - write_poscar','dim - %s %s'%(\
			np.shape(self.atom), np.shape(self.nspec),))
         f.write(''.join(['%4s'%(i) for i in self.atom])+'\n')
      f.write(''.join(['%4i'%(i) for i in self.nspec])+'\n')
      if self.isselc: f.write('Selective Dynamics\n')
      f.write('Direct\n')
      for i in range(self.natoms):
         f.write(''.join(['%20.15f '%(j) for j in self.coord[i][:]]))
         if self.isselc: f.write(''.join(['%3s'%(j) for j in self.selc[i][:]]))
         f.write('\n')
      f.close()

