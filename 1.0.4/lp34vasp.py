"""
   lp34vasp
   ------------------------------------------------------------------
   python3 module for vasp in/out files written by Yang Jin-Hoon
   ------------------------------------------------------------------
   usage : 

   import lp34vasp as lp
   lp.FUNCTION

   more information can be seen by
   print lp.FUNCINON.__doc__
   ------------------------------------------------------------------
   - functions -
   .call_error    .assert_num    .assert_file   .assert_dat
   .get_str       .get_path      .chk_len       .angle      
   .dir_to_cart   .cart_to_dir   .lvec_rebuild

   - classes -
   .get_outcar    .get_doscar    .get_eigenval  
   .get_poscar    .write_poscar  .get_xsf       .write_xsf
"""

import os
import numpy as np

def call_error(str1, str2, isexit=True):
   """
   call_error(str1, str2, isexit=True)
   ------------------------------------------------------------------
   print error message and exit program. if isexit=False, do not exit
   ------------------------------------------------------------------
   if isexit is True,  print ' !!! ERROR : 'str1' : 'str2'
   if isexit is False, print ' !!! WARNING : 'str1' : 'str2'
   """
   if isexit: 
      print(' !!! ERROR : '+str1+' : '+str2)
      exit(1)
   else:
      print(' !!! WARNING : '+str1+' : '+str2)

def assert_num(dat):
   """
   assert_num(dat):
   ------------------------------------------------------------------
   check whether given data 'dat' is a number or not.
   if a number, return 'True', else, return 'False'.
   """
   try:
      float(dat)
   except ValueError:
      return False
   else:
      return True

def assert_file(f_name, str='', isopen=True):
   """
   assert_file(f_name, str, isopen = True)
   ------------------------------------------------------------------
   check whether file 'f_name' exists or not.
   if there is no file, call error with 'str' for error check.
   if there is a file and 'isopen=True', return 'open(f_name)'.
   """
   if os.path.isfile(f_name):
      if isopen: return open(f_name)
      return True
   else:
      if isopen: call_error(f_name+' is missing',str)
      return False

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

def get_str(str, f_name, n=0):
   """
   get_str(str, f_name, n=0)
   ------------------------------------------------------------------
   open file 'f_name' and take a line that contains string 'str'
   and return. if there is any line contins 'str', return 'None'.
   ------------------------------------------------------------------
   - options for 'n' -
   n = 0 (default) : return last one line.
   n > 0           : return 'n'-th line.
                     if 'n' is larger than number of lines, return 
                     last one line.
   n = 'all'       : return all lines.
   """
   assert_file(f_name,'get_str',isopen=False)
   comm = 'grep \'%s\' %s '%(str, f_name)
   nLine = int(os.popen(comm+'| wc -l').read())
   if (nLine == 0): return None
   if (type(n) == int) or (type(n) == float):
      iLine = int(n)
      if iLine <= 0 or iLine >= nLine:
         option = '| tail -1'
      else:
         option = '| head -%i | tail -1'%(iLine)
      line = os.popen(comm+option).read().split()
      return line
   elif type(n) == str:
      if (n[0] != 'a'): call_error('get_str - wrong input for \'n\'',' n : '+n)
      buf = os.popen(comm).read().split('\n')
      lines = []
      for line in buf[:-1]: lines.append(line.split())
      return lines
   else:
      call_error('get_str - wrong input for \'n\'','n : '+n)

def get_path(path, n=0):
   """
   get_path(path, n=0)
   ------------------------------------------------------------------
   if input 'path' include filename at the end, return 'path' without 
   filename. if integer n is larger than 0, substract 'n' directories
   from the end in 'path' and return.
   """
   isfile = False
   if '~' in path: path = os.path.expanduser(path)
   if type(n) is not int: n=0
   if n < 0: n = 0
   isfile = assert_file(path, 'get_path', isopen=False)
   if isfile:
      path_buf = path.split('/')
      path = ''
      for dir in path_buf[:-1]: path += dir+'/'
   if n == 0: return path
   if (len(path.split('/')) - n - 1 < 0): return path
   path_buf = path.replace('/',' ').split()
   path = ''
   for dir in path_buf[:-n]: path += dir+'/'
   return path

def chk_len(list, n, fin=True, bool=False):
   """ 
   chk_len(list, n, fin=True, bool=False)
   ------------------------------------------------------------------
   check whether 'm' (the first dimension of input array 'list') 
   is larger than input integer 'n'.
   ------------------------------------------------------------------
   m > n              : return list[n]
   m < n && fin = T   : exit program
   m < n && fin = F   : return 'm' 
   bool = T && m == n : return True
   bool = T && m != n : return False
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
   len1 = np.float_(np.sqrt(np.dot(vec1, vec1)))
   len2 = np.float_(np.sqrt(np.dot(vec2, vec2)))

   if len1 != 0: vec1 = vec1/len1
   if len2 != 0: vec2 = vec2/len2
   if len1 == 0 or len2 == 0: return np.float_(0)
   return np.arccos(np.sqrt(np.dot(vec1, vec2)))
       
def dir_to_cart(lvec, coo):
   """
   dir_to_cart(lvec, coo)
   ------------------------------------------------------------------
   convert fractional coordinate to cartesian and return.
   ------------------------------------------------------------------
   - input -
   lvec : lattice vector (dim: 3*3)
   coo  : fractional coordinate in xyz (dim: n*3)
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

   natoms  = np.shape(coord)[0]
   convm   = np.zeros([3,3],dtype=np.float_)
   coo_out = coo

   convm[0][:]   = [latconts[0], latconts[1]*cosg, latconts[2]*cosb]
   convm[1][1:2] = [latconts[1]*sing, latconts[2]*(cosa - cosb*cosg)/sing]
   convm[2][2]   = vol/(latconts[0]*latconts[1]*sing)

   for i in range(natoms):
      for j in range(3):
         while coo_out[i][j] >= 1: coo_out[i][j] -= np.float_(1.)
         while coo_out[i][j] <  0: coo_out[i][j] += np.float_(1.)
      coo_out[i][:] = np.matmul(convm, coo_out[i][:])
   return coo_out

def cart_to_dir(lvec, coo):
   """
   cart_to_dir(lvec, coo)
   ------------------------------------------------------------------
   convert cartesian coordinate to fractional and return.
   ------------------------------------------------------------------
   - input -
   lvec : lattice vector (dim: 3*3)
   coo  : cartesian coordinate in xyz (dim: n*3)
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
   vol  = np.prod(latconts)*np.sqrt(np.float_(1.) - cosa**2 - cosb**2 - cosg**2 + \
          np.float_(2)*cosa*cosb*cosg) 

   natoms  = np.shape(coord)[0]
   convm   = np.zeros([3,3],dtype=np.float_)
   coo_out = coo

   convm[0][:]   = [1./latconts[0], -cosg/(latconts[0]*sing), \
                    latconts[1]*latconts[2]*(cosg*cosa - cosb)/(sing*vol)]
   convm[1][1:2] = [1./(latconts[1]*sing), \
                    latconts[0]*latconts[2]*(cosb*cosg - cosa)/(sing*vol)]
   convm[2][2]   = latconts[0]*latconts[1]*sing/vol

   for i in range(natoms):
      coo_out[i][:] = np.matmul(convm, coo[i][:])
      for j in range(3):
         while coo_out[i][j] >= 1: coo_out[i][j] -= np.float_(1.)
         while coo_out[i][j] <  0: coo_out[i][j] += np.float_(1.)
   return coo_out

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

   mat[0][0]   = latconts[0]
   mat[1][0:1] = [latconts[1]*cosg, np.sqrt(latconts[1]**2 - mat[1][0]**2)]
   mat[2][:]   = [latconts[2]*cosb, \
                  (cosa*latconts[1]*latconts[2] - mat[1][0]*mat[2][0])/mat[1][1], \
                  np.sqrt(latconts[2]**2 - mat[2][0]**2 - mat[2][1]**2)]
   return mat

class get_outcar(object):
   """ 
   get_outcar(f_name)
   ------------------------------------------------------------------
   read VASP OUTCAR and return various information.
   ------------------------------------------------------------------
   - output -
   .f_name  : filename of OUTCAR
   .input   : return INCAR variables in PyDic
         keys (integer) : nsw, nions, ispin, nkpts, nbands
         keys (string)  : lorbit, lwave, lcharg, lvtot, lsorbit
   .is_band : 'True' if band structure calculation
   .is_stat : 'True' if calculation is static self-consistant
   .loop    :  Elapsed loop
   .is_end  : 'True' if calculation finished

   .get_energy      : return calculated VASP total energy in PyDic
         keys : E_fermi, E_free, E_w/oS, E_sigma0
   .get_lattice()   : return 2*3*3 lattice vector, real and reciprocal
   .get_magnetize() : return n*3 magnetization vector if LORBIT >= 10
   .get_total_magnetize() : return total magnetization
   .get_force() : return n*3 force vector
   """
   def __init__(self, f_name):
      self.f_name = f_name
      assert_file(self.f_name,'get_outcar', isopen=False)
      f = open(self.f_name,'a')
      self.eof = f.tell() - 1
      f.close()
      self.input  = self.get_inputs()
      self.energy = self.get_energy()
      self.is_band = False
      self.is_stat = False
      if get_str('interpolating k-points', self.f_name, n=1) is not None: self.is_band = True
      if not self.is_band and self.input['nsw'] == 0: self.is_stat = True
      self.loop   = 0 
      self.is_end = self.get_ended()

   def get_inputs(self):
      input = {}
      input['nsw']     = int(chk_len(get_str('  NSW', self.f_name, n=1), 2))
      input['nedos']   = int(chk_len(get_str('  NEDOS', self.f_name, n=1), 5))
      input['nions']   = int(chk_len(get_str('  NIONS', self.f_name, n=1), 11))
      input['ispin']   = int(chk_len(get_str('  ISPIN', self.f_name, n=1), 2))
      input['lorbit']  = int(chk_len(get_str('  LORBIT', self.f_name, n=1), 2))
      input['lwave']   = chk_len(get_str('  LWAVE', self.f_name, n=1), 2)
      input['lcharg']  = chk_len(get_str('  LCHARG', self.f_name, n=1), 2)
      input['lvtot']   = chk_len(get_str('  LVTOT', self.f_name, n=1), 2)
      input['lsorbit'] = chk_len(get_str('  LSORBIT', self.f_name, n=1), 2)
      input['ncollin'] = chk_len(get_str('  LNONCOLLINEAR', self.f_name, n=1), 2)
      input['nkpts']   = int(chk_len(get_str('  NKPTS', self.f_name, n=1), 3))
      input['nbands']  = int(chk_len(get_str('  NBANDS=', self.f_name, n=1), 14))
      input['nelect']  = float(chk_len(get_str('  NELECT =', self.f_name, n=1), 2))
      input['atoms']   = self.get_atoms()
      input['nspec']   = get_str('ions per type =', self.f_name, n=1)[4:]
      for idx,num in enumerate(input['nspec']): input['nspec'][idx] = int(num)
      return input

   def get_energy(self):
      enes = {}
      enes['E_fermi']  = np.float_(chk_len(get_str('E-fermi', self.f_name, n=0), 2))
      enes['E_free']   = np.float_(chk_len(get_str('free  energy', self.f_name, n=0), 4))
      line = get_str('energy  without', self.f_name, n=0)
      enes['E_w/oS']   = np.float_(chk_len(line, 3))
      enes['E_sigma0'] = np.float_(chk_len(line, 6))
      return enes

   def get_ended(self):
      isend = False
      self.loop = int(os.popen('grep \'EDIFF is reached\' %s | wc -l'%(sefl.f_name)).read())
      if self.is_stat:
         if self.loop == 1: isend = True
      else:
         line = get_str('reached required accuracy', self.f_name, n=1)
         if line is not None: isend = True
      return isend

   def get_atoms(self):
      lines = get_str('VRHFIN', self.f_name, n='a')
      atoms = []
      for line in lines:
         atoms.append(line[1].replace(':','').replace('=',''))
      return atoms

   def get_lattice(self):
      """
      .get_lattice()
      ------------------------------------------------------------------
      return 2*3*3 optimized real and reciprocal lattice vectors
      ------------------------------------------------------------------
      get_lattice()[0] : 3*3 real lattice vector
      get_lattice()[1] : 3*3 reciprocal lattice vector
      """
      f = open(self.f_name)
      loc_lvec = None
      while (f.tell() <= self.eof):
         line = f.readline()
         if 'direct lattice vectors' in line: loc_lvec = f.tell()
      if loc_lvec is None: call_error('get_lattice - reading error', \
                                      'lattice information is missing')
      f.seek(loc_lvec)
      lvec = np.zeros([2,3,3],dtype=np.float_)
      for i in range(3):
         line = f.readline().split()
         for j in range(3):
            lvec[0][i][j] = np.float_(line[j])
            lvec[1][i][j] = np.float_(line[j+3])
      f.close()
      return lvec

   def get_magnetize(self):
      """
      .get_magnetize()
      ------------------------------------------------------------------
      return n*3 atomic magnetization if spin-polarized or non-collinear
      calculation is performed with LORBIT >= 10
      ------------------------------------------------------------------
      n : number of total atoms in POSCAR order
      - collinear calculation with ISPIN = 2 -
      get_magnetize()[N][0:1] : 0, no meaning
      get_magnetize()[N][2]   : + is spin up and - is spin down
      - non collinear calculation
      get_magnetize()[N][0:2] : magnetization vector
      """
      iscollin = True
      mag = np.zeros([self.input['nions'], 3],dtype=np.float_)
      if self.input['ncollin'] == 'T': iscollin = False
      if not (self.input['ispin'] == 2 or (not iscollin)): 
         call_error('input error - get_magnetize','spin polarization required')
      if self.input['lorbit'] < 10: 
         call_error('input error - get_magnetize','LORBIT >= 10 required')
      f = open(self.f_name)
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
      """
      .get_total_magnetize()
      ------------------------------------------------------------------
      return total magnetization (or vector) if spin-polarized or non-collinear
      calculation is performed
      ------------------------------------------------------------------
      collinear with ISPIN = 2 : total magnetization 
      non-collinear            : total magnetization vector
      """
      iscollin = True
      if self.input['ncollin'] == 'T': iscollin = False
      if not (self.input['ispin'] == 2 or (not iscollin)):
         call_error('input error - get_total_magnetize','spin polarization required')
      if not iscollin: mag = np.zeros([3], dtype=np.float_)
      f = open(self.f_name)
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
      """
      .get_force()
      ------------------------------------------------------------------
      return n*3 atomic force 
      """
      force = np.zeros([self.input['nions'], 3], dtype=np.float_)
      f = open(self.f_name)
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
   get_doscar(f_name, input)
   ------------------------------------------------------------------
   read VASP DOSCAR and return dataset.
   input : INCAR tags get from lp.get_outcar().input
   ------------------------------------------------------------------
   - output -
   .tdos()           : total DOS with E-EF
   .get_pdos()       : lm and site projected DOS (PDOS) if LORBIT is not 0
   .write_dos(n,d,h) : write dos datafile. requires 'f_name' f,
                       'dos data' d, and 'header' h
   .nedos    : number of grid
   .natoms   : number of atoms, required by PDOS
   .ef       : fermi level
   .header_T : header of total DOS
   .header_P : header of PDOS
   """
   def __init__(self, f_name, inp):
      self.f_name = f_name
      self.f = assert_file(self.f_name, 'get_doscar', isopen=True)
      self.natoms = int(chk_len(self.f.readline().split(), 0))
      for i in range(4): self.f.readline() # read buffer lines
      line = self.f.readline().split() # data line
      self.nedos = int(chk_len(line, 2)) - 1
      self.ef    = np.float_(chk_len(line, 3))
      if not assert_dat(self.nedos, inp['nedos']):
         call_error('input error - NEDOS','%i %i'%(inp['nedos'],self.nedos))
      if not assert_dat(self.natoms, inp['nions']):
         call_error('input error - NIONS','%i %i'%(inp['nions'],self.natoms))
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
      if self.iscollin and self.ispin == 1:
         self.header_P += header_tmp
      elif self.ispin == 2:
         for str in header_tmp:
            self.header_P.append(str+'_U')
            self.header_P.append(str+'_D')
      else:
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
      .get_pdos() 
      ------------------------------------------------------------------
      return PDOS with dimension of NATOMS*NEDOS*DATA
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
      .write_dos(f,d,h)
      ------------------------------------------------------------------
      write dos datafile with filename 'f' using data 'd' and 
      data description header 'h'
      """
      f = open(fname_out,'w')
      f.write('#%9s'%(header[0]))
      for str in header[1:]: f.write(' %10s'%(str))
      f.write('\n')
      for idos in range(self.nedos):
         for dat in dos[idos]:
            f.write('%10.4f '%(dat))
         f.write('\n')

class get_eigenval(object):
   """
   get_eigenval(f_name)
   ------------------------------------------------------------------
   read VASP EIGENVAL and return dataset.
   ------------------------------------------------------------------
   - output -
   .f_name   : filename of EIGENVAL
   .ispin    : ISPIN
   .nelect   : number of electrons
   .nkpts    : number of kpoints
   .nbands   : number of bands
   .kpts     : kpoint coordination in xyz (dim: .nkpts*3)
   .weight   : weight of kpoint (dim: .nkpts)
   .eigens   : eigenvalues (dim: .nkpts*.nbands*.ispin)
   """
   def __init__(self, f_name, inp):
      self.f_name = f_name
      f = assert_file(self.f_name, 'get_eigenval', isopen=True)

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

      if not assert_dat(inp['nelect'], self.nelect):
         call_error('input error - NELECT','%i %i'%(inp['nelect']), self.nelect)
      if not assert_dat(inp['ispin'], self.ispin):
         call_error('input error - ISPIN','%i %i'%(inp['ispin']), self.ispin)

# get eigenvalues
      for j in range(self.nkpts):
         f.readline() # empty line
         line = f.readline().split()
         self.weight[j] = np.float_(chk_len(line,3))
         for i in range(3): self.kpts[j][i] = np.float_(line[i])
         for i in range(self.nbands):
            line = f.readline().split()
            if len(line) < self.ispin+1:
               call_error('reading error - get_eigenval',self.f_name)
            for k in range(self.ispin): self.eigens[j][i][k] = np.float_(line[k+1])
      f.close()
   
#   def __isfile_eigenval(self):
      
class get_poscar(object):
   """
   get_poscar(f_name)
   ------------------------------------------------------------------
   read VASP POCSCAR 'f_name' and return dataset.
   ------------------------------------------------------------------
   - output -
   .f_name   : filename of POSCAR
   .des      : description of POSCAR at the first line
   .lvec     : lattice vector
   .ispbc    : periodic boundary condition, True
   .style    : POSCAR style, 4 or 5
   .natoms   : total number of atoms
   .atoms    : atom species (empty array returned if .style = 4)
   .nspec    : number of each atom species
   .isselc   : True if selective dynamics used
   .iscart   : True if coordinate is cartesian
   .coord    : atomic coordination in xyz (dim : .natoms*3)
   .vector   : vector in xyz (dim : .natoms*3)
   .selc     : if isselc is True, return selective dynamics flags
   """
   def __init__(self, f_name):
      self.f_name = f_name
      f = assert_file(self.f_name, 'get_poscar', isopen = True)
      self.des = f.readline()
# get lattice constant
      buf = f.readline().split()
      if not assert_num(buf[0]): call_error('reading error - get_poscar','lcont - '+buf[0])
      latcont = np.float_(buf[0])
# get lattice vector
      self.lvec = np.zeros((3,3), dtype=np.float_)
      for i in range(3):
         buf = f.readline().split()
         for j in range(3):
            if not assert_num(buf[j]): call_error('reading error - get_poscar','lvec - '+buf[j])
            self.lvec[i][j] = np.float_(buf[j])*latcont
      self.ispbc = True
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
      for i, nspec in enumerate(self.nspec):
         if not assert_num(nspec): call_error('reading error - get_poscar','nspec - '+nspec)
         self.nspec[i] = int(nspec)
         self.natoms  += self.nspec[i]
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
         if len(buf) < 3: call_error('reading error - get_poscar','coo - %i'%(i))
         for j in range(3):
            if not assert_num(buf[j]): call_error('reading error - get_poscar','coo - %i'%(i)+buf[j])
            self.coord[i][j] = np.float_(buf[j])
         if self.isselc:
            if len(buf) < 6: call_error('reading error - get_poscar','selc - %i'%(i))
            self.selc.append(buf[3:6])
# get vectors if exists
      f.readline()
      self.vector = np.zeros((self.natoms,3), dtype=np.float_)
      for i in range(self.natoms):
         buf = f.readline().split()
         if len(buf) == 0 or len(buf) != 3:
            self.vector = None
            return
         for j in range(3):
            if assert_num(buf[j]): 
               self.vector[i][j] = np.float_(buf[j])

class write_poscar(object):
   """
   write_poscar(f_name, dat, islong=True)
   ------------------------------------------------------------------
   get dataset 'dat' then write VAPS coordination POSCAR named 'filname'.
   ------------------------------------------------------------------
   - required input -
   .dat.des    : description of file at first line
   .dat.lvec   : lattice vector
   .dat.natoms : total number of atoms
   .dat.atoms  : atom species
   .dat.nspec  : number of each atom species
   .dat.isselc : whether write selective dynamics or not
   .dat.iscart : if dat.coord is fractional, False
   .dat.coord  : atomic coordination in xyz (dim : .natoms*3)
   .dat.vector : vector in xys (dim : .natoms*3) 
   .dat.selc   : selective dynamics flags
   """
   def __init__(self, f_name, dat, islong=True):
      if not (('POSCAR' in f_name) or ('CONTCAR' in f_name)):
         f_name = 'CONTCAR'+f_name
      f = open(f_name,'w')
      if dat.iscart: dat.coord = cart_to_dir(dat.lvec, dat.coord)
      if dat.isselc and (np.shape(dat.selc) != np.shape(dat.coord)):
         call_error('input error - write_poscar','dim - %s %s'%(\
			np.shape(dat.selc), np.shape(dat.coord),))
# write description
      f.write(dat.des)
# write lattice inform.
      f.write('1.000000000\n')
      for i in dat.lvec: 
         f.write(' ')
         if     islong: f.write(''.join(['%22.16f'%(j) for j in i])+'\n')
         if not islong: f.write(''.join(['%12.6f'%(j) for j in i])+'\n')
# write atom info 
      f.write(''.join(['%4s'%(i) for i in dat.atoms])+'\n')
      f.write(''.join(['%4i'%(i) for i in dat.nspec])+'\n')
      if dat.isselc: f.write('Selective Dynamics\n')
      f.write('Direct\n')
# write coordinate
      for idx,coo in enumerate(dat.coord):
         if     islong: f.write(''.join(['%20.16f'%(j) for j in coo]))
         if not islong: f.write(''.join(['%10.6f'%(j) for j in coo]))
         if dat.isselc: f.write(''.join(['%3s'%(j) for j in dat.selc[idx]]))
         f.write('\n')
      if dat.vector is None:
         f.close()
         return
# write vector
      f.write('\n')
      for vec in dat.vector:
         if     islong: f.write(''.join(['%20.16f'%(j) for j in vec])+'\n')
         if not islong: f.write(''.join(['%10.6f'%(j) for j in vec])+'\n')
      f.close()

class get_xsf(object):
   """
   get_xsf(f_name, dat)
   ------------------------------------------------------------------
   read XCrysDen XSF 'f_name' and return dataset 'dat'
   ------------------------------------------------------------------
   - output -
   .f_name   : filename of XSF
   .des      : description of XSF at the first line
   .lvec     : lattice vector
   .iscart   : cartesian coordinate, True
   .ispbc    : periodic boundary condition
   .natoms   : total number of atoms
   .atoms    : atom species
   .nspec    : number of each atom species
   .coord    : atomic coordination in xyz (dim : .natoms*3)
   .vector   : vector, drawn in VESTA (dim : .natoms*3)
   """
   def __init__(self, f_name):
      self.f_name = f_name
      f = assert_file(self.f_name, 'get_xsf', isopen=True)
      self.des = f.readline()
      self.iscart = True
      cnt = 0
# get lattice inform (pbc or not)
      while True:
         buf  = f.readline().split()
         cnt += 1
         if buf[0] == 'CRYSTAL': 
            self.ispbc = True
            break
         elif buf[0] == 'ATOMS':
            self.ispbc = False
            break
         if cnt > 100: call_error('reading error - get_xsf','STRUCTURE TYPE not found')
      self.lvec  = np.zeros((3,3), dtype=np.float_)
      cnt = 0
# if pbc, get lattice
      if self.ispbc:
         while True:
            buf  = f.readline().split()
            cnt += 1
            if (buf[0] == 'PRIMVEC') or (buf[0] == 'CONVVEC'): break
            if cnt > 100: call_error('reading error - get_xsf','LATTICE TYPE not found')
         for i in range(3):
            buf = f.readline().split()
            for j in range(3):
               if not assert_num(buf[j]): call_error('reading error - get_xsf','lvec - '+buf[j])
               self.lvec[i][j] = np.float_(buf[j])
         cnt = 0
# if pbc, read up to coord info
         while True:
            buf  = f.readline().split()
            cnt += 1
            if buf[0] == 'PRIMCOORD': break
            if cnt > 100: call_error('reading error - get_xsf','PRIMCOORD not found')
         buf = f.readline().split()
         if not assert_num(buf[0]): call_error('reading error - get_xsf','NATOMS - '+buf[0])
         self.natoms = int(buf[0])
      else:
# if not pbc, read up to coord info
         f.readline()
         bufloc = f.tell()
         self.natoms = 0
         while True:
            buf = f.readline().split()
            if len(buf) < 4 or buf[0] == '#': break
            self.natoms += 1
         self.seek(bufloc)
# data init
      self.coord   = np.zeros((self.natoms,3), dtype=np.float_)
      self.vector  = np.zeros((self.natoms,3), dtype=np.float_)
      coord        = np.zeros((self.natoms,3), dtype=np.float_)
      vector       = np.zeros((self.natoms,3), dtype=np.float_)
      self.nspec   = [0]
      self.atoms   = []
      self.atoms_l = ['' for j in range(self.natoms)]
      getvec = False
# read coord and vec
      for i in range(self.natoms):
         buf = f.readline().split()
         if len(buf) < 4: call_error('reading error - get_xsf','coord - %i %i'%(i, len(buf)))
         self.atoms_l[i] = buf[0]
         for j in range(3): coord[i][j] = np.float_(buf[j+1])
         if len(buf) > 6: getvec = True
         if getvec and len(buf) < 7: call_error('reading error - get_xsf','vector - %i %i'%(i, len(buf)))
         if getvec: 
            for j in range(3): vector[i][j] = np.float_(buf[j+4])
      self.atoms.append(self.atoms_l[0])
      for atom_ref in self.atoms_l:
         isdup = False
         for atom in self.atoms:
            if atom == atom_ref: isdup = True
            if isdup: break
         if not isdup:
            self.atoms.append(atom_ref)
            self.nspec.append(0)
      cnt = 0
# arrange coord
      if not getvec: self.vector = None
      for idx_ref,atom_ref in enumerate(self.atoms):
         for idx,atom in enumerate(self.atoms_l):
            if atom == atom_ref:
               cnt += 1
               self.nspec[idx_ref] += 1
               self.coord[cnt] = coord[idx]
               if getvec: self.vector[cnt] = vector[idx]

class write_xsf(object):
   """
   write_xsf(f_name, dat, iimg=1)
   ------------------------------------------------------------------
   get dataset 'dat' then write XCrysDen coordination file (.xsf) 
   named 'f_name'.
   ------------------------------------------------------------------
   - required input -
   .dat.des    : description of XSF at first line
   .dat.lvec   : lattice vector
   .dat.ispbc  : periodic boundary or not
   .dat.natoms : total number of atoms
   .dat.atoms  : atom species
   .dat.nspec  : number of each atom species
   .dat.iscart : if dat.coord is fractional, False
   .dat.coord  : atomic coordination in xyz (dim : .natoms*3)
   .dat.vector : vector, drawn in VESTA (dim : .natoms*3)
   """
   def __init__(self, f_name, dat, iimg=1):
      if not '.xsf' in f_name: f_name = f_name+'.xsf'
      f = open(f_name,'w')
      atoms_l = [at for ii,at in enumerate(dat.atoms) for i in range(dat.nspec[ii])]
      if not dat.iscart: dat.coord = dir_to_cart(dat.lvec, dat.coord)
      f.write('# %s\n'%(dat.des))
      if (dat.ispbc):
         f.write('CRYSTAL\n')
         f.write('PRIMVEC\n')
         for i in dat.lvec: f.write(' '+''.join(['%22.16f'%(j) for j in i])+'\n')
         f.write('PRIMCOORD\n')
         f.write('%i %i\n'%(dat.natoms, iimg))
      else:
         f.write('ATOMS\n')
      for i in range(dat.natoms):
         f.write('%2s'%(atoms_l[i]))
         f.write(''.join([' %23.16f'%(j) for j in dat.coord[i][:]]))
         if dat.vector is not None:
            f.write(''.join([' %13.6f'%(j) for j in dat.vector[i][:]]))
         f.write('\n')
      f.close()

