def read_RT_Hmean(dir,iter):

  tmp = np.fromfile(dir+'RT_mean1D.'+inttostring(iter,ts_size=6),dtype=np.float32)

  Nbands = tmp[0].astype(int); 
  Nvar = tmp[1].astype(int);
  Nz = tmp[2].astype(int);
  time = tmp[3];

  return (tmp[8:].reshape([Nbands,Nvar,Nz]).swapaxes(0,1)).swapaxes(1,2),Nbands,Nvar,Nz,time

def read_Iout(dir,iter):

  tmp = np.fromfile(dir+'I_out.'+inttostring(iter,ts_size=6),dtype=np.float32)

  size = tmp[1:3].astype(int)
  time = tmp[3]

  return tmp[4:].reshape([size[1],size[0]]).swapaxes(0,1),size,time

def read_header(dir,iter):

    header = np.loadtxt(dir+'Header.'+inttostring(iter,ts_size=6))

    return header

def read_var_3d(dir,var,iter,layout=None):

  h = read_header(dir,iter)

  size = h[0:3].astype(int)
  dx   = h[3:6]
  time = h[6]

  tmp = np.fromfile(dir+var+'.'+ inttostring(iter,ts_size=6),dtype=np.float32)
  tmp = tmp.reshape([size[2],size[1],size[0]])

  if layout != None :
      tmp = tmp.transpose(layout)

  return tmp,dx,size,time


def read_dem(path,dir,iter,max_bins=None):

  tmp = np.fromfile(path+'corona_emission_adj_dem_'+dir+'.'+ inttostring(iter,ts_size=6),dtype=np.float32)

  bins   = tmp[0].astype(int)
  size   = tmp[1:3].astype(int)
  time   = tmp[3]
  lgTmin = tmp[4]
  dellgT = tmp[5]

  dem = tmp[6:].reshape([bins,size[1],size[0]]).transpose(2,1,0)

  taxis = lgTmin+dellgT*np.arange(0,bins+1)

  X_H = 0.7
  dem = dem*X_H*0.5*(1+X_H)*3.6e19

  if max_bins != None:
    if bins > max_bins :
      dem = dem[:,:,0:max_bins]
    else :
      tmp=dem
      dem=np.zeros([size[0],size[1],max_bins])
      dem[:,:,0:bins]=tmp

    taxis = lgTmin+dellgT*np.arange(0,max_bins+1)

  return dem,taxis,time

def read_vlos(path,dir,iter,max_bins=None):

  tmp = np.fromfile(path+'corona_emission_adj_vlos_'+dir+'.'+ inttostring(iter,ts_size=6),dtype=np.float32)

  bins   = tmp[0].astype(int)
  size   = tmp[1:3].astype(int)
  time   = tmp[3]
  lgTmin = tmp[4]
  dellgT = tmp[5]

  vlos = tmp[6:].reshape([bins,size[1],size[0]]).transpose(2,1,0)

  taxis = lgTmin+dellgT*np.arange(0,bins+1)

  if max_bins != None:
    if bins > max_bins :
      vlos = vlos[:,:,0:max_bins]
    else :
      tmp=vlos
      vlos=np.zeros([size[0],size[1],max_bins])
      vlos[:,:,0:bins]=tmp

    taxis = lgTmin+dellgT*np.arange(0,max_bins+1)

  return vlos,taxis,time


def read_vrms(path,dir,iter,max_bins=None):

  tmp = np.fromfile(path+'corona_emission_adj_vrms_'+dir+'.'+ inttostring(iter,ts_size=6),dtype=np.float32)

  bins   = tmp[0].astype(int)
  size   = tmp[1:3].astype(int)
  time   = tmp[3]
  lgTmin = tmp[4]
  dellgT = tmp[5]

  vlos = tmp[6:].reshape([bins,size[1],size[0]]).transpose(2,1,0)

  taxis = lgTmin+dellgT*np.arange(0,bins+1)

  if max_bins != None:
    if bins > max_bins :
      vlos = vlos[:,:,0:max_bins]
    else :
      tmp=vlos
      vlos=np.zeros([size[0],size[1],max_bins])
      vlos[:,:,0:bins]=tmp

    taxis = lgTmin+dellgT*np.arange(0,max_bins+1)

  return vlos,taxis,time

def read_fil(path,dir,iter,max_bins=None):

  tmp = np.fromfile(path+'corona_emission_adj_fil_'+dir+'.'+ inttostring(iter,ts_size=6),dtype=np.float32)

  bins   = tmp[0].astype(int)
  size   = tmp[1:3].astype(int)
  time   = tmp[3]
  lgTmin = tmp[4]
  dellgT = tmp[5]

  vlos = tmp[6:].reshape([bins,size[1],size[0]]).transpose(2,1,0)

  taxis = lgTmin+dellgT*np.arange(0,bins+1)

  if max_bins != None:
    if bins > max_bins :
      vlos = vlos[:,:,0:max_bins]
    else :
      tmp=vlos
      vlos=np.zeros([size[0],size[1],max_bins])
      vlos[:,:,0:bins]=tmp

    taxis = lgTmin+dellgT*np.arange(0,max_bins+1)

  return vlos,taxis,time

def read_slice(dir,var,depth,iter):

  tmp = np.fromfile(dir+var+'_slice_'+depth+'.'+inttostring(iter,ts_size=6),dtype=np.float32)

  nslices = tmp[0].astype(int)
  size = tmp[1:3].astype(int)
  time = tmp[3]

  return tmp[4:].reshape([nslices,size[1],size[0]]).swapaxes(1,2),nslices,size,time

def inttostring(ii,ts_size=6):

  str_num = str(ii)

  for bb in range(len(str_num),ts_size,1):
    str_num = '0'+str_num
  
  return str_num

def MURaM_output(arr,dir,output,iter = '000000',prim=True, precision = 'single'):

  if output is 'vx':
      filename = 'result_2.'
  if output is 'vy':
      filename = 'result_3.'
  if output is 'vz':
      filename = 'result_1.'
  if output is 'bx':
      filename = 'result_6.'
      arr = arr/np.sqrt(4.0*np.pi)
  if output is 'by':
      filename = 'result_7.'
      arr = arr/np.sqrt(4.0*np.pi)
  if output is 'bz':
      filename = 'result_5.'
      arr = arr/np.sqrt(4.0*np.pi)
  if output is 'rho':
      filename = 'result_0.'
  if output is 'eps':
      filename = 'result_4.'
  if output is 'sflx':
      filename = 'result_8.'

  if prim: filename=filename.replace('_','_prim_')

  filename = dir + filename + iter

  ## Revert to z,x,y ordering

  if (precision is 'double'):
    arr.astype(np.double)
  if (precision is 'single'):
    arr.astype(np.single)

  arr.swapaxes(0,1).ravel().tofile(filename)

class MURaM_snap:
  
  def __init__(self,dir,filename="parameters.dat"):
    ## Get simulation parameters from a particular parameters.dat (or similar) file
    ## Define toloads
    self.dir = dir
    
  def load(self,iter,primative=True,tooload=['rho','vx','vy','vz','bx','by','bz']):

    def loadrho(iter,primative):
      if primative:
        str = '3D/result_prim_0.'
      else:
        str = '3D/result_0.'

      self.rho = np.zeros([self.nx,self.ny,self.nz],dtype = np.float32)
      self.rho[:,:,:] = np.fromfile(self.dir + str+iter,dtype = np.float32).reshape([self.ny,self.nx,self.nz]).swapaxes(0,1)
   
    def loadvx(iter,primative):

      if primative:
        str = '3D/result_prim_2.'
      else:
        str = '3D/result_2.'

      self.vx = np.zeros([self.nx,self.ny,self.nz],dtype = np.float32)
      self.vx[:,:,:] =  np.fromfile(self.dir + str+iter,dtype = np.float32).reshape([self.ny,self.nx,self.nz]).swapaxes(0,1)
      if not primative:
        self.vx = self.vx/self.rho

    def loadvy(iter,primative):

      if primative:
        str = '3D/result_prim_3.'
      else:
        str = '3D/result_3.'

      self.vy = np.zeros([self.nx,self.ny,self.nz],dtype = np.float32)
      self.vy[:,:,:] =  np.fromfile(self.dir + str+iter,dtype = np.float32).reshape([self.ny,self.nx,self.nz]).swapaxes(0,1)
      if not primative:
        self.vy = self.vy/self.rho

    def loadvz(iter,primative):

      if primative:
        str = '3D/result_prim_1.'
      else:
        str = '3D/result_1.'

      self.vz = np.zeros([self.nx,self.ny,self.nz],dtype = np.float32)
      self.vz[:,:,:] =  np.fromfile(self.dir + str+iter,dtype = np.float32).reshape([self.ny,self.nx,self.nz]).swapaxes(0,1)
      if not primative:
        self.vz = self.vz/self.rho

    def loadeps(iter,primative):
      if primative:
        str = '3D/result_prim_4.'
      else:
        str = '3D/result_4.'

      self.eps = np.zeros([self.nx,self.ny,self.nz],dtype = np.float32)
      self.eps[:,:,:] =  np.fromfile(self.dir + str+iter,dtype = np.float32).reshape([self.ny,self.nx,self.nz]).swapaxes(0,1) 
      self.eps = self.eps/self.rho

    def loadbx(iter,primative):

      if primative:
        str = '3D/result_prim_6.'
      else:
        str = '3D/result_6.'

      self.bx = np.zeros([self.nx,self.ny,self.nz],dtype = np.float32)
      self.bx[:,:,:] =  np.fromfile(self.dir + str+iter,dtype = np.float32).reshape([self.ny,self.nx,self.nz]).swapaxes(0,1)*np.sqrt(4.0*np.pi)

    def loadby(iter,primative):
      
      if primative:
        str = '3D/result_prim_7.'
      else:
        str = '3D/result_7.'

      self.by = np.zeros([self.nx,self.ny,self.nz],dtype = np.float32)
      self.by[:,:,:] =  np.fromfile(self.dir + str+iter,dtype = np.float32).reshape([self.ny,self.nx,self.nz]).swapaxes(0,1)*np.sqrt(4.0*np.pi)

    def loadbz(iter,primative):
      
      if primative:
        str = '3D/result_prim_5.'
      else:
        str = '3D/result_5.'

      self.bz = np.zeros([self.nx,self.ny,self.nz],dtype = np.float32)
      self.bz[:,:,:] =  np.fromfile(self.dir + str+iter,dtype = np.float32).reshape([self.ny,self.nx,self.nz]).swapaxes(0,1)*np.sqrt(4.0*np.pi)

    def loadtem(iter,primative):

      str = '3D/eosT.' 
      self.tem = np.zeros([self.nx,self.ny,self.nz],dtype = np.float32)
      self.tem[:,:,:] = np.fromfile(self.dir + str+iter,dtype = np.float32).reshape([self.ny,self.nx,self.nz]).swapaxes(0,1)

    def loadpre(iter,primative):

      str = '3D/eosP.'
      self.pre = np.zeros([self.nx,self.ny,self.nz],dtype = np.float32)
      self.pre[:,:,:] = np.fromfile(self.dir + str+iter,dtype = np.float32).reshape([self.ny,self.nx,self.nz]).swapaxes(0,1)

    def loadQtot(iter,primative):
      
      str = '3D/Qtot.'
      self.Qtot = np.zeros([self.nx,self.ny,self.nz],dtype = np.float32)
      self.Qtot[:,:,:] = np.fromfile(self.dir + str+iter,dtype = np.float32).reshape([self.ny,self.nx,self.nz]).swapaxes(0,1)
    
    def loadJtot(iter,primative):
      
      str = '3D/Jtot.'
      self.Jtot = np.zeros([self.nx,self.ny,self.nz],dtype = np.float32)
      self.Jtot[:,:,:] = np.fromfile(self.dir + str+iter,dtype = np.float32).reshape([self.ny,self.nx,self.nz]).swapaxes(0,1)

    def loadStot(iter,primative):
      
      str = '3D/Stot.'
      self.Stot = np.zeros([self.nx,self.ny,self.nz],dtype = np.float32)
      self.Stot[:,:,:] = np.fromfile(self.dir + str+iter,dtype = np.float32).reshape([self.ny,self.nx,self.nz]).swapaxes(0,1)

    def loadtau(iter,primative):
      
      str = '3D/tau.'
      self.tau = np.zeros([self.nx,self.ny,self.nz],dtype = np.float32)
      self.tau[:,:,:] = np.fromfile(self.dir + str+iter,dtype = np.float32).reshape([self.ny,self.nx,self.nz]).swapaxes(0,1)

    def loadQxH(iter,primative):
      
      str = '3D/QxH.'
      self.QxH = np.zeros([self.nx,self.ny,self.nz],dtype = np.float32)
      self.QxH[:,:,:] = np.fromfile(self.dir + str+iter,dtype = np.float32).reshape([self.ny,self.nx,self.nz]).swapaxes(0,1)
    
    def loadne(iter,primative):
      
      str = '3D/eosne.'
      self.ne = np.zeros([self.nx,self.ny,self.nz],dtype = np.float32)
      self.ne[:,:,:] = np.fromfile(self.dir + str+iter,dtype = np.float32).reshape([self.ny,self.nx,self.nz]).swapaxes(0,1)

    def loadrhoi(iter,primative):
      
      str = '3D/eosrhoi.'
      self.rhoi = np.zeros([self.nx,self.ny,self.nz],dtype = np.float32)
      self.rhoi[:,:,:] = np.fromfile(self.dir + str+iter,dtype = np.float32).reshape([self.ny,self.nx,self.nz]).swapaxes(0,1)

    def loadamb(iter,primative):
      
      str = '3D/eosamb.'
      self.amb = np.zeros([self.nx,self.ny,self.nz],dtype = np.float32)
      self.amb[:,:,:] = np.fromfile(self.dir + str+iter,dtype = np.float32).reshape([self.ny,self.nx,self.nz]).swapaxes(0,1)

    def loadQxCa(iter,primative):
      
      str = '3D/QxCa.'
      self.QxCa = np.zeros([self.nx,self.ny,self.nz],dtype = np.float32)
      self.QxCa[:,:,:] = np.fromfile(self.dir + str+iter,dtype = np.float32).reshape([self.ny,self.nx,self.nz]).swapaxes(0,1)

    def loadQxMg(iter,primative):
      
      str = '3D/QxMg.'
      self.QxMg = np.zeros([self.nx,self.ny,self.nz],dtype = np.float32)
      self.QxMg[:,:,:] = np.fromfile(self.dir + str+iter,dtype = np.float32).reshape([self.ny,self.nx,self.nz]).swapaxes(0,1)
     
    def loadQxCor(iter,primative):
      
      str = '3D/QxCor.'
      self.QxCor = np.zeros([self.nx,self.ny,self.nz],dtype = np.float32)
      self.QxCor[:,:,:] = np.fromfile(self.dir + str+iter,dtype = np.float32).reshape([self.ny,self.nx,self.nz]).swapaxes(0,1)
    
    def loadQxChr(iter,primative):
      
      str = '3D/QxChr.'
      self.QxChr = np.zeros([self.nx,self.ny,self.nz],dtype = np.float32)
      self.QxChr[:,:,:] = np.fromfile(self.dir + str+iter,dtype = np.float32).reshape([self.ny,self.nx,self.nz]).swapaxes(0,1)

    def loadQres(iter,primative):
      
      str = '3D/Qres.'
      self.Qres = np.zeros([self.nx,self.ny,self.nz],dtype = np.float32)
      self.Qres[:,:,:] = np.fromfile(self.dir + str+iter,dtype = np.float32).reshape([self.ny,self.nx,self.nz]).swapaxes(0,1)
     
    def loadQvis(iter,primative):
      
      str = '3D/Qvis.'
      self.Qvis = np.zeros([self.nx,self.ny,self.nz],dtype = np.float32)
      self.Qvis[:,:,:] = np.fromfile(self.dir + str+iter,dtype = np.float32).reshape([self.ny,self.nx,self.nz]).swapaxes(0,1)


    def loadsflx(iter,primative):
      
      if primative:
        str = '3D/result_prim_8.'
      else:
        str = '3D/result_8.'
      self.sflx = np.zeros([self.nx,self.ny,self.nz],dtype = np.float32)
      self.sflx[:,:,:] = np.fromfile(self.dir + str+iter,dtype = np.float32).reshape([self.ny,self.nx,self.nz]).swapaxes(0,1)

    options = {'rho' : loadrho,
               'bx' : loadbx,
               'by' : loadby,
               'bz' : loadbz,
               'vx' : loadvx,
               'vy' : loadvy,
               'vz' : loadvz,
               'eps' : loadeps,
               'tem' : loadtem,
               'pre' : loadpre,
               'Qtot' : loadQtot,
               'Jtot' : loadJtot,
               'Stot' : loadStot,
               'tau' : loadtau,
               'QxH' : loadQxH,
               'QxCa' : loadQxCa,
               'QxMg' : loadQxMg,
               'QxCor' : loadQxCor,
               'QxChr' : loadQxChr,
               'Qres' : loadQres,
               'Qvis' : loadQvis,
               'ne' : loadne,
               'amb' : loadamb,
               'rhoi': loadrhoi,
               'sflx' : loadsflx,
               }

    self.iter = np.int(iter)

    h = np.loadtxt(self.dir+'3D/Header.'+inttostring(iter,ts_size=6))
    
    self.nz = np.int(h[0])
    self.nx = np.int(h[1])
    self.ny = np.int(h[2])
    
    self.dz = h[3]
    self.dx = h[4]
    self.dy = h[5]
    
    self.gxmin = [0.0,0.0,0.0]
    self.gxmax = [self.dx*self.nx,self.dy*self.ny,self.dz*self.nz]
    
    self.time = h[6]
    self.dt = h[7]
    self.vamax = h[8]

    if self.dx == 1:
        self.NDIM = 1
    elif self.dy == 1:
        self.NDIM = 2
    else:
        self.NDIM = 3
        
    if self.NDIM > 1:
      self.xlength = self.nz*self.dz
    else:
      self.xlength = 0.0
      self.dx = 0.0


    if self.NDIM == 3:
      self.ylength = self.gxmax[1]-self.gxmin[1]
      self.dy = self.ylength/(self.ny)
    else:
      self.ylength = 0.0
      self.dy = 0.0

    self.xax = self.gxmin[0] + self.dx*np.arange(self.nx)
    self.yax = self.gxmin[1] + self.dy*np.arange(self.ny)
    self.zax = self.gxmin[2] + self.dz*np.arange(self.nz)
    
    # If not primative, and vx,vy,vz are needed we must include rho in tooload
    if primative and bool(set(tooload) & set(['vx','vy','vz'])):
      tooload.insert(0,'rho')
    
    # Make sure rho is loaded first so that vx, vy, vz can be calculated
    for pos in range(np.size(tooload)):
      if tooload[pos] is 'rho':
        tooload.insert(0, tooload.pop(pos))

    # Now load everything in too load
    for var in tooload:
      options[var](iter,primative)

  def zax_tau_correct(self):

    tau_av = np.mean(self.tau,axis=(0,1))
    for i in range(self.nz):
      if (tau_av[i] <= 1):
        z_tau_cor = self.zax[i-1] + (self.zax[i]-self.zax[i-1])/(tau_av[i] - tau_av[i-1])*(1.0-tau_av[i-1])
        print(self.zax[i],tau_av[i],z_tau_cor)
        self.zax = self.zax-z_tau_cor
        self.gxmin[2]= self.gxmin[2]-z_tau_cor
        self.gxmax[2]= self.gxmax[2]-z_tau_cor
        break

    return z_tau_cor

  def poynting(self):

    ex = self.vy*self.bz - self.vz*self.by
    ey = - self.vx*self.bz + self.vz*self.bx
    ez = self.vx*self.by - self.vy*self.bx

    fmag=np.zeros([3,self.nx,self.ny,self.nz],dtype = np.float64)

    fmag[0,:,:,:]= ( self.by*ez-self.bz*ey)
    fmag[1,:,:,:]= (-self.bx*ez+self.bz*ex)
    fmag[2,:,:,:]= ( self.bx*ey-self.by*ex)

    fmag = fmag/4.0/np.pi

    return fmag

  def calc_beta(self):
      beta = (8.0*np.pi*(self.pe))/np.clip((self.bx)**2 + (self.by)**2+(self.bz)**2,1.0e-20,None)
      return beta

  def calc_VA(self):
    va = np.sqrt(((self.bx)**2 + (self.by)**2+(self.bz)**2)/(4.0*np.pi*self.rho))

    return va

  def rotor(self):
    self.rotor = deriv_3d_O4(self.vx,1,self.dy,periodic=True) - deriv_3d_O4(self.vy,0,self.dx,periodic=True)
    return self.rotor

  def calc_jotaperp(self):

    jota_perp = np.zeros([3,self.nx,self.ny,self.nz],dtype=np.float64)
    jota = np.zeros([3,self.nx,self.ny,self.nz],dtype=np.float64)

    bi = np.sqrt(self.bx*self.bx + self.by*self.by + self.bz*self.bz)

    dbxdy = deriv_3d_O4(self.bx,1,delta=self.dy,periodic=True)
    dbzdy = deriv_3d_O4(self.bz,1,delta=self.dy,periodic=True)
    dbydx = deriv_3d_O4(self.by,0,delta=self.dx,periodic=True)
    dbzdx = deriv_3d_O4(self.bz,0,delta=self.dx,periodic=True)
    dbxdz = deriv_3d_O4(self.bx,2,delta=self.dz,periodic=False)
    dbydz = deriv_3d_O4(self.by,2,delta=self.dz,periodic=False)

    jota[0,:,:,:] =  (dbzdy-dbydz)
    jota[1,:,:,:] = (-dbzdx+dbxdz)
    jota[2,:,:,:] =  (dbydx-dbxdy)

    bj = self.bx*jota[0,:,:,:] + self.by*jota[1,:,:,:] + self.bz*jota[2,:,:,:]

    jota_perp[0,:,:,:] = jota[0,:,:,:] - bj*self.bx/bi
    jota_perp[1,:,:,:] = jota[1,:,:,:] - bj*self.by/bi
    jota_perp[2,:,:,:] = jota[2,:,:,:] - bj*self.bz/bi

    return jota_perp, jota

  def calc_KE(self):
    ekin=(self.rho)*(self.vx**2 + self.vy**2 + self.vz**2)/2.0

    return ekin

  def calc_bE(self):
    emag=((self.bx)**2 + (self.by)**2 + (self.bz)**2)/8.0/np.pi

    return emag

def bilinear_interpolate(im, x, y):
    x = np.asarray(x, dtype = 'float64')
    y = np.asarray(y, dtype = 'float64')

    x0 = np.floor(x).astype(int)
    x1 = x0 + 1
    y0 = np.floor(y).astype(int)
    y1 = y0 + 1

    x0 = np.clip(x0, 0, im.shape[0]-1);
    x1 = np.clip(x1, 0, im.shape[0]-1);
    y0 = np.clip(y0, 0, im.shape[1]-1);
    y1 = np.clip(y1, 0, im.shape[1]-1);

    Ia = im[x0,y0] # y0, x0 ]
    Ib = im[x0,y1] #y1, x0 ]
    Ic = im[x1,y0] #y0, x1 ]
    Id = im[x1,y1] #y1, x1 ]

    wa = (x1-x) * (y1-y)
    wb = (x1-x) * (y-y0)
    wc = (x-x0) * (y1-y)
    wd = (x-x0) * (y-y0)

    return wa*Ia + wb*Ib + wc*Ic + wd*Id

def deriv_3d_O4(arr,dir,delta=1.0,periodic=True):
    ## Do a derivative in 3 dim array
    ## Dir = array indices
    ## delta = grid spacing

    dim = arr.ndim
    if dim > 3:
        print("Currently coded for only a 3-dim array")
        print(arr.ndim)
        return None
    
    weights=[-1.0/12.0,8.0/12.0,-8.0/12.0,1.0/12.0]
    
    derr = np.zeros(arr.shape,dtype = np.float64)
    derr += weights[0]*np.roll(arr,2,dir)
    derr += weights[1]*np.roll(arr,1,dir)
    derr += weights[2]*np.roll(arr,-1,dir)
    derr += weights[3]*np.roll(arr,-2,dir)

    if not periodic:
        if dir == 0:
            derr[0,:,:] = arr[1,:,:] - arr[0,:,:]
            derr[1,:,:] = 0.5*(arr[2,:,:] - arr[2,:,:])
            derr[-2,:,:] = 0.5*(arr[-1,:,:] - arr[-3,:,:])
            derr[-1,:,:] = arr[-1,:,:] - arr[-2,:,:]
        if dir == 1:
            derr[:,0,:] = arr[:,1,:] - arr[:,0,:]
            derr[:,1,:] = 0.5*(arr[:,2,:] - arr[:,0,:])
            derr[:,-2,:] = 0.5*(arr[:,-1,:] - arr[:,-3,:])
            derr[:,-1,:] = arr[:,-1,:] - arr[:,-2,:]
        if dir == 2:
            derr[:,:,0] = arr[:,:,1] - arr[:,:,0]
            derr[:,:,1] = 0.5*(arr[:,:,2] - arr[:,:,0])
            derr[:,:,-2] = 0.5*(arr[:,:,-1] - arr[:,:,-3])
            derr[:,:,-1] = arr[:,:,-1] - arr[:,:,-2]
   
    derr /= delta

    return derr

def congrid(a, newdims, method='linear', centre=False, minusone=False):
    '''Arbitrary resampling of source array to new dimension sizes.
    Currently only supports maintaining the same number of dimensions.
    To use 1-D arrays, first promote them to shape (x,1).
    
    Uses the same parameters and creates the same co-ordinate lookup points
    as IDL''s congrid routine, which apparently originally came from a VAX/VMS
    routine of the same name.

    method:
    neighbour - closest value from original data
    nearest and linear - uses n x 1-D interpolations using
                         scipy.interpolate.interp1d
    (see Numerical Recipes for validity of use of n 1-D interpolations)
    spline - uses ndimage.map_coordinates

    centre:
    True - interpolation points are at the centres of the bins
    False - points are at the front edge of the bin

    minusone:
    For example- inarray.shape = (i,j) & new dimensions = (x,y)
    False - inarray is resampled by factors of (i/x) * (j/y)
    True - inarray is resampled by(i-1)/(x-1) * (j-1)/(y-1)
    This prevents extrapolation one element beyond bounds of input array.
    '''
    if not a.dtype in [np.float64, np.float32]:
        a = np.cast[float](a)

    m1 = np.cast[int](minusone)
    ofs = np.cast[int](centre) * 0.5
    old = np.array( a.shape )
    ndims = len( a.shape )
    if len( newdims ) != ndims:
        print("[congrid] dimensions error. " \
              "This routine currently only support " \
              "rebinning to the same number of dimensions.")
        return None
    newdims = np.asarray( newdims, dtype=float )
    dimlist = []

    if method == 'neighbour':
        for i in range( ndims ):
            base = np.indices(newdims)[i]
            dimlist.append( (old[i] - m1) / (newdims[i] - m1) \
                            * (base + ofs) - ofs )
        cd = np.array( dimlist ).round().astype(int)
        newa = a[list( cd )]
        return newa

    elif method in ['nearest','linear']:
        # calculate new dims
        for i in range( ndims ):
            base = np.arange( newdims[i] )
            dimlist.append( (old[i] - m1) / (newdims[i] - m1) \
                            * (base + ofs) - ofs )
        # specify old dims
        olddims = [np.arange(i, dtype = np.float) for i in list( a.shape )]

        # first interpolation - for ndims = any
        mint = scipy.interpolate.interp1d( olddims[-1], a, kind=method )
        newa = mint( dimlist[-1] )

        trorder = [ndims - 1] + list(range( ndims - 1))
        for i in range( ndims - 2, -1, -1 ):
            newa = newa.transpose( trorder )

            mint = scipy.interpolate.interp1d( olddims[i], newa, kind=method )
            newa = mint( dimlist[i] )

        if ndims > 1:
            # need one more transpose to return to original dimensions
            newa = newa.transpose( trorder )

        return newa
    elif method in ['spline']:
        oslices = [ slice(0,j) for j in old ]
        oldcoords = np.ogrid[oslices]
        nslices = [ slice(0,j) for j in list(newdims) ]
        newcoords = np.mgrid[nslices]

        newcoords_dims = list(range(n.rank(newcoords)))
        #make first index last
        newcoords_dims.append(newcoords_dims.pop(0))
        newcoords_tr = newcoords.transpose(newcoords_dims)
        # makes a view that affects newcoords

        newcoords_tr += ofs

        deltas = (np.asarray(old) - m1) / (newdims - m1)
        newcoords_tr *= deltas

        newcoords_tr -= ofs

        newa = scipy.ndimage.map_coordinates(a, newcoords)
        return newa
    else:
        print("Congrid error: Unrecognized interpolation type.\n", \
              "Currently only \'neighbour\', \'nearest\',\'linear\',", \
              "and \'spline\' are supported.")
        return None


import numpy as np
import struct
import scipy.interpolate
import scipy.ndimage

print ("Importing read_muram.py")
