import numpy as np

def inttostring(ii,ts_size=6):

  str_num = str(ii)

  for bb in range(len(str_num),ts_size,1):
    str_num = '0'+str_num
  
  return str_num


def read_Iout(dir,iter):

  tmp = np.fromfile(dir+'I_out.'+inttostring(iter,ts_size=6),dtype=np.float32)

  size = tmp[1:3].astype(int)
  time = tmp[3]

  return tmp[4:].reshape([size[1],size[0]]).swapaxes(0,1),size,time


def read_slice(dir,var,depth,iter):

  tmp = np.fromfile(dir+var+'_slice_'+depth+'.'+inttostring(iter,ts_size=6),dtype=np.float32)

  nslices = tmp[0].astype(int)
  size = tmp[1:3].astype(int)
  time = tmp[3]

  return tmp[4:].reshape([nslices,size[1],size[0]]).swapaxes(1,2),nslices,size,time


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

