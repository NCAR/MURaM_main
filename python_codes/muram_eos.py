class mu_eos:

  def __init__(self,eos_file,axis="d",table="f"):

    ## load in the merged EOS tables as used in the new extended chromosphere/coronal MURaM
    ## Filename is needed, only the new eos format is supported, with:
    ## rho = log, pre = log, eps = log, tem = log, s = lin, rhoi = log, nu_in (amb) = log, ne = log
 
    with open(eos_file, "rb") as file:
      filedata=file.read()

    ind_low = np.int(0)
    ind_high = np.int(14*4)
    temp = struct.unpack("f" * ((ind_high-ind_low)//4),filedata[ind_low:ind_high])

    self.n_eps = np.int(temp[0])
    self.n_rho = np.int(temp[1])
    self.n_p = np.int(temp[2])
    self.n_s = np.int(temp[3])

    print (self.n_eps,self.n_rho,self.n_p,self.n_s)
    
    self.eps0 = temp[4]
    self.eps1 = temp[5]
    self.lr0 = temp[6]
    self.lr1 = temp[7]
    self.lp0 = temp[8]
    self.lp1 = temp[9]
    self.s0 = temp[10]
    self.s1 = temp[11]

    self.eps_off = temp[12]
    self.ss_off = temp[13]

    self.del_lr=(self.lr1-self.lr0)/(self.n_rho-1)
    self.del_eps=(self.eps1-self.eps0)/(self.n_eps-1)
    self.del_p=(self.lp1-self.lp0)/(self.n_p-1)
    self.del_s=(self.s1-self.s0)/(self.n_s-1)
    
    if (axis is "f"):
      axis_size = 4
      axis_type = np.dtype(np.float32)
    elif (axis is "d"):
      axis_size = 8
      axis_type = np.dtype(np.float64)

    self.xrho = np.asarray(self.lr0+np.arange(self.n_rho)*self.del_lr,dtype=axis_type)
    self.xeps = np.asarray(self.eps0+np.arange(self.n_eps)*self.del_eps,dtype=axis_type)
    self.xp = np.asarray(self.lp0+np.arange(self.n_p)*self.del_p,dtype=axis_type)
    self.xs = np.asarray(self.s0+np.arange(self.n_s)*self.del_s,dtype=axis_type)

    if (table is "f"):
      table_size = 4
      table_type = np.dtype(np.float32)
    elif (table is "d"):
      table_size = 8
      table_type = np.dtype(np.float64)

    ind_low=ind_high
    ind_high+= self.n_eps*self.n_rho*table_size

    temp = struct.unpack(table * ((ind_high-ind_low)//table_size),filedata[ind_low:ind_high])
    
    self.ptbl = np.asarray(temp[:],dtype=table_type).reshape([self.n_eps,self.n_rho],order="F")

    ind_low=ind_high
    ind_high+= self.n_eps*self.n_rho*table_size

    temp = struct.unpack(table * ((ind_high-ind_low)//table_size),filedata[ind_low:ind_high])

    self.ttbl = np.asarray(temp[:],dtype=table_type).reshape([self.n_eps,self.n_rho],order="F")

    ind_low=ind_high
    ind_high+= self.n_eps*self.n_rho*table_size

    temp = struct.unpack(table * ((ind_high-ind_low)//table_size),filedata[ind_low:ind_high])

    self.stbl = np.asarray(temp[:],dtype=table_type).reshape([self.n_eps,self.n_rho],order="F")
    
    ind_low=ind_high
    ind_high+= self.n_eps*self.n_rho*table_size

    temp = struct.unpack(table * ((ind_high-ind_low)//table_size),filedata[ind_low:ind_high])

    self.netbl = np.asarray(temp[:],dtype=table_type).reshape([self.n_eps,self.n_rho],order="F")

    ind_low=ind_high
    ind_high+= self.n_eps*self.n_rho*table_size

    temp = struct.unpack(table * ((ind_high-ind_low)//table_size),filedata[ind_low:ind_high])

    self.rhoitbl = np.asarray(temp[:],dtype=table_type).reshape([self.n_eps,self.n_rho],order="F")

    ind_low=ind_high
    ind_high+= self.n_eps*self.n_rho*table_size

    temp = struct.unpack(table * ((ind_high-ind_low)//table_size),filedata[ind_low:ind_high])

    self.ambtbl = np.asarray(temp[:],dtype=table_type).reshape([self.n_eps,self.n_rho],order="F")

    ind_low=ind_high
    ind_high+= self.n_p*self.n_s*table_size

    temp = struct.unpack(table * ((ind_high-ind_low)//table_size),filedata[ind_low:ind_high])

    self.epstbl = np.asarray(temp[:],dtype=table_type).reshape([self.n_p,self.n_s],order="F")

    ind_low=ind_high
    ind_high+= self.n_p*self.n_s*table_size

    temp = struct.unpack(table * ((ind_high-ind_low)//table_size),filedata[ind_low:ind_high])

    self.rhotbl = np.asarray(temp[:],dtype=table_type).reshape([self.n_p,self.n_s],order="F")
    
    

  def interp_T(self,eps_ax,rho_ax):
    
    ee = np.log(eps_ax+self.eps_off) 
    rr = np.log(rho_ax)

    ee = np.interp(ee,self.xeps,np.arange(self.n_eps))
    rr = np.interp(rr,self.xrho,np.arange(self.n_rho))
    
    output = bilinear_interpolate(self.ttbl, ee,rr)
    output = np.exp(output)
    
    return output
  
  def interp_p(self,eps_ax,rho_ax):
    
    ee = np.log(eps_ax+self.eps_off) 
    rr = np.log(rho_ax)

    ee = np.interp(ee,self.xeps,np.arange(self.n_eps))
    rr = np.interp(rr,self.xrho,np.arange(self.n_rho))
    
    output = bilinear_interpolate(self.ptbl, ee,rr)
    output = np.exp(output)
    
    return output

  def interp_s(self,eps_ax,rho_ax):
    
    ee = np.log(eps_ax+self.eps_off) 
    rr = np.log(rho_ax)

    ee = np.interp(ee,self.xeps,np.arange(self.n_eps))
    rr = np.interp(rr,self.xrho,np.arange(self.n_rho))
    
    output = bilinear_interpolate(self.stbl, ee,rr)
    output = output-self.ss_off    

    return output
  
  def interp_ne(self,eps_ax,rho_ax):
    
    ee = np.log(eps_ax+self.eps_off) 
    rr = np.log(rho_ax)

    ee = np.interp(ee,self.xeps,np.arange(self.n_eps))
    rr = np.interp(rr,self.xrho,np.arange(self.n_rho))
    
    output = bilinear_interpolate(self.netbl, ee,rr)
    output = np.exp(output)
    
    return output
  
  def interp_amb(self,eps_ax,rho_ax):
    
    ee = np.log(eps_ax+self.eps_off) 
    rr = np.log(rho_ax)

    ee = np.interp(ee,self.xeps,np.arange(self.n_eps))
    rr = np.interp(rr,self.xrho,np.arange(self.n_rho))
    
    output = bilinear_interpolate(self.ambtbl, ee,rr)
    output = np.exp(output)
    
    return output
  
  def interp_rhoi(self,eps_ax,rho_ax):
    
    ee = np.log(eps_ax+self.eps_off) 
    rr = np.log(rho_ax)

    ee = np.interp(ee,self.xeps,np.arange(self.n_eps))
    rr = np.interp(rr,self.xrho,np.arange(self.n_rho))
    
    output = bilinear_interpolate(self.ambtbl, ee,rr)
    output = np.exp(output)
    
    return output

  def interp_eps(self,pre_ax,ss_ax):
    
    ss = ss_ax + self.ss_off
    pp = np.log(pre_ax) 

    pp = np.interp(pp,self.xp,np.arange(self.n_p))
    ss = np.interp(ss,self.xs,np.arange(self.n_s))
    
    output = bilinear_interpolate(self.epstbl, pp,ss)
    output = np.exp(output)-self.eps_off

    return output

  def interp_rho(self,pre_ax,ss_ax):
    
    ss = ss_ax + self.ss_off
    pp = np.log(pre_ax) 

    pp = np.interp(pp,self.xp,np.arange(self.n_p))
    ss = np.interp(ss,self.xs,np.arange(self.n_s))
    
    output = bilinear_interpolate(self.rhotbl, pp,ss)
    output = np.exp(output)

    return output


def bilinear_interpolate(im, x, y):
    x = np.asarray(x, dtype = 'float64')
    y = np.asarray(y, dtype = 'float64')

    x0 = np.floor(x).astype(int)
    x1 = x0 + 1
    y0 = np.floor(y).astype(int)
    y1 = y0 + 1

    x0 = np.clip(x0, 0, im.shape[0]-2);
    x1 = np.clip(x1, 1, im.shape[0]-1);
    y0 = np.clip(y0, 0, im.shape[1]-2);
    y1 = np.clip(y1, 1, im.shape[1]-1);

    Ia = im[x0,y0] # y0, x0 ]
    Ib = im[x0,y1] #y1, x0 ]
    Ic = im[x1,y0] #y0, x1 ]
    Id = im[x1,y1] #y1, x1 ]

    wa = (x1-x) * (y1-y)
    wb = (x1-x) * (y-y0)
    wc = (x-x0) * (y1-y)
    wd = (x-x0) * (y-y0)

    output = wa*Ia + wb*Ib + wc*Ic + wd*Id   
     
    return output

import numpy as np
import struct

print ("Importing muram_eos.py")
