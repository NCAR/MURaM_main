def usage():
    print ("")
    print ("Damien P.'s plotting tools")
    print ("This code has no support, or guarantees it won't ruin your data")
    print ("computer or life.")
    print ("I am bad at python. But in the interests of the death of IDL I am switching.")
    print (" Please let me know anything I can improve.")

def flines3d(bx,by,bz,xx,yy,zz,xx0,yy0,nlines,res_incr=1):

    dx = xx[2]-xx[1]
    dy = yy[2]-yy[1]
    dz = zz[2]-zz[1]

    nx = xx.size
    ny = yy.size
    nz = zz.size

    zz0 = np.zeros(nlines,dtype='float64')

    ds = dz/res_incr
    ns = int(nz*res_incr+1)

    xline = np.zeros([nlines,ns],dtype='float64')
    yline = np.zeros([nlines,ns],dtype='float64')
    zline = np.zeros([nlines,ns],dtype='float64')

    xline[:,0] = xx0
    yline[:,0] = yy0
    zline[:,0] = zz0

    for j in range(nlines):
        for i in range(1,ns):

            bx1 = trilinear_interpolate(bx,(xline[j,i-1]-xx[0])/dx,(yline[j,i-1]-yy[0])/dy,(zline[j,i-1]-zz[0])/dz)
            by1 = trilinear_interpolate(by,(xline[j,i-1]-xx[0])/dx,(yline[j,i-1]-yy[0])/dy,(zline[j,i-1]-zz[0])/dz)
            bz1 = trilinear_interpolate(bz,(xline[j,i-1]-xx[0])/dx,(yline[j,i-1]-yy[0])/dy,(zline[j,i-1]-zz[0])/dz)

            b0=np.clip(np.sqrt(bx1*bx1 + by1*by1+bz1*bz1),1.0e-10,None)

            xline[j,i]=xline[j,i-1]+ds*bx1/b0
            yline[j,i]=yline[j,i-1]+ds*by1/b0
            zline[j,i]=zline[j,i-1]+ds*bz1/b0

            if xline[j,i] < xx[0]:
                xline[j,i] = xx[0]
                yline[j,i] = yline[j,i-1]
                zline[j,i] = zline[j,i-1]

            if xline[j,i] > xx[nx-1]:
                xline[j,i] = xx[nx-1]
                yline[j,i] = yline[j,i-1]
                zline[j,i] = zline[j,i-1]

            if zline[j,i] < zz[0]:
                zline[j,i] = zz[0]
                yline[j,i] = yline[j,i-1]
                xline[j,i] = xline[j,i-1]

            if zline[j,i] > zz[nz-1]:
                zline[j,i] = zz[nz-1]
                yline[j,i] = yline[j,i-1]
                xline[j,i] = xline[j,i-1]

            if yline[j,i] < yy[0]:
                yline[j,i] = yy[0]
                zline[j,i] = zline[j,i-1]
                xline[j,i] = xline[j,i-1]

            if yline[j,i] > yy[ny-1]:
                zline[j,i] = zline[j,i-1]
                yline[j,i] = yy[ny-1]
                xline[j,i] = xline[j,i-1]

    return xline, yline, zline

def flines2d(bx,bz,xx,zz,xx0,res_incr=1):

    dx = xx[2]-xx[1]
    dz = zz[2]-zz[1]

    nx = xx.size
    nz = zz.size

    nlines = xx0.size
    zz0 = np.zeros(nlines,dtype='float64')

    ds = dz/res_incr
    ns = int(nz*res_incr+1)*4

    xline = np.zeros([nlines,ns],dtype='float64')
    zline = np.zeros([nlines,ns],dtype='float64')

    xline[:,0] = xx0
    zline[:,0] = zz0
    
    for j in range(nlines):
        for i in range(1,ns):
             
            bx1 = bilinear_interpolate(bx,(xline[j,i-1]-xx[0])/dx,(zline[j,i-1]-zz[0])/dz)
            bz1 = bilinear_interpolate(bz,(xline[j,i-1]-xx[0])/dx,(zline[j,i-1]-zz[0])/dz)
            b0 = np.clip(np.sqrt(bx1*bx1 + bz1*bz1),1.0e-10,None)

            xline[j,i]=xline[j,i-1]+ds*bx1/b0
            zline[j,i]=zline[j,i-1]+ds*bz1/b0

            if xline[j,i] < xx[0]:
                xline[j,i] = xx[0]
                zline[j,i] = zline[j,i-1]

            if xline[j,i] > xx[nx-1]:
                xline[j,i] = xx[nx-1]
                zline[j,i] = zline[j,i-1]

            if zline[j,i] < zz[0]:
                zline[j,i] = zz[0]
                xline[j,i] = xline[j,i-1]
            
            if zline[j,i] > zz[nz-1]:
                zline[j,i] = zz[nz-1]
                xline[j,i] = xline[j,i-1]

    return xline, zline

def deriv_nd_O2(arr,dir,delta=1.0):
    ## Do a derivative in 3 dim array
    ## Dir = array indices
    ## delta = grid spacing

    dim = arr.ndim
    if dim > 3:
        print("Currently coded for only a 3-dim array")
        print(arr.ndim)
        return None

    derr = np.zeros(arr.shape,dtype = np.float64)
    derr = -0.5*np.roll(arr,1,dir) + 0.5*np.roll(arr,-1,dir)

    if dir == 0:
        derr[0,:,:] = -1.5*arr[0,:,:] + 2.0*arr[1,:,:] - 0.5*arr[2,:,:]
        derr[-1,:,:] = 1.5*arr[-1,:,:] - 2.0*arr[-2,:,:] + 0.5*arr[-3,:,:]
    if dir == 1:
        derr[:,0,:] = -1.5*arr[:,0,:] + 2.0*arr[:,1,:] - 0.5*arr[:,2,:]
        derr[:,-1,:] = 1.5*arr[:,-1,:] - 2.0*arr[:,-2,:] + 0.5*arr[:,-3,:]
    if dir == 2:
        derr[:,:,0] = -1.5*arr[:,:,0] + 2.0*arr[:,:,1] - 0.5*arr[:,:,2]
        derr[:,:,-1] = 1.5*arr[:,:,-1] - 2.0*arr[:,:,-2] + 0.5*arr[:,:,-3]
   
    derr = derr/delta

    return derr

def plotvhslice(a,arr,xyslice,xzslice,xpos,ypos,xax,yax,zax,xrange,yrange,zrange, br = None,title=None,btitle='',beta=None,flines=None,cmap='jet',sym=False,animated=False,cbar=True,xtitle=True,ytitle=True,asp=[1,1]):
    
    [nx,ny,nz] = arr.shape

    r_x = 0
    r_y = 0
    r_z = xyslice
    d_x = nx
    d_y = ny
    d_z = r_z + 1
    x = xax
    y = yax
    xr = xrange
    yr = yrange
    xttl = None
    if ytitle:
        yttl = 'Y (Mm)'
    else:
        yttl = None

    ttl = title

    if beta is not None:
        bslice = beta[r_x:d_x,r_y:d_y,r_z:d_z]
        Y,X = np.meshgrid(np.linspace(y[0],y[-1],bslice.squeeze().shape[1]),
np.linspace(x[0],x[-1],bslice.squeeze().shape[0]))

    if flines is not None:
        nlines = flines[0].shape[0]

    if br is None:
      arr_min = arr[r_x:d_x,r_y:d_y,r_z:d_z].min()
      arr_max = arr[r_x:d_x,r_y:d_y,r_z:d_z].max()
    else:
      arr_min = br[0]
      arr_max = br[1]

    if cbar:
       if cbar is 'Top':
         cbar1 = 'Top'
         cbar2 = False
       elif cbar is 'Bottom':
         cbar1 = False
         cbar2 = True
       elif cbar is 'Right':
         cbar1 = 'Right'
         cbar2 = 'Right'
       elif cbar is 'Left':
         cbar1 = 'Left'
         cbar2 = 'Left'
       else:
         cbar1 = True
         cbar2 = True
    else:
      cbar1 = False
      cbar2 = False

    [im1,ax,cbt1,cb,cax] = plotsnapshot(a,arr[r_x:d_x,r_y:d_y,r_z:d_z].squeeze(),title=title,
xax=[x[r_x],x[d_x-1]],yax=[y[r_y],y[d_y-1]],xr=xr,yr=yr,xtitle = xttl, ytitle = yttl,arr_min = arr_min, arr_max = arr_max,btitle=btitle, sym=sym,pos=xpos,asp = asp[0],cmap=cmap,animated=animated,cbar=cbar1)
    if beta is not None:
        CF = ax.contour(bslice.squeeze().T,(1.0,),colors = "black",extent=[x[0],x[-1],y[0],y[-1]]
,linestyles='dashed',linewidth=0.3)
    ax.plot(xr,[y[xzslice],y[xzslice]])   
   
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%02.1f'))
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%02.1f'))

    r_x = 0
    r_y = xzslice
    r_z = 0
    d_x = nx
    d_y = r_y + 1
    d_z = nz
    x = xax
    y = zax
    xr = xrange
    yr = zrange
    if xtitle:
        xttl = 'X (Mm)'
    else:
        xttl = None
    if ytitle:
        yttl = 'Z (Mm)'
    else:
        yttl = None
    ttl = 'y ' + "{0:.3f}".format(y[r_y]*1.0e-8)

    if beta is not None:
        bslice = beta[r_x:d_x,r_y:d_y,r_z:d_z]
        Y,X = np.meshgrid(np.linspace(y[0],y[-1],bslice.squeeze().shape[1]),
np.linspace(x[0],x[-1],bslice.squeeze().shape[0]))

    if br is None:
      arr_min = arr[r_x:d_x,r_y:d_y,r_z:d_z].min()
      arr_max = arr[r_x:d_x,r_y:d_y,r_z:d_z].max()
    else:
      if (np.size(br) == 2):
        arr_min = br[0] 
        arr_max = br[1]
      elif (np.size(br) == 4):
        arr_min = br[2]
        arr_max = br[3]
      else:
        arr_min = arr[r_x:d_x,r_y:d_y,r_z:d_z].min()
        arr_max = arr[r_x:d_x,r_y:d_y,r_z:d_z].max()


    [im2,ax,cbt2, cb,cax] = plotsnapshot(a,arr[r_x:d_x,r_y:d_y,r_z:d_z].squeeze(),
xax=[x[r_x],x[d_x-1]],yax=[y[r_z],y[d_z-1]],xr=xr,yr=yr,title = None, xtitle = xttl, ytitle = yttl,btitle=btitle, arr_min = arr_min, arr_max = arr_max,sym=sym,pos=ypos,asp = asp[1],cmap=cmap,animated=animated,cbar=cbar2)
    if beta is not None:
        CF = ax.contour(bslice.squeeze().T,(1.0,),colors = "black",extent=[x[0],x[-1],y[0],y[-1]]
,linestyles='dashed',linewidth=0.3)
    if flines is not None:
        for i in range(nlines):
           ax.plot(flines[0][i,:]*1.0e-8,flines[1][i,:]*1.0e-8,'w')

    ax.plot(xr,[y[xyslice],y[xyslice]])
    
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%02.1f'))
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%02.1f'))

    return im1,im2,cbt1,cbt2

def plotsnapshot(fig,arr,xax=None,yax=None,xax2=None,yax2=None,pos=0,xr=None,yr=None, title = None , xtitle = None ,ytitle = None,btitle = '', arr_min = None, arr_max = None,nticks=5,sym = False,asp = 'auto',cmap='jet',animated=False,cbar=True,bsize=None,interp='None'):
    ## Add some kind of command to do the aspect

    if arr.ndim != 2: 
        print ("Array not 2-dimensional! arr.ndim = ", arr.ndim)
        return None

    nx = arr.shape[0]
    ny = arr.shape[1]

    if xax is None:
        xax = np.arange(nx)
    elif len(xax)==2:
        xax = np.arange(nx)/(nx-1)*(xax[1]-xax[0])+xax[0]
    elif len(xax) == nx:
        xax = xax
    else:
        print ("Problem, xax != nx")
        xax = np.arange(nx)
        
    if yax is None:
        yax = np.arange(ny)
    elif len(yax)==2:
        yax = np.arange(ny)/(ny-1)*(yax[1]-yax[0])+yax[0]
    elif len(yax) == ny:
        yax = yax
    else:
        print ("Problem, yax != ny")
        yax = np.arange(ny)
   
    if xr is None:
       ind_x = np.arange(nx)
    else:
       ind_x = np.where((xax[:] > xr[0]) & (xax[:] < xr[1]))[0]
            

    if yr is None:
        ind_y = np.arange(ny)
    else:
        ind_y = np.where((yax[:] > yr[0]) & (yax[:] < yr[1]))[0]

    if arr_min == None: 
        arr_min = np.amin(arr[ind_x[0]:ind_x[-1],ind_y[0]:ind_y[-1]])
        if arr_min == 0:
            arr_min = 0.0
        if arr_min != 0:
            a_mino = np.floor(np.log10(np.abs(arr_min)))
            arr_min = np.floor(arr_min/10.0**(a_mino-1))*10.0**(a_mino-1)

    if arr_max == None:
        arr_max = np.amax(arr[ind_x[0]:ind_x[-1],ind_y[0]:ind_y[-1]])
        if arr_max == 0:
            arr_max = 0.0
        if arr_max != 0:
            a_maxo = np.floor(np.log10(np.abs(arr_max)))
            arr_max = np.ceil(arr_max/10.0**(a_maxo-1))*10.0**(a_maxo-1)

    if sym:
        arr_max = max(abs(arr_max),abs(arr_min))
        arr_min = -arr_max

    if arr_max == arr_min:
        cbarticks = [arr_min,arr_max]
    if arr_max != arr_min:
        cbarticks = np.linspace(arr_min, arr_max, nticks, endpoint=True)

    ax = fig.add_subplot(pos)
    if title is not None:
        ax.set_title(title)
    if xtitle is not None:
        ax.set_xlabel(xtitle)
    if ytitle is not None:
        ax.set_ylabel(ytitle)

    img = ax.imshow(arr.T,interpolation=interp, vmin=arr_min, vmax=arr_max,origin='lower',
extent=[xax[0],xax[-1],yax[0],yax[-1]], aspect=asp,cmap=plt.get_cmap(cmap),animated=animated)
    if xr:
        ax.set_xlim(xr)
    if yr:
        ax.set_ylim(yr)
        
            
    ax2 = 0
    ax3 = 0
    
    if xax2 is not None:     
      ax2 = ax.twiny()
      ax2.set_xlim(ax.get_xlim())
      ax2.set_xticks(xax2[0])
      ax2.set_xticklabels(xax2[1])
      ax2.set_xlabel(xax2[2])
     
    if yax2 is not None:
      ax3 = ax.twinx()
      ax3.set_ylim(ax.get_ylim())
      ax3.set_yticks(yax2[0])
      ax3.set_yticklabels(yax2[1])
      ax3.set_ylabel(yax2[2],labelpad=4)
      plt.draw() 
        
    divider = make_axes_locatable(ax)
               
    if cbar:
      if (cbar is True) or (cbar is 'Right'):
        cax = divider.append_axes("right", size="2%", pad=0.1)        
        cb = fig.colorbar(img,cax=cax,ticks = cbarticks, format=ticker.FuncFormatter(fmt))
      if cbar is 'Left':
        cax = divider.append_axes("left", size="2%", pad=0.1)
        cb = fig.colorbar(img,cax=cax,ticks = cbarticks, format=ticker.FuncFormatter(fmt))
      if cbar is 'Top':
        cax = divider.append_axes("top", size="5%", pad=0.1)
        cb = fig.colorbar(img,cax=cax,ticks = cbarticks,orientation='horizontal', format=ticker.FuncFormatter(fmt))
      if cbar is 'Bottom':
        cax = divider.append_axes("bottom", size="5%", pad=0.1)
        cb = fig.colorbar(img,cax=cax,ticks = cbarticks,orientation='horizontal', format=ticker.FuncFormatter(fmt))
      cb.set_label(btitle, labelpad=4)
    else:
      cax = []
      cb = []
    
    return img,[ax,ax2,ax3],cbarticks,cb,cax

def trilinear_interpolate(im,x,y,z):
    
  x = np.asarray(x, dtype = 'float64')
  y = np.asarray(y, dtype = 'float64')
  z = np.asarray(z, dtype = 'float64')

  x0 = np.floor(x).astype(int)
  x1 = x0 + 1
  y0 = np.floor(y).astype(int)
  y1 = y0 + 1
  z0 = np.floor(z).astype(int)
  z1 = z0+1    

  x0 = np.clip(x0, 0, im.shape[0]-1);
  x1 = np.clip(x1, 0, im.shape[0]-1);
  y0 = np.clip(y0, 0, im.shape[1]-1);
  y1 = np.clip(y1, 0, im.shape[1]-1);
  z0 = np.clip(z0, 0, im.shape[2]-1);
  z1 = np.clip(z1, 0, im.shape[2]-1);

  Iaaa = im[x0,y0,z0] # y0, x0 ]
  Iaab = im[x0,y0,z1] #y1, x0 ]
  Iaba = im[x0,y1,z0] #y0, x1 ]
  Iabb = im[x0,y1,z1] #y1, x1 ]
  Ibaa = im[x1,y0,z0] # y0, x0 ]
  Ibab = im[x1,y0,z1] #y1, x0 ]
  Ibba = im[x1,y1,z0] #y0, x1 ]
  Ibbb = im[x1,y1,z1] #y1, x1 ]

  waaa = (x1-x) * (y1-y) * (z1-z)
  waab = (x1-x) * (y1-y) * (z-z0)
  waba = (x1-x) * (y-y0) * (z1-z)
  wabb = (x1-x) * (y-y0) * (z-z0)
  wbaa = (x-x0) * (y1-y) * (z1-z)
  wbab = (x-x0) * (y1-y) * (z-z0)
  wbba = (x-x0) * (y-y0) * (z1-z)
  wbbb = (x-x0) * (y-y0) * (z-z0)

  return waaa*Iaaa + waab*Iaab + waba*Iaba + wabb*Iabb + wbaa*Ibaa + wbab*Ibab + wbba*Ibba + wbbb*Ibbb

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

def fmt(x, pos):
  a, b = '{:.2e}'.format(x).split('e')
  b = int(b)

  if (abs(x) == 0.0):
    outfmt = '{:03.2f}'.format(x)
  elif((abs(x) >= 1000.0) or (abs(x) <= 1.0e-3)):
    outfmt = r'${}\times10^{{{}}}$'.format(a, b)
  elif (abs(x) >= 100.0):
    outfmt = '{:04.1f}'.format(x)
  elif (abs(x) >= 10.0):
    outfmt = '{:03.1f}'.format(x)
  elif (abs(x) < 1.0e-2):
    outfmt = '{:04.3f}'.format(x)
  else:
    outfmt = '{:03.2f}'.format(x)
  return outfmt

def fig_open(figsize=None,numx=1,numy=1,hr = None, wr = None):
  fig = plt.figure(figsize=figsize)
  if not hr:
    hr = np.tile(1,numy)
  if not wr:
    wr = np.tile(1,numx)
  gs = gridspec.GridSpec(numy,numx,height_ratios=hr,width_ratios=wr)
  return fig,gs

def fig_saveandcls(fig,name,dir,dpi=20):
  fig.savefig(dir+name,dpi=dpi)
  fig.clf()
  plt.close()
  return

def test_MHS(bkg,ptb = None):
  if ptb:
      bx = bkg.bx + ptb.bx
      by = bkg.by + ptb.by
      bz = bkg.bz + ptb.bz
      rho = bkg.rho + ptb.rho
      pe  = bkg.pe + ptb.pe
  else:
      bx = bkg.bx
      by = bkg.by
      bz = bkg.bz
      rho = bkg.rho
      pe  = bkg.pe

  dbxdy = deriv_nd_O2(bx,1,delta=bkg.y[1]-bkg.y[0])
  dbzdy = deriv_nd_O2(bz,1,delta=bkg.y[1]-bkg.y[0])
  dbydx = deriv_nd_O2(by,0,delta=bkg.xax[1]-bkg.xax[0])
  dbzdx = deriv_nd_O2(bz,0,delta=bkg.xax[1]-bkg.xax[0])
  dbxdz = deriv_nd_O2(bx,2,delta=bkg.zax[1]-bkg.zax[0])
  dbydz = deriv_nd_O2(by,2,delta=bkg.zax[1]-bkg.zax[0])

  dpdx = deriv_nd_O2(pe,0,delta=bkg.xax[1]-bkg.xax[0])
  dpdy = deriv_nd_O2(pe,1,delta=bkg.yax[1]-bkg.yax[0])
  dpdz = deriv_nd_O2(pe,2,delta=bkg.zax[1]-bkg.zax[0])

  jx =  (dbzdy-dbydz)
  jy = (-dbzdx+dbxdz)
  jz =  (dbydx-dbxdy)

  eq_x = (dpdx - (jy*bz - jz*by))/abs(dpdx+(jy*bz - jz*by))
  eq_y = (dpdy - (-jx*bz + jz*bx))/abs(dpdy+(-jx*bz - jz*bx))
  eq_z = (dpdz - (jx*by - jy*bx) + rho*2.74e4)/abs(dpdz+(jx*by - jy*bx)-rho*2.74e4)

def inttostring(ii,ts_size=4):
## For an integer ii, and a timestamp size ts_size, create the string needed
## for loading the file, i.e. Alfven_0000

  str_num = str(ii)

  for bb in range(len(str_num),4,1):
        str_num = '0' + str_num

  return str_num

import math
import numpy as np
import sys
from subprocess import call
import time
import getopt
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.animation as animation
from matplotlib import gridspec
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.interpolate import griddata

print ("Importing dp_plot_tools.py")

