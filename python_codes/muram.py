"""MURaM output analysis

This package contains functions and classes for reading in MURaM
simulation output.  Data objects subclass numpy.ndarray in a way that
represnets the underlying data directly, but provied named attributes
for each physical variable, so that users are not required to remember
indexes.  2D data is loaded in-memory, while 3D datacubes use the
numpy.memmap to keep data on disk.  Data objects have meta-data
attributes such as iteration number and simulation time.
"""

import numpy as np
import os.path
from collections import OrderedDict

def read_hmean(datapath, iteration):
    filename = f"hmean1D.{iteration:06d}"
    filepath = os.path.join(datapath, filename)
    data = np.fromfile(filepath, dtype=np.float32)
    Nvar = data[0].astype(int) # this is 1
    shape = tuple(data[1:3].astype(int))
    time = data[3]
    Iout = data[4:].reshape([Nvar,shape[0]]).swapaxes(0,1)
    return Iout, Nvar, shape[0], time

def read_Iout(datapath, iteration):
    filename = f"I_out.{iteration:06d}"
    filepath = os.path.join(datapath, filename)
    data = np.fromfile(filepath, dtype=np.float32)
    Nvar = data[0].astype(int) # this is 1
    shape = tuple(data[1:3].astype(int))
    time = data[3]
    Iout = data[4:].reshape([shape[1], shape[0]]).swapaxes(0,1)
    return Iout, shape, time

def read_slice(datapath, iteration, kind, ix):
    filename = f"{kind}_slice_{ix}.{iteration:06d}"
    filepath = os.path.join(datapath, filename)
    data = np.fromfile(filepath, dtype=np.float32)
    Nvar = data[0].astype(int)
    shape = tuple(data[1:3].astype(int))
    time = data[3]
    slice = data[4:].reshape([Nvar, shape[1], shape[0]]).swapaxes(1,2)
    return slice, Nvar, shape, time

class MuramIntensity(np.ndarray):
    """Horizontal slice of intensity at tau=1
    
    This class inherits from ndarray and may be used as a normal NumPy array.  The 
    dimensions are (y, z) and the values are the intensity at tau=1.
    """    
    def __new__(cls, datapath, iteration):
        """Instantiates a new MuramIntensity object
        
           datapath <str> := path to a MURaM output directory containing I_out files
           iteration <int> :=  iteration number
        """
        Iout, shape, time = read_Iout(datapath, iteration)
        self = np.array(Iout).view(cls)
        self.datapath = datapath
        self.iteration = iteration
        self.time = time
        self.shape = shape
        return self
    
    def __init__(self, datapath, iteration):
        """datapath (str): Path to a MURaM output directory containing xy_slice files
          iteration (int): Iteration number
        """
        pass # this is for documentation only; actual instantiation done in __new__

class MuramSlice(np.ndarray):
    """2D slice through a MURaM domain
    
    This is a generic class that may represent horizontal, vertical, or tau slices.
    
    This object subclasses ndarray and may be used as if it were a normal NumPy array.  
    The index ordering is [var, dim1, dim2]. The variable attributes (e.g. slice.rho) have 
    dimension (dim1, dim2). The transpose option is for the xy_slice files that order the
    data by (var, y, x), which is counter-intuitive for something called "xy".
    
    Args:
      datapath (str): Path to a MURaM output directory containing slice files
      kind (str): Type of slice (xy, yz, or tau)
      ix (str): Index of the dimension orthogonal to the slice for (xy, yz), or tau value
      iteration (int): Iteration number
      transpose (bool): Whether to transpose the slice dimensions on instantiation.
      
    Attributes:
      ix (str): Index of the dimension orthogonal to the slice for (xy, yz), or tau value
      iteration (int): Iteration number
      time (float): Seconds since start of simulation
      rho (ndarray): density
      eint (ndarray): internal energy
      Temp (ndarray): temperature
      Pres (ndarray): pressure
      vx (ndarray): velocity in the x-direction
      vy (ndarray): velocity in the y-direction
      vz (ndarray): velocity in the z-direction
      Bx (ndarray): magnetic field in the x-direction
      By (ndarray): magnetic field in the y-direction
      Bz (ndarray): magnetic field in the z-direction
    """    
    def __new__(cls, datapath, iteration, kind, ix, transpose=False):
        sl, Nvar, shape, time = read_slice(datapath, iteration, kind, ix)
        if transpose:
            sl = np.transpose(sl, [0, 2, 1])
        self = np.array(sl).view(cls)
        self.ix = ix
        self.iteration = iteration
        self.time = time
        self.rho = self[0,:,:]
        self.vx = self[1,:,:]
        self.vy = self[2,:,:]
        self.vz = self[3,:,:]
        self.eint = self[4,:,:]
        self.Bx = self[5,:,:] * np.sqrt(4 * np.pi)
        self.By = self[6,:,:] * np.sqrt(4 * np.pi)
        self.Bz = self[7,:,:] * np.sqrt(4 * np.pi)
        # TODO: bug based on convention.
        #   We assume that the above are always written, and the following are sequentially 
        #   enabled.  The real way to do this is to read the configuration and figure out
        #   what was really written out.
        if Nvar > 8:
            self.Temp = self[8,:,:]
        if Nvar > 9:
            self.Pres = self[9,:,:]
        return self
    
class MuramTauSlice(MuramSlice):
    """Slice at a constant optical depth"""
    def __new__(cls, datapath, iteration, tau):
        if(tau>=1e-3):
            self = MuramSlice.__new__(cls, datapath, iteration,  'tau', f'{tau:05.3f}')
        else:
            self = MuramSlice.__new__(cls, datapath, iteration,  'tau', f'{tau:08.6f}')
        self.tau = tau
        return self
    
    def __init__(self, datapath, iteration, tau):
        """datapath (str): Path to a MURaM output directory containing xy_slice files
          tau (float): tau value for slice
          iteration (int): Iteration number
        """
        pass # this is for documentation only; actual instantiation done in __new__


class MuramColumn:
    """Vertical column from an xy slice"""
    
    def __init__(self, xy_slice, yix):
        """xy_slice (MuramXYSlice): vertical slice
           yix (int): y-index from which to select a column (x-array)
        """
        self.yix = yix
        self.time = xy_slice.time
        self.rho = xy_slice.rho[:, yix]
        self.vx = xy_slice.vx[:, yix]
        self.vy = xy_slice.vy[:, yix]
        self.vz = xy_slice.vz[:, yix]
        self.eint = xy_slice.eint[:, yix]
        self.Bx = xy_slice.Bx[:, yix]
        self.By = xy_slice.By[:, yix]
        self.Bz = xy_slice.Bz[:, yix]
        if xy_slice.shape[0] > 8:
            self.Temp = xy_slice.T[:, yix]
        if xy_slice.shape[0] > 9:
            self.Pres = xy_slice.P[:, yix]
        
class MuramXYSlice(MuramSlice):
    """Vertical slice at a constant z-index"""
    
    def __new__(cls, datapath, z, iteration):
        self = MuramSlice.__new__(cls, datapath, 'xy', f'{z:04d}', iteration, transpose=True)
        self.z = z
        return self
    
    def __init__(self, datapath, z, iteration):
        """datapath (str): Path to a MURaM output directory containing xy_slice files
          z (int): z-index of slice
          iteration (int): Iteration number
        """
        pass # this is for documentation only; actual instantiation done in __new__
    
    def column(self, yix):
        """Returns a column in depth (x) for a given y-index"""
        return MuramColumn(self, yix)


class MuramCube(np.memmap):
    """A cube (x, y, z) of MURaM output
    
    This class inherits from NumPy memmap and uses that to reprresent the data.  The data is 
    not loaded into memory, it is kept on file and read in as required.  For example, 
    indexing cube[:,0,0] will read and return a column of data at the coordinate y=0, z=0, 
    without ever reading the whole cube into memory.  If in-memory data is desired, use 
    .copy().
    
    Note that the first dimension (x) is the vertical component.
    
    Args:
      datapath (str):   Path to a MURaM output directory containing data cubes
      iteration (int):  iteration number
      var (str):        variable name (e.g. 'rho'). See MuramCube.varnames for a list

    Attributes:
        var (str):       Name of the variable the cube represnets (e.g. 'rho')
        filepath (str):  Path t the data file
        iteration (int): iteration number
        dX (ndarray):    Pixel scale in cm for x, y, z dimensions
        X (ndarray):     Domain size in cm for x, y, z dimensions
        time (float):    Seconds since start of simulation
        vamax (float):   Maximum Alfven velocity (TODO: varify)
    """
    _varmap = OrderedDict((('rho','result_prim_0'),
                       ('vx','result_prim_1'),
                       ('vy','result_prim_2'),
                       ('vz','result_prim_3'),
                       ('eint','result_prim_4'),
                       ('Bx','result_prim_5'),
                       ('By','result_prim_6'),
                       ('Bz','result_prim_7'),
                       ('sflx','result_prim_8'),
                       ('E_nt','result_prim_9'),
                       ('E_nt_flx','result_prim_10'),
                       ('Temp','eosT'),
                       ('Pres','eosP'),
                       ('ne','eosne'),
                       ('amb','eosamb'),
                       ('rhoi','eosrhoi'),
                       ('Q','Qtot'),
                       ('J','Jtot'),
                       ('S','Stot'),
                       ('tau','tau'),
                       ('QxH','QxH'),
                       ('QxCa','QxCa'),
                       ('QxMg','QxMg'),
                       ('QxCor','QxCor'),
                       ('QxChr','QxChr'),
                       ('Qres','Qres'),
                       ('Qvis','Qvis'),
                       ('Qamb','Qamb'),
                       ('tvar1','tvar1'),
                       ('tvar2','tvar2'),  
                       ('tvar3','tvar3'),
                       ('tvar4','tvar4'),
                       ('tvar5','tvar5'),  
                       ('tvar6','tvar6'),
                       ('tvar7','tvar7'),
                       ('tvar8','tvar8'),      
                       ('alpha','alpha'),
                       ('dips','dips'),
                       ('kappa_trac','kappa_trac'),
                       ('jabs','jabs'),
                       ('JxB_abs','Fa'),
                       ('vJxB','vFl'),
                       ('JxB_norm','Fn'),    
                       ('decay_ind','decay_ind')))
    
    varnames = tuple(_varmap.keys())
    
    def _read_header(datapath, iteration):
        filename = f"Header.{iteration:06d}"
        filepath = os.path.join(datapath, filename)
        header = np.loadtxt(filepath, dtype=np.float32)
        shape = header[0:3].astype(int)
        dX = header[3:6] # cm; pixel size
        X = shape * dX # cm; domain size
        time = header[6] # s; time since t=0
        dtime = header[7] # s; timestep
        vamax = header[8] # TODO: what is this? Alfven velocity?
        return shape, dX, X, time, dtime, vamax
    
    def _filepath(datapath, iteration, var):
        filepre = MuramCube._varmap[var]
        filename = f"{filepre}.{iteration:06d}"
        filepath = os.path.join(datapath, filename)
        return filepath
    
    def __new__(cls, datapath, iteration, var):
        shape, dX, X, time, dtime, vamax = cls._read_header(datapath, iteration)
        filepath = cls._filepath(datapath, iteration, var)
        # Note: binary data index order is (z, y, x). Here we reverse the true shape 
        # to reflect that, then transpose the data in reverse index order to be (x, y, z)
        self = np.memmap(filepath, mode='r', dtype=np.float32,
                         shape=tuple(shape[::-1])).transpose(2, 1, 0)
        self.var = var
        self.filepath = filepath
        self.iteration = iteration
        self.dX = dX
        self.X = X
        self.time = time
        self.dtime = dtime
        self.vamax = vamax
        return self
    
class MuramSnap():
    """Snapshot of MURaM output; a collection of data cubes

    This class prepares MuramCube objects for each known output type that is found 
    in the data directory for the given iteration.  The MuramCube objects are available 
    as class attributes with the same name as found in MuramCube.varnames.  Variable 
    types which are not found are silently ignored and the attribute will not be available.  
    The 'available' and 'unavailable' attributes contain lists of variables which were 
    found or not found, respectively.
    
    Args:
      datapath (str):   Path to a MURaM output directory containing data cubes
      iteration (int):  iteration number
    """
    def __init__(self, datapath, iteration):
        self.available = []
        self.unavailable = []
        for var in MuramCube._varmap.keys():
            try:
                cube = MuramCube(datapath, iteration, var)
            except FileNotFoundError:
                self.unavailable.append(var)
                continue
            setattr(self, var, cube)
            self.available.append(var) 

class MuramSubCube(np.memmap):
    """A cube (x, y, z) of MURaM output
    
    This class inherits from NumPy memmap and uses that to reprresent the data.  The data is 
    not loaded into memory, it is kept on file and read in as required.  For example, 
    indexing cube[:,0,0] will read and return a column of data at the coordinate y=0, z=0, 
    without ever reading the whole cube into memory.  If in-memory data is desired, use 
    .copy().
    
    Note that the first dimension (x) is the vertical component.
    
    Args:
      datapath (str):   Path to a MURaM output directory containing data cubes
      iteration (int):  iteration number
      var (str):        variable name (e.g. 'rho'). See MuramCube.varnames for a list

    Attributes:
        var (str):       Name of the variable the cube represnets (e.g. 'rho')
        filepath (str):  Path t the data file
        iteration (int): iteration number
        dX (ndarray):    Pixel scale in cm for x, y, z dimensions
        X (ndarray):     Domain size in cm for x, y, z dimensions
        time (float):    Seconds since start of simulation
        vamax (float):   Maximum Alfven velocity (TODO: varify)
    """
    _varmap = OrderedDict((('rho','subdomain_0'),
                       ('vx','subdomain_1'),
                       ('vy','subdomain_2'),
                       ('vz','subdomain_3'),
                       ('eint','subdomain_4'),
                       ('Bx','subdomain_5'),
                       ('By','subdomain_6'),
                       ('Bz','subdomain_7'),
                       ('Temp','subdomain_8'),
                       ('Pres','subdomain_9'),
                       ('ne','subdomain_10'),
                       ('tau','subdomain_11')))
    
    varnames = tuple(_varmap.keys())
    
    def _read_header(datapath, iteration):
        filename = f"Header_Sub.{iteration:06d}"
        filepath = os.path.join(datapath, filename)
        header = np.loadtxt(filepath, dtype=np.float32)
        shape = header[0:3].astype(int)
        dX = header[3:6] # cm; pixel size
        X = shape * dX # cm; domain size
        time = header[6] # s; time since t=0
        dtime = header[7] # s; timestep
        vamax = header[8] # TODO: what is this? Alfven velocity?
        return shape, dX, X, time, dtime, vamax
    
    def _filepath(datapath, iteration, var):
        filepre = MuramSubCube._varmap[var]
        filename = f"{filepre}.{iteration:06d}"
        filepath = os.path.join(datapath, filename)
        return filepath
    
    def __new__(cls, datapath, iteration, var):
        shape, dX, X, time, dtime, vamax = cls._read_header(datapath, iteration)
        filepath = cls._filepath(datapath, iteration, var)
        # Note: binary data index order is (z, y, x). Here we reverse the true shape 
        # to reflect that, then transpose the data in reverse index order to be (x, y, z)
        self = np.memmap(filepath, mode='r', dtype=np.float32,
                         shape=tuple(shape[::-1])).transpose(2, 1, 0)
        self.var = var
        self.filepath = filepath
        self.iteration = iteration
        self.dX = dX
        self.X = X
        self.time = time
        self.dtime = dtime
        self.vamax = vamax
        return self
    
class MuramSubSnap():
    """Snapshot of MURaM output; a collection of data cubes

    This class prepares MuramCube objects for each known output type that is found 
    in the data directory for the given iteration.  The MuramCube objects are available 
    as class attributes with the same name as found in MuramCube.varnames.  Variable 
    types which are not found are silently ignored and the attribute will not be available.  
    The 'available' and 'unavailable' attributes contain lists of variables which were 
    found or not found, respectively.
    
    Args:
      datapath (str):   Path to a MURaM output directory containing data cubes
      iteration (int):  iteration number
    """
    def __init__(self, datapath, iteration):
        self.available = []
        self.unavailable = []
        for var in MuramSubCube._varmap.keys():
            try:
                cube = MuramSubCube(datapath, iteration, var)
            except FileNotFoundError:
                self.unavailable.append(var)
                continue
            setattr(self, var, cube)
            self.available.append(var) 
