import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as inter
from fenics import *
# from fenics import *

def integralPath(x,y,u):
    z = np.sqrt(x**2 + y**2)
  
    f = inter.interp1d(z, u, kind='cubic')
    new_z = np.linspace(z.min(), z.max())
    new_u  = f(new_z)
     
    nds = new_z[1]-new_z[0]
    sol = np.trapz(new_u, dx = nds)
 
    return sol

def check_array(x):
    try:
        x.shape
        return True
    except:
        return False
        
def dbW(signal):
	sdB = 10*np.log10(np.abs(signal)/10e-12)
	return sdB


def RotationMatrix(theta):
    return np.array([[np.cos(theta), np.sin(theta)], [-np.sin(theta), np.cos(theta)]])

def RotatePoint(x,y,theta):
    X = x*np.cos(theta) - y*np.sin(theta) 
    Y = x*np.sin(theta) + y*np.cos(theta)
    return X,Y

def InertiaMatrix(fmesh, rho):   
    dx = Measure("dx")(domain=fmesh)
    V = VectorFunctionSpace(fmesh, 'CG', 1)
    v = interpolate(Expression(('x[0]', 'x[1]'), degree=2), V)
    
    Mx = assemble(rho*v[1]*dx) # integrate over y
    My = assemble(rho*v[0]*dx) # integrate over x
    
    x0 = My/assemble(rho*dx)
    y0 = Mx/assemble(rho*dx)
    
    #% Get eigenvalues
    Ix  = assemble(v[1]**2*rho*dx)
    Iy  = assemble(v[0]**2*rho*dx)
    Ixy = -assemble(v[0]*v[1]*rho*dx)
    
    N =2
    # Construct the inertia matrix
    I = np.zeros((N,N))
    I[0][0] =  Ix
    I[1][1] =  Iy
    I[0][1] = Ixy
    I[1][0] = Ixy
    
    # Get the eigenvalues adn eigenvectors
    eig = np.linalg.eig(I)
    mode = np.c_[eig[1][:,0], eig[1][:,1]]
    
    return I, x0, y0, mode