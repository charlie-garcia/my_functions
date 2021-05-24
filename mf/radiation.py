import numpy as np
from mf.plots import create_plate2D
from mf.dcm import distance3, MatrixR
from IPython.display import clear_output

def CreateIdMatrix(x,y):
    cc = [0]*x*y
    ii=0
    for ix in range(x):
        for iy in range(y):
            cc[ii] = (ix,iy)
            ii=ii+1
    return cc

def SigmaPlate_Analytical(Lx, Ly, N):
    # Material properties: Aluminium
    rho = 2700.                                        # density         (km.m^-3)
    E   = 71e9                                         # Young Modulus   (Pa)
    nu  = 0.23                                         # Poisson coeff
    t   = 3e-3                                         # thickness       (m) 
    eta = 0.1										   # damping 
    D0  = E*t**3/( 12*(1-nu**2) )                      # Bending 

    # Compute theoretical mondes
    M = N**2
    Ne = int(np.sqrt(M))
    dx = Lx/Ne
    dy = Ly/Ne

    cs, x, y = create_plate2D(Lx,Ly,Ne, 'no')
    xx = cs[:,0];           yy = cs[:,1]
    A = dx*dy * np.ones([Ne**2,1])

    #  Modes in SS Plate
    Nm = 16;                                Nn = 15
    m  =  np.r_[1:Nm+1];                    n  =np.r_[1:Nn+1]
    km = np.array([m*np.pi/Lx]);            kn = np.array([n*np.pi/Ly]);
    # m * nx                                # n * ny
    alpha = np.sin(km.T*xx);                beta  = np.sin(kn.T*yy)          

    omega_mn  =np. sqrt(D0/(rho*t)) * (km.T**2+kn**2)
    rij   = distance3(cs, cs)                                   # distance center to center 

    # Radiation
    Nf0   = 650                                                 # N frequencies
    c0    = 343.5                                               # [m/s] c air 
    rho0  = 1.23   

    f0    = np.min(omega_mn)/(2*np.pi)
    fc    = c0**2/(2*np.pi)*np.sqrt(rho*t/D0)

    f     = np.logspace( 0, np.log10(2*10**4), Nf0 )                                           
    f     = np.unique(f.round(decimals=0))
    Nf = len(f)

    omega = 2*np.pi*f
    k0 = omega/c0
    se = np.sum(A)    

    from mf.dcm import MatrixR
    R = []
    for iif in range(Nf):
        R.append(MatrixR(rho0, c0, k0[iif], A, rij))

    fres = (np.sort(omega_mn.reshape(-1)/(2*np.pi))).reshape(Nm*Nn,1)
    omega_res = (np.sort(omega_mn.reshape(-1))).reshape(Nm*Nn,1)
    cc = CreateIdMatrix(Nm,Nn)
    sigma_mode = np.zeros([Nf, Nm*Nn], dtype='complex' )
    wm = np.zeros([1,Nm*Nn])

    for ii in range(Nm*Nn):
        print('Computing mode %.f'%(ii+1))
        id_res = ((np.abs(omega_res[ii] - omega_mn)).argmin())      # find closest freq value
        m = cc[id_res][0]
        n = cc[id_res][1]
        print('m: %.f, n %.f'%(m,n))
        mode = (alpha[m,:] * beta[n,:]).reshape(N,N)
        wm[0,ii] = omega_mn[m,n]

        vn_modes = np.array(mode.reshape(-1))
        v2m = np.mean(np.abs(vn_modes)**2)/2    

        for iif in range(0,Nf):
            # Rapid method
            W = 1/2* vn_modes.conj().T @ R[iif] @ vn_modes           # W = 1/2 * vn'*Q*L*Q'*vn; W = 1/2 * L2'*(abs(y).^2);                                      
            sigma_mode[iif,ii]  = np.real(W / (rho0*c0*se*v2m))               # Radiation Efficiency    
    # Final Reshape
        clear_output(wait=True)
    om = omega.reshape(Nf,1)                                                 # (self.Nm,self.Nn,1)
    vqm = 1 / ( (wm**2 - om**2)**2 + eta**2*wm**4 )

    sigma_mean = np.sum(sigma_mode*vqm,axis=1)/ np.sum(vqm,axis=1)
    clear_output()
    return f, sigma_mode, wm, sigma_mean, vqm

