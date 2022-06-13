import numpy as np
from mf.plots import create_plate2D
from mf.dcm import distance3, MatrixR
from IPython.display import clear_output
import scipy.special as sp

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

# Radiation efficiency in a simply supported plate
def sigma_mn(m,n,k0,Lx,Ly):
    from scipy.integrate import dblquad

    # limits for theta
    t1 = 0
    t2 = np.pi/2
    
    # limits for phi
    p1 = 0
    p2 = np.pi/2

    def wallace(phi, theta):
        alpha = k*Lx*np.sin(theta)*np.cos(phi)
        beta  = k*Ly*np.sin(theta)*np.sin(phi)
        if (m % 2) ==0 and (n % 2) ==0 :
            return ( np.sin(alpha/2)*np.sin(beta/2) /  ( ( (alpha/(m*np.pi))**2-1)*((beta/(n*np.pi))**2-1) ))**2  * np.sin(theta)
        elif (m % 2) ==0 :
            return ( np.sin(alpha/2)*np.cos(beta/2) /  ( ( (alpha/(m*np.pi))**2-1)*((beta/(n*np.pi))**2-1) ))**2  * np.sin(theta)     
        elif (n % 2) ==0 :
            return ( np.cos(alpha/2)*np.sin(beta/2) /  ( ( (alpha/(m*np.pi))**2-1)*((beta/(n*np.pi))**2-1) ))**2  * np.sin(theta)
        else :
            return ( np.cos(alpha/2)*np.cos(beta/2) /  ( ( (alpha/(m*np.pi))**2-1)*((beta/(n*np.pi))**2-1) ))**2  * np.sin(theta)
    
    Nf = len(k0)
    integral = np.zeros((Nf,))
    
    for iif in range(Nf):
        print('Calculating %.f of %.f'%(iif, Nf))
        k= k0[iif]
        integral[iif] = dblquad(wallace,  t1, t2,
                                lambda theta:   p1, lambda theta:   p2)[0]
    
    return 64*k0**2*Lx*Ly / (np.pi**6 * m**2 * n**2) * integral

# Self radiation efficiency in pistons
def sigma_11(ka):  # Beranek ApxII, eq1
    # ka = k*a
    return (1 - sp.jv(1, 2*ka)/ka)

# Mutual radiation efficiency in pistons
def sigma_12a(k, a, d):
    ka = k*a
    kd = k*d

    # Mutual Radiation Beranek
    IM, IN = 10,10
    try:
        L = len(ka)
    except:
        L = len(kd)
     
    s12_mn = np.zeros((L, IM, IN))
    
    for m in range(IM):
        for n in range(IN):            
            s12_mn[:, m,n] = (ka/(kd))**(m+n)   \
                            *np.sqrt(2/(kd))    \
                            *sp.gamma(m+n+1/2)  \
                            *sp.jv(m+1, ka)     \
                            *sp.jv(n+1, ka)     \
                            *sp.jv(m+n+1/2, kd) \
                            /(sp.factorial(m)*sp.factorial(n))

    return np.sum(np.sum(s12_mn, 1),1)

def sigma_12(k, a, d):
    # paper carlos
    ka = k*a
    kd = k*d
    f1 =   sp.jv(1,ka)**2
    f2 =  (a/d)*sp.jv(1,ka)*sp.jv(2,ka)*(1/kd - 1/np.tan(kd)) 
    f3 = 3/4*(a/d)**2*sp.jv(2,ka)**2*( 3/(kd**2) - 3/(kd*np.tan(kd)) -1 )
    return 2*np.sin(kd)/kd *(f1+f2+f3)

def Rmn(m,n, Lx, Ly):
# def Rmn(m,n,d):
    cc = [0]*m*n
    ii=0
    for ix in range(0,m):
        for iy in range(0,n):
            cc[ii] = (ix,iy)
            ii=ii+1
    
    # create coords Ids
    coords = {}
    for ii in range(len(cc)):
        if cc[ii][0] == cc[ii][1]:
            coords[ii] = [cc[ii], '+']
        elif (cc[ii][0]+cc[ii][1]) % 2 ==0:
            coords[ii] = [cc[ii], '+']
            if cc[ii][1] == 0 and cc[ii][0] %2 != 0:
                coords[ii] = [cc[ii], '-']
        else:
            coords[ii] = [cc[ii], '-']
    
    # Index of +- coordinates
    idx_plus  = np.where([coords[i][1]=='+' for i in range(len(cc))])[0]
    idx_minus = np.where([coords[i][1]=='-' for i in range(len(cc))])[0]
    
    Nplus  = len(idx_plus)
    Nminus = len(idx_minus)
    
    r_pp = np.zeros((Nplus, Nplus))
    r_mm = np.zeros((Nminus, Nminus))
    r_pm = np.zeros((Nminus, Nplus))
    
    dx, dy = Lx/m, Ly/n
    plus  = np.float32(np.array([coords[i][0] for i in idx_plus]) )
    minus = np.float32(np.array([coords[i][0] for i in idx_minus]))
    
    pm = [plus, minus]
    
    if m>1 or n>1:
        plus[:,0] = plus[:,0]*dx
        plus[:,1] = plus[:,1]*dy
        minus[:,0] = minus[:,0]*dx
        minus[:,1] = minus[:,1]*dy
    
    def calculateDistance(array, point):
        # x is a vector, y is a point
        x1 = array[:,0]
        y1 = array[:,1]
        x2 = point[0]
        y2 = point[1]
        
        dist = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
        return dist
    
    for iim in range(len(plus)):
        r_pp[iim,:] = calculateDistance(plus, plus[iim])
    
    for iim in range(len(minus)):
        r_mm[iim,:] = calculateDistance(minus, minus[iim])
    
    for iim in range(len(minus)):
        r_pm[iim,:] = calculateDistance(plus, minus[iim])
    
    # for ii in range
    R_pp = r_pp[r_pp != 0.]
    R_mm = r_mm[r_mm != 0.]
    R_pm = r_pm[r_pm != 0.]
    
    return R_pp, R_mm, R_pm, pm

def sigma_mnp(m,n,k0,Lx,Ly):
    a0 = np.sqrt((Lx/m)*(Ly/n)/np.pi)            # equivalent radius for same AREA (hypothesis, but simmetrization?)
    perimeter_circle = 2*np.pi*a0
    perimeter_complex   = 2*(Lx/m + Ly/n)

    ratio = perimeter_circle/perimeter_complex     # correction
    print(ratio)

    # Plots
    coef_piston_ss = 0.87
    a = a0*coef_piston_ss*ratio**(3/5)

    R_pp, R_mm, R_pm, coords_pm = Rmn(m,n, Lx, Ly)

    mono      = (m*n)*sigma_11(k0*a).reshape(-1)    
    dipole_plus  = np.sum([sigma_12(k0,a,di) for di in R_pp], 0)
    dipole_minus = np.sum([sigma_12(k0,a,di) for di in R_mm], 0)
    dipole_plus_minus = 2*np.sum([sigma_12(k0,a,di) for di in R_pm], 0)

    multipole = (mono + dipole_plus + dipole_minus - dipole_plus_minus)/(m * n)
    return multipole#, coords_pm, a

# Radition efficiency simply supported radiator Greenspan
def Zs(ka):
    y = 2*ka
    F1s = (10-y**2)*sp.jv(1,y) - 5*y*sp.jv(0,y) - y**3/8
    return  1 - 96/((2*ka)**5) * F1s  