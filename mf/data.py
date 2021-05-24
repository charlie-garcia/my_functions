import numpy as np
import time

# Load Modes
def LoadComplexPlateModes(Nmode, Nfx_points):
    # Initialize table to save computations
    import tables as tb
    path = '/home/carlos/dataPython/ComplexPlate/' 
    filename_modes = path + str(Nmode)+'_modes_'+str(Nfx_points)+'_fixed_points.h5'
    h5 = tb.open_file(filename_modes, 'r')
    # Initialize arrays 
    mode_p   = np.zeros((h5.root.modes_ctr.shape))
    omega_p  = np.zeros((h5.root.omega_p.shape))
    mass_p   = np.zeros((h5.root.mass_p.shape))
    phi_x0y0 = np.zeros((h5.root.phi_x0y0.shape))

    tt1 = time.time()
    print('Loading %.f Modes...' %Nmode)
    for jj in range(Nmode):
        mode_p[:,jj]   = h5.root.modes_ctr[:,jj]
        omega_p[0,jj]  = h5.root.omega_p[0,jj]
        mass_p[0,jj]   = h5.root.mass_p[0,jj]
        phi_x0y0[:,jj] = h5.root.phi_x0y0[:,jj]

    h5.close()
    print('Loaded in : %.2f seconds\n' %(time.time() - tt1))
    return mode_p , omega_p, mass_p, phi_x0y0
    