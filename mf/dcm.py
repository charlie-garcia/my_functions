import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as inter
# from fenics import *
import warnings
#ignore by message
warnings.filterwarnings("ignore", message="invalid value encountered in true_divide")

def MatrixR(rho,c0,k0,A,rij):
    import scipy.special as sp
    import numpy as np
    #Compute the Acoustical Resistance Radiation Matrix via Hashimoto (discrete calculation method)
    ai  = np.sqrt(A/np.pi)
    ak  = np.sqrt(A.T/np.pi)
    AA  = A*A.T
    
    R = 2*rho*c0*(k0**2)*AA /(np.pi)* ( sp.jv(1, (k0*ai) )/ (k0*ai))* ( sp.jv(1, (k0*ak) )/ (k0*ak))*  np.sin(k0*rij)/(k0*rij)
    Rii = rho*c0*A* ( 1 - sp.jv(1, (2*k0*ai) ) /(k0*ai) )
    R[np.isnan(R)] = 0
    np.fill_diagonal(R, Rii)
    R = np.real(R)
    return R

def MatrixR32(rho0,c0,k0,A,rij):
    import scipy.special as sp
    import numpy as np
    #Compute the Acoustical Resistance Radiation Matrix via Hashimoto (discrete calculation method)
    ai  = np.float32( np.sqrt(A/np.pi))
    ak  = np.float32( np.sqrt(A.T/np.pi) )
    A = np.float32(A)
    AA  = A*A.T
    
    R = 2*rho*c0*(k0**2)*AA /(np.pi)* ( sp.jv(1, (k0*ai) )/ (k0*ai))* ( sp.jv(1, (k0*ak) )/ (k0*ak))*  np.sin(k0*rij)/(k0*rij)
    Rii = rho*c0*A* ( 1 - sp.jv(1, (2*k0*ai) ) /(k0*ai) )
    R[np.isnan(R)] = 0
    np.fill_diagonal(R, Rii)
    R = np.real(R)
    return R
    
    return R

def MatrixRVEC(rho0,c0,k0,A,rij):
    import scipy.special as sp
    import numpy as np
    M = len(rij);   Nf = len(k0)
    #Compute the Acoustical Resistance Radiation Matrix via Hashimoto (discrete calculation method)
    ai  = np.sqrt(A/np.pi).reshape(M,1,1)
    ak  = np.sqrt(A.T/np.pi).reshape(1,M,1)
    AA  = (A*A.T).reshape(M,M,1)
    k0 = k0.reshape(1,1,Nf)
    rij = rij.reshape(M,M,1)
    
    R = 2*rho0*c0*(k0**2)*AA / (np.pi)  *\
        ( sp.jv(1, (k0*ai) ) / (k0*ai)) *\
        ( sp.jv(1, (k0*ak) ) / (k0*ak)) *\
          np.sin(k0*rij)/(k0*rij)  
    
    Rii = rho0*c0*A* ( 1 - sp.jv(1, (2*k0.reshape(1,Nf)*ai.reshape(M,1)) ) /(k0.reshape(1,Nf)*ai.reshape(M,1)) )      # Get only array
    R[np.isnan(R)]     = 0
    Rii[np.isnan(Rii)] = 0
    
    np.einsum('iij->ij', R)[...] = Rii
    
    R = np.real(R)
    return R

def MatrixRVEC32(rho0,c0,k0,A,rij):
    import scipy.special as sp
    import numpy as np
    M = len(rij);   Nf = len(k0)
    #Compute the Acoustical Resistance Radiation Matrix via Hashimoto (discrete calculation method)
    ai  = np.sqrt(A/np.pi, dtype=np.float32).reshape(M,1,1)
    ak  = np.sqrt(A.T/np.pi, dtype=np.float32).reshape(1,M,1)
    
    A = np.float32(A)
    AA  = (A*A.T).reshape(M,M,1)
    
    k0 = np.float32(k0).reshape(1,1,Nf)
    rij = np.float32(rij).reshape(M,M,1)

    R = 2*rho0*c0*(k0**2)*AA / (np.pi)  *\
        ( sp.jv(1, (k0*ai) ) / (k0*ai)) *\
        ( sp.jv(1, (k0*ak) ) / (k0*ak)) *\
          np.sin(k0*rij)/(k0*rij)  
    
    Rii = rho0*c0*A* ( 1 - sp.jv(1, (2*k0.reshape(1,Nf)*ai.reshape(M,1)) ) /(k0.reshape(1,Nf)*ai.reshape(M,1)) )      # Get only array
    R[np.isnan(R)]     = 0
    Rii[np.isnan(Rii)] = 0
    
    np.einsum('iij->ij', R)[...] = Rii
    
    R = np.real(R)
    return R


def MatrixRHDF5(rho0,c0,k0,A,rij, filename):
    import scipy.special as sp
    import numpy as np
    #Compute the Acoustical Resistance Radiation Matrix via Hashimoto (discrete calculation method)
    ai  = np.sqrt(A/np.pi)
    ak  = np.sqrt(A.T/np.pi)
    AA  = A*A.T
    lk = len(k0)
    
    import tables as tb
    l, m, n = rij.shape[0], rij.shape[1], lk
    file = tb.open_file(filename, 'w')
    filters = tb.Filters(complevel=3, complib='blosc:lz4')
    # filters = tb.Filters(complevel=5)
    out = file.create_carray(file.root, 'data', tb.Float32Atom(), shape=(l, m, n), filters=filters)
    
    for iif in range(0,len(k0)):
        print('Compution ARM, freq %.f from %.f calculated' %( iif, lk ) )
        Rd = 2*rho0*c0*(k0[iif]**2)*AA /(np.pi)* ( sp.jv(1, (k0[iif]*ai) )/ (k0[iif]*ai))* ( sp.jv(1, (k0[iif]*ak) )/ (k0[iif]*ak))*  np.sin(k0[iif]*rij)/(k0[iif]*rij)

        Rii = rho0*c0*A* ( 1 - sp.jv(1, (2*k0[iif]*ai) ) /(k0[iif]*ai) )
        Rd[np.isnan(Rd)] = 0
        np.fill_diagonal(Rd, Rii)
        
        out[:,:,iif]  = np.real(Rd)
    
    file.close()

def MatrixRHDF5Vect(rho0,c0,k0,A,rij, filename):
    import scipy.special as sp
    import numpy as np
    M = len(rij);   Nf = len(k0)
    #Compute the Acoustical Resistance Radiation Matrix via Hashimoto (discrete calculation method)
    ai  = np.sqrt(A/np.pi).reshape(M,1,1)
    ak  = np.sqrt(A.T/np.pi).reshape(1,M,1)
    AA  = (A.reshape(M,1)*A.reshape(1,M)).reshape(M,M,1)
    k0 = k0.reshape(1,1,Nf)
    rij = rij.reshape(M,M,1)
    
    import tables as tb
    l, m, n = rij.shape[0], rij.shape[1], Nf
    file = tb.open_file(filename, 'w')
    filters = tb.Filters(complevel=3, complib='blosc:lz4')
    # filters = tb.Filters(complevel=5)
    out = file.create_carray(file.root, 'data', tb.Float32Atom(), shape=(l, m, n), filters=filters)
    
    R = 2*rho0*c0*(k0**2)*AA / (np.pi) * ( sp.jv(1, (k0*ai) ) / (k0*ai)) * ( sp.jv(1, (k0*ak) ) / (k0*ak)) * np.sin(k0*rij)/(k0*rij)  
    Rii = rho0*c0*A* ( 1 - sp.jv(1, (2*k0.reshape(1,Nf)*ai.reshape(M,1)) ) /(k0.reshape(1,Nf)*ai.reshape(M,1)) )      # Get only array
    R[np.isnan(R)]     = 0
    Rii[np.isnan(Rii)] = 0
    
    np.einsum('iij->ij', R)[...] = Rii
    # R = np.real(R)
    
    for iif in range(Nf):
        print('Saving ARM, freq %.f from %.f ' %( iif, Nf ) )
        out[:,:,iif]  = np.real(R)

    file.close()
    
def MatrixRSparse(rho0,c0,k0,area,rij):
    import scipy.special as sp
    from scipy.sparse import csr_matrix, triu, eye, diags
    import numpy as np
    #Compute the Acoustical Resistance Radiation Matrix via Hashimoto (discrete calculation method)
    ai  = np.sqrt(area/np.pi).reshape( len(area), 1)
    ak  = np.sqrt(area/np.pi).reshape( 1, len(area) )
    A = csr_matrix(area)
    AA  = np.dot(A, A.T)
    
    Rii,R_diag = ([], [])
    for iif in range(0,len(k0)):
        print('freq %.f from %.f calculated' %( iif,len(k0) ))
        # R0 = 2*rho0*c0*(k0[iif]**2)*AA /(np.pi)
        # R1 = csr_matrix( sp.jv(1, (k0[iif]*ai) )/ (k0[iif]*ai))
        # R2 = csr_matrix( sp.jv(1, (k0[iif]*ak) )/ (k0[iif]*ak))
        # R3 = csr_matrix( np.sin(k0[iif]*(rij+rij.T))/(k0[iif]*(rij+rij.T)))
        
        # Rii =  eye(len(area)).multiply(rho0*c0*A.multiply(csr_matrix( 1 - sp.jv(1, (2*k0[iif]*ai) ) /(k0[iif]*ai) ))) 
        
        # # R_diag = triu((R0.multiply(R1*R2)).multiply(R3),k=1)
        
        # R_diag = triu(((2*rho0*c0*(k0[iif]**2)*AA /(np.pi)).multiply(csr_matrix( sp.jv(1, (k0[iif]*ai) )/ (k0[iif]*ai))\
        #                          *csr_matrix( sp.jv(1, (k0[iif]*ak) )/ (k0[iif]*ak)))).multiply(csr_matrix( np.sin(k0[iif]*(rij+rij.T))/(k0[iif]*(rij+rij.T)))),k=1)
        # R.append(triu((2*rho0*c0*(k0[iif]**2)*AA /(np.pi)).multiply((csr_matrix( sp.jv(1, (k0[iif]*ai) )/ (k0[iif]*ai)))\
        # *(csr_matrix( sp.jv(1, (k0[iif]*ak) )/ (k0[iif]*ak)))).multiply(csr_matrix( np.sin(k0[iif]*(rij+rij.T))/(k0[iif]*(rij+rij.T)))),k=1))
        
        # R.append(R_diag+R_diag.T+Rii)
        Rii.append(eye(len(area)).multiply(rho0*c0*A.multiply(csr_matrix( 1 - sp.jv(1, (2*k0[iif]*ai) ) /(k0[iif]*ai) ))) )
        
        # R_diag = triu((R0.multiply(R1*R2)).multiply(R3),k=1)
        
        R_diag.append( triu(((2*rho0*c0*(k0[iif]**2)*AA /(np.pi)).multiply(csr_matrix( sp.jv(1, (k0[iif]*ai) )/ (k0[iif]*ai))\
                                 *csr_matrix( sp.jv(1, (k0[iif]*ak) )/ (k0[iif]*ak)))).multiply(csr_matrix( np.sin(k0[iif]*(rij+rij.T))/(k0[iif]*(rij+rij.T)))),k=1) )
    return Rii, R_diag

def VectorR(rho,c0,k0,A,rij):
    import scipy.special as sp
    import numpy as np
    #Compute the Acoustical Resistance Radiation Vector
    ai  = np.sqrt(A/np.pi)
    Rii = rho*c0*A* ( 1 - sp.jv(1, (2*k0*ai) ) /(k0*ai) )

    return Rii

def distance3(cs,FieldPoint):
    r = np.zeros([cs.shape[0], FieldPoint.shape[0]], dtype=np.float32)
    for ii in range(FieldPoint.shape[0]):
        r[:,ii] = np.sqrt( np.sum( (cs-FieldPoint[ii,:])**2 ,1) )
    return r

def distance3Sparse(cs,FieldPoint):
    from scipy.sparse import triu
    r = []
    for ii in range(FieldPoint.shape[0]):
        r0 = np.sqrt( np.sum( (cs-FieldPoint[ii,:])**2 ,1) )
        r.append(triu(r0, k=1).toarray())
    rr = triu(np.array(r).reshape(len(cs), len(cs)), k=1)
    return rr

def RadModes(R, n_f, n_modes):
    import scipy.linalg as la            

    N = int(np.sqrt( R.shape[0] ))
    
    ll, Q = la.eig(R)                                          # RADIATION MODES 
    mode_i = [0]*n_modes
    efficiency_i = [0]*n_modes
    
    for irm in range(0, n_modes):
        # mQ[irm] = Q.T[irm,:]
        mode_i[irm] = np.real( Q.T[irm,:] ).reshape(N,N) 
        efficiency_i[irm] = ll[irm]
    return mode_i, efficiency_i

