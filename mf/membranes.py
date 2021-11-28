from dolfin import *
import numpy as np

def EigenSolveMembrane(N, u, v, T, rho_s, bcs_w):
    c2 = T/rho_s
    a = inner(grad(u), grad(v))*dx 
    b = 1/c2*inner(u,v)*dx
    
    M_ = PETScMatrix()
    assemble(L, tensor = M_)
    
    # Assemble matrices
	K_ = PETScMatrix()
	assemble(a, tensor = K_)

	B = PETScMatrix()
	assemble(b, tensor = M_)

	# Apply boundary conditions
	bc.apply(A)

    # object is a list or not 
    if isinstance(bcs_w, list): 
        bcs_w  = bcs_w
    else:        
        bcs_w  = [bcs_w]
        
    for bcs in bcs_w:
        bcs.apply(K_)
        bcs.zero(M_)                                            # avoid spurius eig vals
    
    tt1 = time.time()
    print('Computing %.f Modes' %N)           
    solver = SLEPcEigenSolver(K_,M_)                            #[𝐾]{𝑈} = 𝜆[𝑀]{𝑈}
    solver.parameters["solver"]             = "krylov-schur"
    solver.parameters["problem_type"]       = "gen_hermitian"
    solver.parameters['spectral_transform'] = 'shift-and-invert'
    solver.parameters['spectral_shift']     = 1e-14
    solver.solve(N)
    k = solver.get_number_converged()
    print(' %.f Modes have converged in %.2f secs!' %(k, time.time() - tt1))
    
    return solver