import gmsh, sys, meshio 
import numpy as np
from mf.fem import gmsh2dolfin, gmsh2dolfin_subd
from fenics import *
from matplotlib import pyplot as plt

def CreatePlate(mesh_name, Lx, Ly, pxy, h1, dim, plot_info):
    # Mesh generation with GMSH
    gmsh.initialize(sys.argv)
    my_path = 'meshes/'
    # Ask GMSH to display information in the terminal
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add(mesh_name)
    model = gmsh.model
    model.add("MyPlate")
    
    # We start by defining some points and some lines
    p,l = ([0]*4 for line in range(2))
    
    p[0] = model.geo.addPoint(0,     0,   0,  h1);
    p[1] = model.geo.addPoint(Lx,    0,   0,  h1);
    p[2] = model.geo.addPoint(Lx,    Ly,  0,  h1);
    p[3] = model.geo.addPoint(0,     Ly,  0,  h1);
    
    l[0] = model.geo.addLine(1, 2);
    l[1] = model.geo.addLine(2, 3);
    l[2] = model.geo.addLine(3, 4);
    l[3] = model.geo.addLine(4, 1);

    p0 = model.geo.addPoint(pxy[0], pxy[1], 0, h1);
    p1 = model.geo.addPoint(pxy[0]+Lx**2/300, pxy[1], 0, h1);
    l0 = model.geo.addLine(p0, p1)
    gmsh.model.geo.mesh.setTransfiniteCurve(l0, int(np.round(1/(2*h1))))  # call later the mesh.embed command to work    
    
    curve = model.geo.addCurveLoop( l )
    surface = model.geo.addPlaneSurface([curve])
    
    model.geo.synchronize()
    
    model.mesh.embed(0, [p0], 1, l0)                                            # dim =1, line 
    model.mesh.embed(0, [p1], 1, l0)                                            # dim =1, line 
    model.mesh.embed(1, [l0], 2, surface)                                       # dim =1, line    
    
    
    # Mesh (2D)
    model.mesh.generate(2)
    
    line_string_tag = "my_point"
    tag_inside = gmsh.model.addPhysicalGroup(1, [p0])                           # REMEMBER THIS TAG FOR BC
    gmsh.model.setPhysicalName(0, tag_inside, line_string_tag)                  # dim, tag, name

    bord_string_tag = "my_bords"
    tag_bc = [0]*4
    for ii in range(0,4):
        tag_bc[ii] = gmsh.model.addPhysicalGroup(1, [l[ii]])                    # REMEMBER THIS TAG FOR BC
        gmsh.model.setPhysicalName(1, tag_bc[ii], bord_string_tag)              # dim, tag, name
    
    # # # surface phisical group
    surface_string_tag = "my_surface"
    tag_dom= model.geo.addPhysicalGroup(2, [surface])
    gmsh.model.setPhysicalName(2, tag_dom, surface_string_tag)                  # dim, tag, name
    
    model.geo.synchronize()
    # # =============================================================================
    # # the command line arguments:
    if plot_info == 'plot':
        gmsh.fltk.run()
    # # =============================================================================
    
    gmsh.write(my_path+mesh_name)
    gmsh.finalize()
    
    # From gmsh to dolfin
    fmesh, mf = gmsh2dolfin(my_path, mesh_name, dim, bord_string_tag, surface_string_tag)

    return fmesh, mf, tag_bc, tag_inside

def TestPlate(my_path, mesh_name, Lx, Ly, h1, dim, plot_info):
    # Mesh generation with GMSH
    gmsh.initialize(sys.argv)
    # Ask GMSH to display information in the terminal
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add(mesh_name)
    
    model = gmsh.model
    model.add("MyPlate")
    
    # We start by defining some points and some lines
    p,l = ([0]*4 for line in range(2))
    
    p[0] = model.geo.addPoint(0,     0,   0,  h1);
    p[1] = model.geo.addPoint(Lx,    0,   0,  h1);
    p[2] = model.geo.addPoint(Lx,    Ly,  0,  h1);
    p[3] = model.geo.addPoint(0,     Ly,  0,  h1);
    
    l[0] = model.geo.addLine(1, 2);
    l[1] = model.geo.addLine(2, 3);
    l[2] = model.geo.addLine(3, 4);
    l[3] = model.geo.addLine(4, 1);
    
    curve = model.geo.addCurveLoop( l )
    surface = model.geo.addPlaneSurface([curve])
    
    model.geo.synchronize()
    
    
    # Mesh (2D)
    model.mesh.generate(2)

    bord_string_tag = "my_bords"
    tag_bc = [0]*4
    for ii in range(0,4):
        tag_bc[ii] = gmsh.model.addPhysicalGroup(1, [l[ii]])                    # REMEMBER THIS TAG FOR BC
        gmsh.model.setPhysicalName(1, tag_bc[ii], bord_string_tag)              # dim, tag, name
    
    # # # surface phisical group
    surface_string_tag = "my_surface"
    tag_dom= model.geo.addPhysicalGroup(2, [surface])
    gmsh.model.setPhysicalName(2, tag_dom, surface_string_tag)                  # dim, tag, name
    
    model.geo.synchronize()
    # # =============================================================================
    # # the command line arguments:
    if plot_info == 'plot':
        gmsh.fltk.run()
    # # =============================================================================
    
    gmsh.write(my_path+mesh_name)
    gmsh.finalize()
    
    # From gmsh to dolfin
    fmesh, mf = gmsh2dolfin(my_path, mesh_name, dim, bord_string_tag, surface_string_tag)

    return fmesh, mf, tag_bc

def TestComplexPlate(my_path, mesh_name, Lx, Ly, h1, dim, line_string_tag, surface_string_tag, plot_info):
        
    # Mesh generation with GMSH
    gmsh.initialize(sys.argv)
    
    # Ask GMSH to display information in the terminal
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add(mesh_name)
    
    model = gmsh.model
    model.add("MyPlate")
    # parameters from manip
    Lx = 0.5
    Ly = 0.6
    addRectangle = gmsh.model.occ.addRectangle
    
    rectangle = addRectangle(0, 0, 0, Lx, Ly)
    obstacle = addRectangle(Lx/3, 0.1*Ly, 0, Lx*0.05, 0.8*Ly)
    gmsh.model.occ.fragment([(2,rectangle)], [(2, obstacle)] )                  # Fusion 2 surfaces, no tag needed
    
    gmsh.model.occ.synchronize()
    
    # Mesh (2D)
    gmsh.model.mesh.setSize(gmsh.model.getEntities(0), h1)
    gmsh.model.mesh.generate(2)
    #% Add physical groups
    # Bord phisical group
    
    tag_border = gmsh.model.addPhysicalGroup(1, [1,2,3,4])                # rectangle is the tag of the rectangle, not of the lines!
    gmsh.model.setPhysicalName(1, tag_border, line_string_tag)              # dim, tag, name
    
    # Stiffener phisical group
    tag_stiff_bord  = gmsh.model.addPhysicalGroup(1, [5,6,7,8])                # REMEMBER THIS TAG FOR BC
    gmsh.model.setPhysicalName(1, tag_stiff_bord, line_string_tag)          # dim, tag, name
    
    # Surface phisical group
    # surface_string_tag = "surface"
    tag_dom = gmsh.model.addPhysicalGroup(2, [rectangle+2])      # Delete original tag when fragment
    gmsh.model.setPhysicalName(2, tag_dom, surface_string_tag)              # dim, tag, name
    
    tag_stiff = gmsh.model.addPhysicalGroup(2, [obstacle])      # Delete original tag when fragment
    gmsh.model.setPhysicalName(2, tag_stiff, surface_string_tag)              # dim, tag, name
    
    # # =============================================================================
    # # the command line arguments:
    if plot_info == 'plot':
        gmsh.fltk.run()
    # # =============================================================================
    
    gmsh.write(my_path+mesh_name)
    gmsh.finalize()
    
    return tag_border, tag_stiff_bord, tag_dom, tag_stiff


def GetRandomPointExcitation(Npts, Lx, Ly, d2boundary, d2source):
    from matplotlib import pyplot as plt
    from numpy import random as npr

    d2source = 2*d2boundary
    print('Minimal distance: %.2f' %d2boundary)
    
    px = []
    py = []
    circles = []
    
    for ii in range(Npts):
        px.append(d2boundary + (Lx - 2*d2boundary)*npr.rand())
        py.append(d2boundary + (Ly - 2*d2boundary)*npr.rand())
        dist = np.sqrt( (px[ii] - px[:-1])**2 + (py[ii] - py[:-1])**2 )
        if ii>=1:
            while any(t < d2source for t in dist):
                # print(dist)
                px[ii] = d2boundary + (Lx - 2*d2boundary)*npr.rand()
                py[ii] = d2boundary + (Ly - 2*d2boundary)*npr.rand()
                dist = np.sqrt( (px[ii] - px[:-1])**2 + (py[ii] - py[:-1])**2 )
        
    	#Draw a circle to see if min_dist is satisfied
        circles.append( plt.Circle((px[ii], py[ii]), d2boundary, color='grey', alpha=0.3) )
    return px, py, circles


def ComplexPate(my_path, mesh_name, Lx, Ly, h1, loc_x, loc_y, plot_info):
    import gmsh, sys
    from mf.fem import write_gmsh
    from mf.fem import gmsh2dolfin
    gmsh.initialize(sys.argv)
    gmsh.option.setNumber("General.Terminal", 1)        # Ask GMSH to display information in the terminal
    
    model = gmsh.model
    model.add("complex_plate")
    
    # Create Rectangle
    dim    = '2D'                              
    a_mass = 0.003                         # blocked surface radius
    
    # Create Polygon
    N      = 8
    alpha  = np.deg2rad(360/N)
    phi    = alpha/2+np.pi/2
    beta   = np.deg2rad(90) - alpha/2
    l2     = a_mass*np.sin(alpha/2)/np.sin(beta)
    
    # Create plate
    rect         = gmsh.model.occ.addRectangle(0,0,0, Lx, Ly)
    border_plate = [1,2,3,4]
    
    # all added circular of added mass
    Nfx = len(loc_x)
    c      =  [0]*Nfx
    vertex, borders = ( [] for i in range(2))
    
    for ii in range(Nfx):
        c[ii] = gmsh.model.occ.addPoint(loc_x[ii], loc_y[ii], 0, h1)
        
        vx, bs = ( [] for i in range(2))
        for jj in range(N):
            vx.append(model.occ.addPoint( a_mass*np.cos(2*np.pi*jj/N - phi) + loc_x[ii],\
                                          a_mass*np.sin(2*np.pi*jj/N - phi) + loc_y[ii], 0, h1))
        # connect borders
        for j in range(N):
            if j < N:
                bs.append( gmsh.model.occ.addLine(vx[j-1], vx[j]))
        
        # connect inside triangles            
        for j in range(N):
            bs.append( gmsh.model.occ.addLine(vx[j], c[ii]))
    
        vertex = vertex+vx
        borders = borders+bs
    
    gmsh.model.occ.synchronize()
    
    gmsh.model.mesh.embed(0, vertex,  2, rect)             # dim =0, point in dim=2, surface
    gmsh.model.mesh.embed(1, borders, 2, rect)             # dim =1, curve in dim=2, surface  
    
    # Append conectric points to redefine mesh
    ixs = [3]
    vx_out  = []
    for ii in range(Nfx):
        for jj in range(N):
            for ix in ixs:
                vx_out.append(model.occ.addPoint( ix*a_mass*np.cos(2*np.pi*jj/N - phi) + loc_x[ii],\
                                                  ix*a_mass*np.sin(2*np.pi*jj/N - phi) + loc_y[ii], 0, h1))
    
    gmsh.model.occ.synchronize()
    gmsh.model.mesh.embed(0, vx_out,  2, rect)              # dim =0, point in dim=2, surface
    
    gmsh.model.mesh.setAlgorithm(2, rect, 8)             # algorithm "Packing of parallelograms" (experimental =9)
    gmsh.model.mesh.setSize(gmsh.model.getEntities(0), h1)
    # gmsh.option.setNumber('Mesh.MeshSizeMin', h1)
    gmsh.model.mesh.generate(2)                                 # 2D mesh
    
    # Bord phisical group
    line_string_tag = "border"                                  # all borders even inside
    tag_border = gmsh.model.addPhysicalGroup(1, border_plate+borders)
    gmsh.model.setPhysicalName(1, tag_border, line_string_tag)   # dim, tag, name
    
    # Ponint zones phisical group
    # tag_points = [0]*(Nfx)
    # len0=2*N
    # for ii in range(Nfx):
    #     tag_points[ii] = gmsh.model.addPhysicalGroup(1, borders[len0*(ii):len0*(ii+1)])    
    #     gmsh.model.setPhysicalName(1, tag_points[ii], line_string_tag)   # dim, tag, name
    
    # Surface phisical group
    surface_string_tag = "surface"
    tag_dom = gmsh.model.addPhysicalGroup(2, [rect])       # Delete original tag when fragment
    gmsh.model.setPhysicalName(2, tag_dom, surface_string_tag)   # dim, tag, name
    
    if plot_info == 'plot':
        gmsh.fltk.run()

    write_gmsh(my_path, mesh_name)
    
    gmsh.finalize()
    
    fmesh, mf_boundary = gmsh2dolfin(my_path, mesh_name, dim, line_string_tag, surface_string_tag)
    
    return fmesh, mf_boundary, tag_border

def Clamped_Plate(W, w_, mesh, E_, nu_, t_, ds, tag_bords):
    # We take constant material properties throughout the domain::
    E = Constant(E_)
    nu = Constant(nu_)
    t = Constant(t_)
    
    # Rotation \theta = \nabla w, which can be expressed in UFL as::
    theta = grad(w_)
    
    # The bending tensor can then be calculated from the derived rotation field
    #     k = \frac{1}{2}(\nabla \theta + (\nabla \theta)^T) 
    k = variable(sym(grad(theta)))
    
    # Again, identically to the Reissner-Mindlin model we can calculate the bending
    # energy density as::
    D = (E*t**3)/(12.0*(1.0 - nu**2))
    psi_M = 0.5*D*((1.0 - nu)*tr(k*k) + nu*(tr(k))**2)
    
    # Clamped Plate
    # For the definition of the CDG stabilisation terms and the (weak) enforcement of
    # the Dirichlet boundary conditions on the rotation field, we need to explicitly
    # derive the moment tensor :math:`M`. Following standard arguments in elasticity,
    #     M = \frac{\partial \psi_M}{\partial k}
    M = diff(psi_M, k)
    
    # The Lagrangian formulation of the CDG stabilisation term is then:
    #     L_{\mathrm{CDG}}(w) = \sum_{E \in \mathcal{E}_h^{\mathrm{int}}} \int_{E} - [\!\![ \theta ]\!\!]  \cdot \left< M \cdot (n \otimes n) \right > + \frac{1}{2} \frac{\alpha}{\left< h_E \right>} \left< \theta \cdot n \right> \cdot \left< \theta \cdot n \right> \; \mathrm{d}s
    
    # We choose the penalty parameter to be on the order of the norm of the bending stiffness matrix :math:`\dfrac{Et^3}{12}`.
    
    alpha = E*t**3
    h = CellDiameter(mesh)
    h_avg = (h('+') + h('-'))/2.0 
    
    n = FacetNormal(mesh)
    
    M_n = inner(M, outer(n, n))
    
    L_CDG = -inner(jump(theta, n), avg(M_n))*dS + \
                (1.0/2.0)*(alpha('+')/h_avg)*inner(jump(theta, n), jump(theta, n))*dS
    

    # In this example, we would like :math:`\theta_d = 0` everywhere on the boundary:: 
     
    theta_d = Constant((0.0, 0.0))                                                  # To impose non-homogeneous BC
    # The definition of the exterior facets and Dirichlet rotation field were trivial
    # in this demo, but you could extend this code straightforwardly to
    # non-homogeneous Dirichlet conditions.
    # 
    # The weak boundary condition enforcement term can be written:
    #
    # .. math::
    #     L_{\mathrm{BC}}(w) = \sum_{E \in \mathcal{E}_h^{\mathrm{D}}} \int_{E} - \theta_e  \cdot (M \cdot (n \otimes n))  + \frac{1}{2} \frac{\alpha}{h_E} (\theta_e \cdot n)  \cdot (\theta_e \cdot n)  \; \mathrm{d}s
    # 
    # where :math:`\theta_e = \theta - \theta_d` is the effective rotation field, and
    # :math:`\mathcal{E}_h^{\mathrm{D}}` is the set of all exterior facets of the triangulation
    # :math:`\mathcal{T}` where we would like to apply Dirichlet boundary conditions, or in UFL::
    
    theta_effective = theta - theta_d 
    L_BC = -inner(inner(theta_effective, n), M_n)*ds(tag_bords) + \
            (1.0/2.0)*(alpha/h)*inner(inner(theta_effective, n), inner(theta_effective, n))*ds(tag_bords) 
    
    # The remainder of the demo is as usual::
    L = psi_M*dx + L_CDG + L_BC   
    
    return L

def SS_Plate(W, w_, mesh, E_, nu_, t_):
   
    # We take constant material properties throughout the domain::
    E = Constant(E_)
    nu = Constant(nu_)
    t = Constant(t_)
    
    # Rotation \theta = \nabla w, which can be expressed in UFL as::
    theta = grad(w_)
    
    # The bending tensor can then be calculated from the derived rotation field
    #     k = \frac{1}{2}(\nabla \theta + (\nabla \theta)^T) 
    k = variable(sym(grad(theta)))
    
    # Again, identically to the Reissner-Mindlin model we can calculate the bending
    # energy density as::
    D = (E*t**3)/(12.0*(1.0 - nu**2))
    psi_M = 0.5*D*((1.0 - nu)*tr(k*k) + nu*(tr(k))**2)
    
    # Clamped Plate
    # For the definition of the CDG stabilisation terms and the (weak) enforcement of
    # the Dirichlet boundary conditions on the rotation field, we need to explicitly
    # derive the moment tensor :math:`M`. Following standard arguments in elasticity,
    #     M = \frac{\partial \psi_M}{\partial k}
    M = diff(psi_M, k)
    
    # The Lagrangian formulation of the CDG stabilisation term is then:
    #     L_{\mathrm{CDG}}(w) = \sum_{E \in \mathcal{E}_h^{\mathrm{int}}} \int_{E} - [\!\![ \theta ]\!\!]  \cdot \left< M \cdot (n \otimes n) \right > + \frac{1}{2} \frac{\alpha}{\left< h_E \right>} \left< \theta \cdot n \right> \cdot \left< \theta \cdot n \right> \; \mathrm{d}s
    
    # We choose the penalty parameter to be on the order of the norm of the bending stiffness matrix :math:`\dfrac{Et^3}{12}`.
    
    alpha = E*t**3
    h = CellDiameter(mesh)
    h_avg = (h('+') + h('-'))/2.0 
    
    n = FacetNormal(mesh)
    
    M_n = inner(M, outer(n, n))
    
    L_CDG = -inner(jump(theta, n), avg(M_n))*dS + \
                (1.0/2.0)*(alpha('+')/h_avg)*inner(jump(theta, n), jump(theta, n))*dS
    

    # In this example, we would like :math:`\theta_d = 0` everywhere on the boundary:: 
     
    # theta_d = Constant((0.0, 0.0))                                                  # To impose non-homogeneous BC
    # The definition of the exterior facets and Dirichlet rotation field were trivial
    # in this demo, but you could extend this code straightforwardly to
    # non-homogeneous Dirichlet conditions.
    # 
    # The weak boundary condition enforcement term can be written:
    #
    # .. math::
    #     L_{\mathrm{BC}}(w) = \sum_{E \in \mathcal{E}_h^{\mathrm{D}}} \int_{E} - \theta_e  \cdot (M \cdot (n \otimes n))  + \frac{1}{2} \frac{\alpha}{h_E} (\theta_e \cdot n)  \cdot (\theta_e \cdot n)  \; \mathrm{d}s
    # 
    # where :math:`\theta_e = \theta - \theta_d` is the effective rotation field, and
    # :math:`\mathcal{E}_h^{\mathrm{D}}` is the set of all exterior facets of the triangulation
    # :math:`\mathcal{T}` where we would like to apply Dirichlet boundary conditions, or in UFL::
    
    # theta_effective = theta - theta_d 
    # L_BC = -inner(inner(theta_effective, n), M_n)*ds(1) + \
    #         (1.0/2.0)*(alpha/h)*inner(inner(theta_effective, n), inner(theta_effective, n))*ds(1) 
    
    # The remainder of the demo is as usual::
    L = psi_M*dx + L_CDG 
    
    return L

def ComputeModesPlates(N, which_eig, Vh, L_, w, u, v , rho, t, bcs_w, *arg):
    import time 
    F = derivative(L_, w, v)
    J = derivative(F, w, u)
       
    K_ = PETScMatrix()
    assemble(J, tensor = K_)
    
    L  = rho*t*inner(u,v)*dx
    
    M_ = PETScMatrix()
    assemble(L, tensor = M_)
    
    # object is a list or not 
    if isinstance(bcs_w, list): 
        bcs_w  = bcs_w
    else:        
        bcs_w  = [bcs_w]
        
    for bcs in bcs_w:
        bcs.apply(K_)
        bcs.zero(M_) 

    solver = SLEPcEigenSolver(K_,M_)         #[𝐾]{𝑈} = 𝜆[𝑀]{𝑈}
    solver.parameters["problem_type"]       = "gen_hermitian"
    solver.parameters['spectral_transform'] = 'shift-and-invert'
    solver.parameters['spectral_shift']     = 1e-14
    solver.solve(N)
    k = solver.get_number_converged()
    
    # PETScOptions.set("eps_monitor_all")                   # monitor bugs
    from mf.plots import PlotSettingsSmall, get_axes_coord, PlotSettings, MidpointNormalize
    my_map = plt.get_cmap('RdBu')
    
    if which_eig == 'all':
        nn = int(N/25)
    else:
        nn = 1
    
    for ifig in range(nn):
        
        fig, myax = plt.subplots(5,5)
        fig.tight_layout(h_pad=0)
        # PlotSettings(fig, ax)
        cc =get_axes_coord(myax) 
        
        eig = Function(Vh)
        eig_vec = eig.vector()
        for jj in range(0,25):
            nn = jj + ifig*25
            r, c, rx, cx = solver.get_eigenpair(nn)
            eig_vec[:] = rx
            f = np.sqrt(np.real(r))/(2*np.pi)
            ax = myax[cc[jj]]
            
            plt.sca(ax)
            # PlotSettingsSmall(fig, ax)

            if np.abs(eig_vec.min())< 1e-10:
                vvin = -eig_vec.max()
            else:
                vvin = eig_vec.min()
            norm = MidpointNormalize(vmin=vvin, vmax=eig_vec.max(), midpoint=0)

            im = plot(eig, cmap=my_map, norm=norm, title=r'$\mathbf{Mode ~%.f}$'', %.2f Hz'%(nn+1, f))
            ax.set_aspect('equal', 'box')
            Lx = np.max(Vh.mesh().coordinates()[:,0]);          Ly = np.max(Vh.mesh().coordinates()[:,1])
            plt.xticks([0, Lx]);   plt.yticks([0, Ly])
            plt.draw()
            
            if len(arg) !=0:
                im.ax.set_aspect(arg[0])
            
            cb = plt.colorbar(im,  ax=ax)

        PlotSettings(fig, fig.axes)


def EigenSolvePlate(N, L_, w, u, v, rho, t, bcs_w):
    import time 
    F = derivative(L_, w, v)
    J = derivative(F, w, u)
       
    K_ = PETScMatrix()
    assemble(J, tensor = K_)
    
    L  = rho*t*inner(u,v)*dx
    
    M_ = PETScMatrix()
    assemble(L, tensor = M_)
    
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

def AverageSigma():
    path = '/home/carlos/dataPython/sigma_mn/' 
     
    # Initialize table to save computations
    import tables as tb
    M = 20
    N = 20

    Lx = 0.5
    Ly = 0.6
    t_ = 3e-3    

    # We take constant material properties throughout the domain::
    E_   = 7.1e10
    nu_  = 0.23
    D_   = E_*t_**3/(12*(1-nu_**2))
    eta  = 0.1
    rho  = 2700
    mass = rho*Lx*Ly*t_

    filename = path+'sigma_mn_'+str(M)+'_'+str(N)+'.h5'
    h5 = tb.open_file(filename, 'r')

    f  = (h5.root.frequency[:]).squeeze()
    Nf = len(f)
    c0 = 343.5  
    k0 = 2*np.pi*f/c0
    force_point = 1

    #  Analytical Solution : Modes in SS Plate
    Nm = M;                                 Nn = N
    m  =  np.r_[1:Nm+1];                    n  =np.r_[1:Nn+1]
    km = np.array([m*np.pi/Lx]);            kn = np.array([n*np.pi/Ly]);

    mass_a   = rho*Lx*Ly*t_
    wmn      = np. sqrt(D_/(rho*t_)) * (km.T**2+kn**2)

    imn      = 0
    omega_n  = np.zeros((M*N, 1))
    sigma_mn = np.zeros((Nf,M*N))

    for iim in range(M):
        m = iim+1
        for iin in range(N):
            n = iin+1
            omega_n[imn]    = wmn[iim, iin]
            integral        = h5.root.integral_mn[:,iim,iin]
            sigma_mn[:,imn] = 64*k0**2*Lx*Ly / (np.pi**6 * m**2 * n**2) * integral
            imn+=1

    print(sigma_mn.shape)  

    om = 2*np.pi*f.reshape(1,Nf)

    v4m = force_point**2*om**2/(2*mass**2) *1/((omega_n**2 - om**2)**2 + eta**2*omega_n**4)
    v2m = np.sum(v4m,axis=0)
    sigma = np.real( np.sum(sigma_mn*v4m.T,axis=1) / np.sum(v4m,axis=0))
    h5.close()
    return f, sigma


















