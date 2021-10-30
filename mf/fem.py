import numpy as np
import scipy.interpolate as inter
from dolfin import *
import meshio
import matplotlib.pyplot as plt

def fem2cart(V, u, mesh, npoints, offset):
    
    nodal_values_u = u.vector()
    array_u = np.array(nodal_values_u)
    
    xx, yy = node2coord(V,mesh)

    [X, Y] = np.meshgrid (np.linspace(np.min(xx)-offset, np.max(xx) + offset, npoints), 		
                         np.linspace(np.min(yy)-offset, np.max(yy)+offset,npoints) )
    Z = inter.griddata(np.c_[xx,yy], array_u,(X,Y), method = 'cubic')
    
    return X,Y,Z, xx, yy

def node2coord(V, mesh):
    # Acces to coordinates
    n = V.dim()                                                                      
    d = mesh.geometry().dim()                                                        
    dof_coordinates = V.tabulate_dof_coordinates().reshape(n,d)
    
    xx = dof_coordinates[:,0]
    yy = dof_coordinates[:,1]

    return xx, yy

def fem3coord(V, u, mesh):
    nodal_values_u = u.vector()
    array_u = np.array(nodal_values_u)
    
    # Acces to coordinates
    n = V.dim()                                                                      
    d = mesh.geometry().dim()                                                        
    dof_coordinates = V.tabulate_dof_coordinates().reshape(n,d)
    
    xx = dof_coordinates[:,0]
    yy = dof_coordinates[:,1]    
    zz = dof_coordinates[:,2]
    
    return xx, yy, zz

def eig2cart(eig_vec, V, mesh, npoints, offset):
    array_u = eig_vec[:]
    
    # Acces to coordinates
    xx, yy = node2coord(V,mesh)
    
    [X, Y] = np.meshgrid(np.linspace(np.min(xx)-offset, np.max(xx) + offset, npoints), 	np.linspace(np.min(yy)-offset, np.max(yy)+offset,npoints));
    mode = inter.griddata([xx, yy],array_u,(X,Y), method = 'cubic')
    
    return mode,xx,yy

def fem3cart(V, u, mesh, npoints, offset):
    
    nodal_values_u = u.vector()
    array_u = np.array(nodal_values_u)
    
    # Acces to coordinates
    n = V.dim()                                                                      
    d = mesh.geometry().dim()                                                        
    dof_coordinates = V.tabulate_dof_coordinates().reshape(n,d)
    
    xx = dof_coordinates[:,0]
    yy = dof_coordinates[:,1]
    zz = dof_coordinates[:,2]
    

    [X, Y] = np.meshgrid(np.linspace(np.min(xx) - offset, np.max(xx) + offset, npoints),   \
                            np.linspace(np.min(yy) - offset, np.max(yy) + offset, npoints) );
    C = inter.griddata(dof_coordinates, array_u,(X,Y), method = 'linear')
    
    return X,Y,Z, C, xx, yy

def gmsh2dolfin(path, mesh_name, dim, bord_string_tag, surface_string_tag):
    my_mesh = meshio.read(path+mesh_name)
    if dim=='2D':
        set_prune_z = True
    elif dim=='3D':
        set_prune_z = False

    def create_mesh(mesh, cell_type, my_tag, prune_z=False):
        cells = np.vstack([cell.data for cell in mesh.cells if cell.type==cell_type])
        cell_data = np.hstack([mesh.cell_data_dict["gmsh:physical"][key]
                              for key in mesh.cell_data_dict["gmsh:physical"].keys() if key==cell_type])
    
        # Remove z-coordinates from mesh if we have a 2D cell and all points have the same third coordinate
        points= mesh.points
        if prune_z:
            points = points[:,:2]
        mesh_new = meshio.Mesh(points=points, cells={cell_type: cells}, cell_data={my_tag:[cell_data]})
        return mesh_new
    
    border_mesh_name_xdmf = "border_fmv.xdmf"
    dom_mesh_name_xdmf = "domains_fmesh.xdmf"
    
    line_border_mesh = create_mesh(my_mesh, "line", bord_string_tag, prune_z=set_prune_z)
    meshio.write(path + border_mesh_name_xdmf, line_border_mesh)
    
    triangle_mesh = create_mesh(my_mesh, "triangle", surface_string_tag, prune_z=set_prune_z)  # Change to false if wanna 2D
    meshio.write(path + dom_mesh_name_xdmf, triangle_mesh)
    
    # 
    fmesh = Mesh()
    # plot(fmesh)
    with XDMFFile(path + dom_mesh_name_xdmf) as infile:
        infile.read(fmesh)
    
    mvc = MeshValueCollection("size_t", fmesh, 1)
    
    with XDMFFile(path + border_mesh_name_xdmf) as infile:
        print("Reading 1d line data into dolfin mvc")
        infile.read(mvc, bord_string_tag)
    
    print("Constructing MeshFunction from MeshValueCollection")
    mf = MeshFunction('size_t',fmesh, mvc)
    
    return fmesh, mf

def gmsh3dolfin(path, mesh_name):
    
    msh = meshio.read(path+mesh_name)

    def create_mesh(mesh, cell_type):
        cells = np.vstack([cell.data for cell in mesh.cells if cell.type==cell_type])
        cell_data = np.hstack([mesh.cell_data_dict["gmsh:geometrical"][key]
                              for key in mesh.cell_data_dict["gmsh:geometrical"].keys() if key==cell_type])
        
        mesh_new = meshio.Mesh(points=msh.points, cells={cell_type: cells}, cell_data={"name_to_read":[cell_data]})
        return mesh_new
    
    surface_mesh_name_xdmf = "surface_fmesh.xdmf"
    dom_mesh_name_xdmf     = "domains_fmesh.xdmf"
    
    triangle_mesh = create_mesh(msh, "triangle")            
    meshio.write(path + surface_mesh_name_xdmf, triangle_mesh)
    
    tetra_mesh = create_mesh(msh, "tetra")            
    meshio.write(path + dom_mesh_name_xdmf, triangle_mesh)
    
    meshio.write(path+"mesh.xdmf", meshio.Mesh(points=msh.points, cells={"tetra": 
                    np.vstack([cell.data for cell in msh.cells if cell.type=="tetra"])}))
        
    fmesh = Mesh()
    
    with XDMFFile(path+"mesh.xdmf") as infile:
        infile.read(fmesh)
    
    mvc = MeshValueCollection("size_t", fmesh, 2)
    
    with XDMFFile(path + surface_mesh_name_xdmf) as infile:
        print("Reading 1d line data into dolfin mvc")
        infile.read(mvc, "name_to_read")
    
    print("Constructing MeshFunction from MeshValueCollection")
    
    mf = MeshFunction('size_t', fmesh, mvc)
     
    mvc2 = MeshValueCollection("size_t", fmesh, 3)
    
    with XDMFFile(path + dom_mesh_name_xdmf) as infile:
        infile.read(mvc, "name_to_read")
        
    cf = MeshFunction('size_t', fmesh, mvc2)

    return fmesh, mf, cf

def gmsh2dolfin_subd(path, mesh_name, dim, bord_string_tag, surface_string_tag):
    my_mesh = meshio.read(path+mesh_name)
    if dim=='2D':
        set_prune_z = True
    elif dim=='3D':
        set_prune_z = False
    
    def create_mesh(mesh, cell_type, my_tag, prune_z=False):
        cells = np.vstack([cell.data for cell in mesh.cells if cell.type==cell_type])
        cell_data = np.hstack([mesh.cell_data_dict["gmsh:physical"][key]
                              for key in mesh.cell_data_dict["gmsh:physical"].keys() if key==cell_type])
    
        # Remove z-coordinates from mesh if we have a 2D cell and all points have the same third coordinate
        points= mesh.points
        if prune_z:
            points = points[:,:2]
        mesh_new = meshio.Mesh(points=points, cells={cell_type: cells}, cell_data={my_tag:[cell_data]})
        return mesh_new
    
    border_mesh_name_xdmf = "border_fmv.xdmf"
    dom_mesh_name_xdmf = "domains_fmesh.xdmf"
    
    line_border_mesh = create_mesh(my_mesh, "line", bord_string_tag, prune_z=set_prune_z)
    meshio.write(path + border_mesh_name_xdmf, line_border_mesh)
    
    domains_mesh = create_mesh(my_mesh, "triangle", surface_string_tag, prune_z=set_prune_z)  # Change to false if wanna 2D
    meshio.write(path + dom_mesh_name_xdmf, domains_mesh)
    
    #
    fmesh = Mesh()
    
    with XDMFFile(path + dom_mesh_name_xdmf) as infile:
        infile.read(fmesh)
    
    mvc_border = MeshValueCollection("size_t", fmesh, 1)
    
    with XDMFFile(path + border_mesh_name_xdmf) as infile:
        print("Reading 1d line data into dolfin mvc")
        infile.read(mvc_border, bord_string_tag)
    
    # Load Subdomains
    mvc_domains= MeshValueCollection("size_t", fmesh, 2)
    
    with XDMFFile(path + dom_mesh_name_xdmf) as infile:
        print("Reading 2d line data into subdomains")
        infile.read(mvc_domains, surface_string_tag)
        
    print("Constructing MeshFunction from MeshValueCollection")
    mf_boundary = MeshFunction('size_t',fmesh, mvc_border)
    mf_domains = MeshFunction('size_t',fmesh, mvc_domains)
    
    return fmesh, mf_boundary, mf_domains

def gmsh2grid(path, mesh_name, bord_string_tag, surface_string_tag):
    my_mesh = meshio.read(path+mesh_name)

    def create_mesh(mesh, cell_type, my_tag, prune_z=True):  # Prune for 2D FEM
        cells = np.vstack([cell.data for cell in mesh.cells if cell.type==cell_type])
        cell_data = np.hstack([mesh.cell_data_dict["gmsh:physical"][key]
                              for key in mesh.cell_data_dict["gmsh:physical"].keys() if key==cell_type])
    
        # Remove z-coordinates from mesh if we have a 2D cell and all points have the same third coordinate
        points= mesh.points
        if prune_z:
            points = points[:,:2]
        mesh_new = meshio.Mesh(points=points, cells={cell_type: cells}, cell_data={my_tag:[cell_data]})
        return mesh_new
    
    dom_mesh_name_xdmf = "domains_fmesh.xdmf"
    triangle_mesh = create_mesh(my_mesh, "triangle", surface_string_tag)  # Change to false if wanna 2D
    meshio.write(path + dom_mesh_name_xdmf, triangle_mesh)

    grid = Mesh()
    with XDMFFile(path + dom_mesh_name_xdmf) as infile:
        infile.read(grid)

    return grid
    
def connect_triangles_fem(V, u, mesh, element, plot_info):
    if element == 'dof':
        n = V.dim()                                                     # n nodes
        d = mesh.geometry().dim()                                                        
        dof_coordinates = V.tabulate_dof_coordinates().reshape(n,d)
        xx = dof_coordinates[:,0]
        yy = dof_coordinates[:,1]
        coordinates = dof_coordinates
        
    elif element =='mesh':
        xx = mesh.coordinates()[:,0]
        yy = mesh.coordinates()[:,1]
        coordinates = mesh.coordinates()
        
    spl = mesh.cells()
    cs, se = getCentersTriangles(xx,yy,spl)

    if plot_info=='plot':
        import matplotlib.pyplot as plt
        plt.triplot(xx, yy, spl)
        plt.plot(xx, yy, 'o')
        plt.plot(cs[:,0], cs[:,1], 'rx')

    return cs, se

def connect_triangles_grid(mesh, plot_info):
    xx = mesh.points[:,0] 
    yy = mesh.points[:,1] 
    coordinates = mesh.points[:,0:2] 
    spl = mesh.cells_dict['triangle']
    cs, se = getCentersTriangles(xx,yy,spl)

    if plot_info=='plot':
        plt.triplot(xx, yy, tri.simplices)
        plt.plot(xx, yy, 'o')
        plt.plot(cs[:,0], cs[:,1], 'rx')

    return cs, se

def getCentersTriangles(xx, yy, spl):
    from scipy.linalg import det
    # get baricenters
    x, y, z= (np.zeros(( len(spl), 3)) for ii in range(0,3))
    for jj in range(0, len(spl)):
        x[jj] = xx[spl[jj,:]]
        y[jj] = yy[spl[jj,:]]

    cx = np.mean(x,1)
    cy = np.mean(y,1)
    cz = 0*cy
    cs = np.c_[cx, cy, cz]

    # get triangle areas
    se = np.zeros( (x.shape[0], ) )
    for jj in range(0, len(se)):
        uv = np.array([ [x[jj][1] - x[jj][0], y[jj][1] - y[jj][0], z[jj][1]-z[jj][0]], \
                        [x[jj][2] - x[jj][0], y[jj][2] - y[jj][0], z[jj][2]-z[jj][0]] ])

        uv0 = np.array([uv[:,1], uv[:,2]])
        uv1 = np.array([uv[:,0], uv[:,2]])
        uv2 = np.array([uv[:,0], uv[:,1]])

        se[jj] = np.sqrt( det(uv0)**2 + det(uv1)**2 + det(uv2)**2 ) /2
        
    return cs, se.reshape(len(cs),1)

# def getCentersTriangles(xx, yy, spl):
#     from scipy.linalg import det
#     # get baricenters
#     x, y = (np.zeros(( len(spl), 3)) for ii in range(2))
#     for jj in range(0, len(spl)):
#         x[jj] = xx[spl[jj,:]]
#         y[jj] = yy[spl[jj,:]]

#     cx = np.mean(x,1)
#     cy = np.mean(y,1)
#     cz = 0*cy
#     cs = np.c_[cx, cy, cz]

#     # get triangle areas
#     se = np.zeros( (x.shape[0], ) )
#     # A = (x1y2 + x2y3 + x3y1 – x1y3 – x2y1 – x3y2)/2.
#     for jj in range(0, len(se)):
#         se[jj] = np.sqrt( x[jj][0]* y[jj][1] + x[jj][1]* y[jj][2] + x[jj][2]* y[jj][0] -\
#                           x[jj][0]* y[jj][2] - x[jj][1]* y[jj][0] + x[jj][2]* y[jj][1] ) /2 
#     return cs, se

def Clamped_Plate(W, w_, mesh, E_, nu_, t_, force, ds):
    tol = 1E-14

    class Omega_0(SubDomain):
        def inside(self, x, on_boundary):
            return x[1] <= 0.5 + tol
    
    class Omega_1(SubDomain):
        def inside(self, x, on_boundary):
            return x[1] >= 0.5 - tol
        
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
    L_BC = -inner(inner(theta_effective, n), M_n)*ds(1) + \
            (1.0/2.0)*(alpha/h)*inner(inner(theta_effective, n), inner(theta_effective, n))*ds(1) 
    
    # The remainder of the demo is as usual::
    
    f = Constant(force)
    W_ext = f*w_*dx
    
    L = psi_M*dx - W_ext + L_CDG + L_BC   
    
    return L

def SS_Plate(W, w_, mesh, E_, nu_, t_, force):
   
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
    
    f = Constant(force)
    W_ext = f*w_*dx
    
    L = psi_M*dx - W_ext + L_CDG 
    
    return L

def ComputeModesPlates(number_of_requested_eigenpairs, which_eig, J, W, w, w_t, rho, t, bcs_w, *arg):
    
    a = J
    A2 = PETScMatrix()
    assemble(a, tensor = A2)
    
    # # element dkt kirchhoff - love
    L = rho*t*inner(w,w_t)*dx
    B = PETScMatrix()
    assemble(L, tensor = B)
    
    # code to check whether 
    # object is a list or not 
    if isinstance(bcs_w, list): 
        bcs_w  = bcs_w
    else:        
        bcs_w  = [bcs_w]
        
    for bc in bcs_w:
        bc.apply(A2)
        bc.zero(B)                          # avoid spurius eig vals
    # u = Function(V)
    
    # px = 0.1
    # py = 0
    #PointSource(V, Point(px,py), f).apply(B)           # See how to use the delta function instead
    
    # solve(A, u.vector(), b)
    solver = SLEPcEigenSolver(A2,B)                     #[𝐾]{𝑈} = 𝜆[𝑀]{𝑈}
    solver.parameters["problem_type"] = "gen_hermitian"
    solver.parameters['spectral_transform'] = 'shift-and-invert'
    solver.parameters['spectral_shift'] = 1e-14
    solver.solve(number_of_requested_eigenpairs)
    k = solver.get_number_converged()
    
    # PETScOptions.set("eps_monitor_all")                   # monitor bugs
    from mf.plots import PlotSettingsSmall, get_axes_coord, ColorbarSettings
    my_map = plt.get_cmap('seismic')
    #%%
    if which_eig == 'all':
        nn = int(number_of_requested_eigenpairs/25)
    else:
        nn = 1
    
    for ifig in range(nn):
        
        fig,myax = plt.subplots(5,5)
        fig.tight_layout(h_pad=0)
        # PlotSettings(fig, ax)
        cc =get_axes_coord(myax) 
        
        eig = Function(W)
        eig_vec = eig.vector()
        for jj in range(0,25):
            nn = jj + ifig*25
            r, c, rx, cx = solver.get_eigenpair(nn)
            eig_vec[:] = rx
            f = np.sqrt(np.real(r))/(2*np.pi)
            ax = myax[cc[jj]]
            
            plt.sca(ax)
            PlotSettingsSmall(fig, ax)
            im = plot(eig, cmap=my_map, title=r'$\mathbf{Mode ~%.f}$'', %.2f Hz'%(nn, f))
            ax.set_aspect('equal', 'box')
            plt.draw()
            
            if len(arg) !=0:
                im.ax.set_aspect(arg[0])
            
            cb = plt.colorbar(im,  ax=ax)
            ticklabs = cb.ax.get_yticklabels()
            cb.ax.set_yticklabels(ticklabs, fontsize=7)


def CreateMembrane(mesh_name, R, h, dim, plot_info):
    import gmsh
    import numpy as np
    
    # Create membrane mesh with Gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add(mesh_name)
    factory = gmsh.model.geo
    
    # Create geometry
    # Add Points, px, py, pz, resolution, tag 
    factory.addPoint(-R,  0 , 0, h, 1)                          
    factory.addPoint(0,   0 , 0, h, 2)
    factory.addPoint(R,   0 , 0, h, 3)
    
    factory.addPoint(0,  -R, 0, h, 4)
    factory.addPoint(0,   R, 0, h, 5)
    
    # Add arcs, p1, center, p2, tag
    factory.addCircleArc(3, 2, 4, 10)                            
    factory.addCircleArc(4, 2, 1, 11)
    
    factory.addCircleArc(1, 2, 5, 12)
    factory.addCircleArc(5, 2, 3, 13)
    
    # Add curve loop (closes the boundary) [tag_line1,.. tag_lineN], tag 
    id_border = factory.addCurveLoop([10, 11, 12, 13])                  
    
    # create a planar surface with tag = id_surface
    id_surface = factory.addPlaneSurface([id_border])                  
    
    # sync changes
    factory.synchronize()                                       
    # create mesh
    gmsh.model.mesh.generate(2)                                 
    
    # add a group name to the surface
    surface_string_tag = 'my_surface'
    tag_surface = gmsh.model.addPhysicalGroup(2, [id_surface])    
    # set the name to element in this sirface, (dim, tag, name)
    gmsh.model.setPhysicalName(2, tag_surface, surface_string_tag)        
    
    # add a group name to the boundary
    bord_string_tag = 'my_bords'
    tag_bords  = gmsh.model.addPhysicalGroup(1, [10, 11, 12, 13])                    # REMEMBER THIS TAG FOR BC
    # set the name to element in this sirface, (dim, tag, name)
    gmsh.model.setPhysicalName(1, tag_bords, bord_string_tag)                  # dim, tag, name

    
    # Comment/Uncomment to see the mesh in gmsh
    if plot_info == 'plot':
        gmsh.fltk.run()                                             
    # Write mesh
    gmsh.write(mesh_name)                                       
    # Don't forget to finalize gmsh
    gmsh.finalize()           

    fmesh, mf = gmsh2dolfin('', mesh_name, dim, bord_string_tag, surface_string_tag)            

    return fmesh, mf, tag_bords


def write_gmsh(path, mesh_name):
    import os.path, gmsh
    
    if os.path.exists(path+mesh_name):
        while True:
            try:
                overwrite = str(input("Existing file, wanna replace-it? [y/n/c]: "))
            except ValueError:
                print("Sorry, I didn't understand that.")                       # Return to the start of the loop
                continue
            else:
                break                                                           #we're ready to exit the loop.

        if 'y' in overwrite: 
            gmsh.write(path+mesh_name)
            print('Overwriting existing file')
        elif 'c' in overwrite:
            print('Skipping mesh creation, using old mesh: '+mesh_name)

        else:
            raise ValueError("Please change the name of the file")
    
    else:
        print ("Creating new file")
        gmsh.write(path+mesh_name)
