import gmsh, sys, meshio 
import numpy as np
from mf.fem import gmsh2dolfin, gmsh2dolfin_subd

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


def CreateComplexPate(my_path, mesh_name, Lx, Ly, coeff, loc_pts, exc_points):
    import gmsh, sys
    # Mesh generation with GMSH
    gmsh.initialize(sys.argv)

    # Ask GMSH to display information in the terminal
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add(mesh_name)

    model = gmsh.model
    model.add("MyPlate")

    # Create Rectangle
    bx = [0,Lx, Lx, 0]
    by = [0, 0, Ly, Ly]
    h1  = 0.08
    dim = '2D'
    
    vertex, borders, ptA, ptB = ( [] for i in range(4))

    for j in range(len(bx)):
        vertex.append(model.geo.addPoint( bx[j], by[j], 0, h1))

    for j in range(len(bx)):
        if j < len(bx):
            borders.append(gmsh.model.geo.addLine(vertex[j-1],vertex[j]))

    # Curveloop and Surface
    curveloop = model.geo.addCurveLoop(borders)
    id_surface = model.geo.addPlaneSurface([curveloop])
    
        
    # This command is mandatory and synchronize CAD with GMSH Model. The less you launch it, the better it is for performance purpose
    gmsh.model.geo.synchronize()

    # Add fixed points 
    loc_x = loc_pts[:,0]
    loc_y = loc_pts[:,1]

    l = []
    tol = 1e-4
    for j in range(len(loc_x)):
        ptA.append(model.geo.addPoint( loc_x[j],     loc_y[j],     0, h1))
        ptB.append(model.geo.addPoint( loc_x[j]+tol, loc_y[j]+tol, 0, h1))
        l.append(model.geo.addLine(ptA[j], ptB[j]))

    gmsh.model.geo.synchronize()

    for j in range(len(l)):
        gmsh.model.mesh.embed(1, [l[j]], 2, id_surface)   

    # Embedded ponts into the surface
    px = exc_points[:,0]
    py = exc_points[:,1]
    points = []
    
    for j in range(len(px)):
        points.append(model.geo.addPoint(px[j], py[j], 0, h1))

    gmsh.model.geo.synchronize()
    gmsh.model.mesh.embed(0, points, 2, id_surface)             # dim =1, line  
    
    # set algorithm "Packing of parallelograms" (experimental =9)
    gmsh.model.mesh.setAlgorithm(2, id_surface, 9)
    gmsh.option.setNumber('Mesh.MeshSizeFactor', coeff)
    gmsh.model.mesh.generate(2)                                 # 2D mesh

    #% Add physical groups
    # Bord phisical group
    line_string_tag = "border"
    tag_border = gmsh.model.addPhysicalGroup(1, borders)         
    gmsh.model.setPhysicalName(1, tag_border, line_string_tag)   # dim, tag, name

    # # Fixed point phisical group
    tag_points = []
    for ii in range(len(l)):
        tag_points.append(gmsh.model.addPhysicalGroup(1, [l[ii]]))              # REMEMBER THIS TAG FOR BC
        gmsh.model.setPhysicalName(1, tag_points[ii], line_string_tag)   # dim, tag, name

    # Surface phisical group
    surface_string_tag = "surface"
    tag_dom = gmsh.model.addPhysicalGroup(2, [id_surface])       # Delete original tag when fragment
    gmsh.model.setPhysicalName(2, tag_dom, surface_string_tag)   # dim, tag, name

    gmsh.write(my_path+mesh_name)
    # gmsh.fltk.run()
    
    gmsh.finalize()

    fmesh, boundaries = gmsh2dolfin(my_path, mesh_name, dim, line_string_tag, surface_string_tag)

    return fmesh, boundaries, tag_border, tag_points

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


















