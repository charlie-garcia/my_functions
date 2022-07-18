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

def fem2cart_ds(V, u, mesh, ds_cart):

    nodal_values_u = u.vector()
    array_u = np.array(nodal_values_u)
    
    xx, yy = node2coord(V,mesh)
    npts_x = int((np.max(xx) - np.min(xx))/ds_cart )
    npts_y = int((np.max(yy) - np.min(yy))/ds_cart)
    [X, Y] = np.meshgrid (np.linspace(np.min(xx), np.max(xx), npts_x),         
                         np.linspace(np.min(yy), np.max(yy), npts_y) )
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
    
    border_mesh_name_xdmf = mesh_name[:-4]+"_border_fmv.xdmf"
    dom_mesh_name_xdmf = mesh_name[:-4]+"_domains_fmesh.xdmf"
    
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

def ReadFenicsMesh(path, mesh_name, bord_string_tag):
    # ! This does not return tag's boundary 
    border_mesh_name_xdmf = mesh_name[:-4]+"_border_fmv.xdmf"
    dom_mesh_name_xdmf = mesh_name[:-4]+"_domains_fmesh.xdmf"
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

def ReadFenicsMeshMoved(path, mesh_name):
    # to use this, save the moved mesh with 
    # mesh_name = 'some_name.msh9'
    # File(path+mesh_name[:-4]+'_fmesh_moved.xml') << fmesh  # the moved mesh
    # File(path+mesh_name[:-4]+'_mf_moved.xml') << mf        # the moved meshfunction

    fmesh = Mesh(path+mesh_name[:-4]+'_fmesh_moved.xml')
    mf_boundary = MeshFunction('size_t', fmesh, path+mesh_name[:-4]+'_mf_moved.xml')
    return fmesh, mf_boundary

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

def PlotMode(solver, Vh, ax, jj):
    eig = Function(Vh)
    from mf.plots import MidpointNormalize
    r, c, rx, cx = solver.get_eigenpair(jj)
    eig.vector()[:] = rx
    f = np.sqrt(np.real(r))/(2*np.pi)
    
    eig_min, eig_max = eig.vector()[:].min(), eig.vector()[:].max()
    vvmin, vvmax =  eig_min, eig_max
    
    if np.isclose(eig_min, 0.) or np.isclose(eig_max, 0.) :
        vvmin, vvmax = -np.abs((eig_min+eig_max)), np.abs((eig_min+eig_max))
        
    norm = MidpointNormalize(vmin=vvmin, vmax=vvmax, midpoint=0)
    
    plt.sca(ax)
    my_map = plt.get_cmap('RdBu_r')
    im = plot(eig, cmap=my_map, norm=norm, title=r'$\mathbf{Mode ~%.f}$'', %.2f Hz'%(jj+1, f))
    ax.set_aspect('equal', 'box')
    
    Lx = np.round(np.max(Vh.mesh().coordinates()[:,0]),2);  Lx0 = np.round(np.min(Vh.mesh().coordinates()[:,0]),2);          
    Ly = np.round(np.max(Vh.mesh().coordinates()[:,1]),2);  Ly0 = np.round(np.min(Vh.mesh().coordinates()[:,1]),2);          
    plt.xticks([Lx0, Lx]);   plt.yticks([Ly0, Ly])

    im_ratio = np.min([Lx,Ly])/np.max([Lx,Ly])
    cb = plt.colorbar(im,  ax=ax, ticks = [0],fraction=0.046*im_ratio, pad=0.05)
    
    plt.show()

    return f, eig