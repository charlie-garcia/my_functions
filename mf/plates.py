import gmsh
import sys
import numpy as np
import meshio
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












