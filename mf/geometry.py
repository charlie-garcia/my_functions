import gmsh
import numpy as np
from mf.fem import gmsh2dolfin, write_gmsh

def CircularMesh(path, mesh_name, R, h, plot_info):
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
    write_gmsh(path,mesh_name)                                       
    # Don't forget to finalize gmsh
    gmsh.finalize()  

    fmesh, mf_boundary = gmsh2dolfin(path, mesh_name, '2D', bord_string_tag, surface_string_tag)            

    return fmesh, mf_boundary, tag_bords

def RectangularMesh(path, mesh_name, Lx,Ly, h, plot_info):
    # Create membrane mesh with Gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add(mesh_name)
    factory = gmsh.model.geo
    
    # Create geometry
    # Add Points, px, py, pz, resolution, tag 
    factory.addPoint(0,   0,  0, h, 1)                          
    factory.addPoint(Lx,  0,  0, h, 2)
    factory.addPoint(Lx,  Ly, 0, h, 3)
    factory.addPoint(0,   Ly, 0, h, 4)
    
    # Add arcs, p1, center, p2, tag
    factory.addLine(1, 2, 10)
    factory.addLine(2, 3, 11)
    factory.addLine(3, 4, 12)
    factory.addLine(4, 1, 13)

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
    write_gmsh(path,mesh_name)                                                                                                                  
    # Don't forget to finalize gmsh
    gmsh.finalize()  

    fmesh, mf_boundary = gmsh2dolfin(path, mesh_name, '2D', bord_string_tag, surface_string_tag)            

    return fmesh, mf_boundary, tag_bords