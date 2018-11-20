from math import *
import sys
import numpy as np
import random
import config
import collections
import copy
import svgwrite
import matplotlib.tri as tri
import svgutils.transform as st
import matplotlib.pyplot as plt
from openmesh_example import *

## Manage paths
sys.path.append("/Users/jtbondorf/bin/OpenMesh/build/Build/python")             # Add OpenMesh library to Python's search path

## Import OpenMesh
from openmesh import *



#---------------------------------------------------------------------------
# Functions to do a parameter test on number of iterations and step size for
# the equilibriate_vertices method. Moved to separate file to save space
#---------------------------------------------------------------------------


def prepare_mesh(mesh):
    """ Read mesh and add properties
        TriMesh mesh: mesh to be prepared
    """
    read_mesh(mesh, "initial_foam_mesh.obj")
    add_mid_point(mesh)
    add_pressure(mesh)
    remove_sliver_faces(mesh)
    remove_shark_fin_faces(mesh)


def compute_vertex_distances(mesh):
    """ Compute array of distances between all vertices (some distances appear twice)
        TriMesh mesh: dual mesh
        return array distances: distances between vertices
    """

    distances = []

    for vh0 in mesh.vertices():

        for vh in mesh.vv(vh0):
            dist = np.sqrt((mesh.point(vh)[1]-mesh.point(vh0)[1])**2 + (mesh.point(vh)[1]-mesh.point(vh0)[1])**2)
            distances.append(dist)

    return distances


def comprehensive_parameter_test(mesh):
    """ Do parameter test on average_equilibrate, spring_equilibrate
        TriMesh mesh: an existing dual OpenMesh mesh
    """

    # Constants
    par_range = 2

    # Save initial mesh
    write_mesh(mesh, "initial_foam_mesh.obj")


    #-------------------------------------------------------------------------------------------------------------------
    # Original mesh
    #-------------------------------------------------------------------------------------------------------------------
    temp_mesh = TriMesh()
    read_mesh(temp_mesh, "initial_foam_mesh.obj")
    add_mid_point(temp_mesh)
    add_pressure(temp_mesh)
    draw_primary_mesh(temp_mesh, title='parameter_test/original', draw_boundary=False)

    temp_mesh = TriMesh()
    prepare_mesh(temp_mesh)
    draw_primary_mesh(temp_mesh, title='parameter_test/original_nosliver', draw_boundary=False)



    #-------------------------------------------------------------------------------------------------------------------
    # Average equilibriate
    #-------------------------------------------------------------------------------------------------------------------

    iter_list  = [10,100]
    scale_list = [0.05,0.2]

    for i in range(par_range):
        for j in range(par_range):
            temp_mesh = TriMesh()
            prepare_mesh(temp_mesh)
            average_equilibriate_vertices(temp_mesh,scale_list[j],iter_list[i],False)
            #draw_primary_mesh(temp_mesh, title='parameter_test/primary_mesh_after_average_equilibration_scale'+str(scale_list[j])+'_iter'+str(iter_list[i]),draw_boundary=False)
            draw_primary_mesh(temp_mesh, title='parameter_test/ae_' + str((i * par_range) + j), draw_boundary=False)



    # -------------------------------------------------------------------------------------------------------------------
    # Average equilibriate, fix boundaries
    # -------------------------------------------------------------------------------------------------------------------

    # Solve sliver face problem
    for i in range(par_range):
        for j in range(par_range):
            temp_mesh = TriMesh()
            prepare_mesh(temp_mesh)
            average_equilibriate_vertices(temp_mesh, scale_list[j], iter_list[i])
            #draw_primary_mesh(temp_mesh, title='parameter_test/primary_mesh_after_average_equilibration_fixed_scale' + str(scale_list[j]) + '_iter' + str(iter_list[i]), draw_boundary=False)
            draw_primary_mesh(temp_mesh,title='parameter_test/aef_'+str((i * par_range)+j), draw_boundary=False)



    #-------------------------------------------------------------------------------------------------------------------
    # Spring equilibriate
    #-------------------------------------------------------------------------------------------------------------------

    distances = compute_vertex_distances(mesh)

    l_0_list  = [np.mean(distances)-0.1*np.std(distances), np.mean(distances)+0.1*np.std(distances)]
    k_list    = [1, 10]

    print l_0_list

    for i in range(par_range):
        for j in range(par_range):
            temp_mesh = TriMesh()
            prepare_mesh(temp_mesh)
            spring_equilibriate_vertices(temp_mesh, k=k_list[j], l_0=l_0_list[i], delta_t=0.005, max_iterations=10)
            #draw_primary_mesh(temp_mesh, title='parameter_test/primary_mesh_after_spring_equilibration_fixed_l0' + str(l_0_list[j]) + '_deltat' + str(delta_t_list[i]), draw_boundary=False)
            draw_primary_mesh(temp_mesh, title='parameter_test/se_' + str((i * par_range) + j), draw_boundary=False)



    #-------------------------------------------------------------------------------------------------------------------
    # Combined equilibration
    #-------------------------------------------------------------------------------------------------------------------
    temp_mesh = TriMesh()
    prepare_mesh(temp_mesh)
    combine_equilibriate_vertices(temp_mesh)
    draw_primary_mesh(temp_mesh, title='parameter_test/comb',draw_boundary=False)


    temp_mesh = TriMesh()
    prepare_mesh(temp_mesh)
    combine_equilibriate_vertices(temp_mesh)
    add_pressure(temp_mesh, equal=False)
    draw_primary_mesh(temp_mesh, title='parameter_test/comb_nonequal')




def equilibrium_parameter_test():
    
    min = config.length['min']
    max = config.length['max']
    
    xv, yv       = init_grid_rand(min,max,5,False)                              # Initialize random grid
    triang       = tri.Triangulation(xv, yv)                                    # Do Delauney-triangulation on the set of x- and y-coordinates
    faces        = triang.get_masked_triangles()                                # Extracts faces from the Delaunay-triangulation from the matplotlib.tri-package
    numver       = len(xv)
    numfac       = len(faces)
    vertices     = np.array([[xv[i], yv[i], 0] for i in range(numver)])
    
    # Original mesh
    foam_mesh    = TriMesh()
    vhs          = [foam_mesh.add_vertex(TriMesh.Point(float(vertices[i][0]),float(vertices[i][1]),float(vertices[i][2]))) for i in range(numver)]
    for i in range(numfac):
        foam_mesh.add_face(vhs[faces[i][0]], vhs[faces[i][1]], vhs[faces[i][2]])
    # Add mid points to each face
    add_mid_point(foam_mesh)
    add_pressure(foam_mesh)


    # Original mesh 0
    original_mesh_0    = TriMesh()
    vhs                = [original_mesh_0.add_vertex(TriMesh.Point(float(vertices[i][0]),float(vertices[i][1]),float(vertices[i][2]))) for i in range(numver)]
    for i in range(numfac):
        original_mesh_0.add_face(vhs[faces[i][0]], vhs[faces[i][1]], vhs[faces[i][2]])
    # Add mid points to each face
    add_mid_point(original_mesh_0)
    add_pressure(original_mesh_0)
    
    
    # Original mesh 1
    original_mesh_1    = TriMesh()
    vhs                = [original_mesh_1.add_vertex(TriMesh.Point(float(vertices[i][0]),float(vertices[i][1]),float(vertices[i][2]))) for i in range(numver)]
    for i in range(numfac):
        original_mesh_1.add_face(vhs[faces[i][0]], vhs[faces[i][1]], vhs[faces[i][2]])
    # Add mid points to each face
    add_mid_point(original_mesh_1)
    add_pressure(original_mesh_1)
    
    
    # Original mesh 2
    original_mesh_2    = TriMesh()
    vhs                = [original_mesh_2.add_vertex(TriMesh.Point(float(vertices[i][0]),float(vertices[i][1]),float(vertices[i][2]))) for i in range(numver)]
    for i in range(numfac):
        original_mesh_2.add_face(vhs[faces[i][0]], vhs[faces[i][1]], vhs[faces[i][2]])
    # Add mid points to each face
    add_mid_point(original_mesh_2)
    add_pressure(original_mesh_2)
    
    
    # Original mesh 3
    original_mesh_3    = TriMesh()
    vhs                = [original_mesh_3.add_vertex(TriMesh.Point(float(vertices[i][0]),float(vertices[i][1]),float(vertices[i][2]))) for i in range(numver)]
    for i in range(numfac):
        original_mesh_3.add_face(vhs[faces[i][0]], vhs[faces[i][1]], vhs[faces[i][2]])
    # Add mid points to each face
    add_mid_point(original_mesh_3)
    add_pressure(original_mesh_3)


    # Original mesh 4
    original_mesh_4    = TriMesh()
    vhs                = [original_mesh_4.add_vertex(TriMesh.Point(float(vertices[i][0]),float(vertices[i][1]),float(vertices[i][2]))) for i in range(numver)]
    for i in range(numfac):
        original_mesh_4.add_face(vhs[faces[i][0]], vhs[faces[i][1]], vhs[faces[i][2]])
    # Add mid points to each face
    add_mid_point(original_mesh_4)
    add_pressure(original_mesh_4)
    
    
    # Original mesh 5
    original_mesh_5    = TriMesh()
    vhs                = [original_mesh_5.add_vertex(TriMesh.Point(float(vertices[i][0]),float(vertices[i][1]),float(vertices[i][2]))) for i in range(numver)]
    for i in range(numfac):
        original_mesh_5.add_face(vhs[faces[i][0]], vhs[faces[i][1]], vhs[faces[i][2]])
    # Add mid points to each face
    add_mid_point(original_mesh_5)
    add_pressure(original_mesh_5)
    
    
    # Original mesh 6
    original_mesh_6    = TriMesh()
    vhs                = [original_mesh_6.add_vertex(TriMesh.Point(float(vertices[i][0]),float(vertices[i][1]),float(vertices[i][2]))) for i in range(numver)]
    for i in range(numfac):
        original_mesh_6.add_face(vhs[faces[i][0]], vhs[faces[i][1]], vhs[faces[i][2]])
    # Add mid points to each face
    add_mid_point(original_mesh_6)
    add_pressure(original_mesh_6)
    
    
    # Original mesh 7
    original_mesh_7    = TriMesh()
    vhs                = [original_mesh_7.add_vertex(TriMesh.Point(float(vertices[i][0]),float(vertices[i][1]),float(vertices[i][2]))) for i in range(numver)]
    for i in range(numfac):
        original_mesh_7.add_face(vhs[faces[i][0]], vhs[faces[i][1]], vhs[faces[i][2]])
    # Add mid points to each face
    add_mid_point(original_mesh_7)
    add_pressure(original_mesh_7)
    
    
    # Original mesh 8
    original_mesh_8    = TriMesh()
    vhs                = [original_mesh_8.add_vertex(TriMesh.Point(float(vertices[i][0]),float(vertices[i][1]),float(vertices[i][2]))) for i in range(numver)]
    for i in range(numfac):
        original_mesh_8.add_face(vhs[faces[i][0]], vhs[faces[i][1]], vhs[faces[i][2]])
    # Add mid points to each face
    add_mid_point(original_mesh_8)
    add_pressure(original_mesh_8)

    equilibriate_vertices(original_mesh_0,0.05,10)
    equilibriate_vertices(original_mesh_1,0.2,10)
    equilibriate_vertices(original_mesh_2,0.5,10)
    equilibriate_vertices(original_mesh_3,0.05,50)
    equilibriate_vertices(original_mesh_4,0.2,50)
    equilibriate_vertices(original_mesh_5,0.5,50)
    equilibriate_vertices(original_mesh_6,0.05,100)
    equilibriate_vertices(original_mesh_7,0.2,100)
    equilibriate_vertices(original_mesh_8,0.5,100)
    
    draw_primary_mesh(foam_mesh,title='primary_mesh_fixed_original')
    draw_primary_mesh(original_mesh_0,title='primary_mesh_fixed_05_10')
    draw_primary_mesh(original_mesh_1,title='primary_mesh_fixed_2_10')
    draw_primary_mesh(original_mesh_2,title='primary_mesh_fixed_5_10')
    draw_primary_mesh(original_mesh_3,title='primary_mesh_fixed_05_50')
    draw_primary_mesh(original_mesh_4,title='primary_mesh_fixed_2_50')
    draw_primary_mesh(original_mesh_5,title='primary_mesh_fixed_5_50')
    draw_primary_mesh(original_mesh_6,title='primary_mesh_fixed_05_100')
    draw_primary_mesh(original_mesh_7,title='primary_mesh_fixed_2_100')
    draw_primary_mesh(original_mesh_8,title='primary_mesh_fixed_5_100')


