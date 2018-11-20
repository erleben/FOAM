import pymesh
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
import math
import config
import collections
from timeit import *


def init_grid_even(min, max, factor=1):
    """ Initializes square grid, evenly distributed
    double min: minimum value for both x and y in grid
    double max: maximum value for both x and y in grid
    factor int: increases number of points in grid if > 1
    return: Grid (in 1D) for x and y, respectively
    """
    nsteps   = (max - min + 1) * factor
    x        = np.linspace(min, max, nsteps)
    y        = np.linspace(min, max, nsteps)
    xv0, yv0 = np.meshgrid(x, y)
    return np.ravel(xv0), np.ravel(yv0)


def init_grid_rand(min, max, factor=1):
    """ Initializes square grid, randomly distributed
        double min: minimum value for both x and y in grid
        double max: maximum value for both x and y in grid
        factor int: increases number of points in grid if > 1
        return: Grid (in 1D) for x and y, respectively
        """
    nsteps   = max - min + 1                                                    # Number of steps in grid
    len      = nsteps**2*factor                                                 # Number of elements in grid
    xv       = np.random.uniform(min, max, size=(1,len))
    yv       = np.random.uniform(min, max, size=(1,len))
    return xv[0],yv[0]


def plot_del(xv,yv,title='Title'):
    """ Plot the Delauney triangulation of a given set of x- and y-coordinates using matplotlib.tri
        double array xv: 1D array of x-coordinates
        double array yv: 1D array of corresponding y-coordinates
        title string: title of plot
    """
    triang = tri.Triangulation(xv, yv)                                          # Do Delauney-triangulation on the set of x- and y-coordinates

    plt.figure()
    plt.triplot(triang, 'bo-')
    plt.title(title)
    plt.show()


def neighbouring_bubbles(vertices,faces):
    """ Count and name the number of neighbors for each vertex (bubble in foam)
        double array vertices: 2D array of vertex coordinates
        double array faces: 2D array of faces and which vertices they border
        return array: a 2D-array of each vertex (bubble) and which neighbors it has
    """
    v_len  = len(vertices)
    f_len  = len(faces)
    neigh  = []
    
    for i in range(v_len):
        corners = []
        for j in range(f_len):
            edge = faces[j]
            if i in edge:
                corners_temp = [x for x in edge if x!=i]
                corners = corners + corners_temp
        temp = np.unique(corners)
        neigh.append(temp.tolist())

    return neigh



if __name__ == "__main__":
    # Read variables from config-file (Initially for random grid)
    min = config.length['min']
    max = config.length['max']
    
    
    # Create points, form mesh and plot using matplotlib.tri
    xv_e, yv_e   = init_grid_even(min,max)                                      # Initialize even grid
    xv_r, yv_r   = init_grid_rand(min,max,5)                                    # Initialize random grid

    plot_del(xv_e,yv_e,"Grid")
    plot_del(xv_r,yv_r,"Random")

    ## Create points and form mesh using PyMesh and Matplotlib libraries
    # Even
    # Matplotlib
    numver_e     = len(xv_e)
    tic()
    triang_e     = tri.Triangulation(xv_e, yv_e)                                # Do Delauney-triangulation on the set of x- and y-coordinates
    toc()
    faces_e      = triang_e.get_masked_triangles()                              # Extracts faces from the Delaunay-triangulation from the matplotlib.tri-package
    vertices_e   = np.array([[xv_e[i], yv_e[i], 0] for i in range(numver_e)])   # Create 2D array containing vertex-coordinates
    
    # PyMesh
    segments     = np.array([[0,10],[10,120],[120,110],[110,0]])                # Indicies of bounding corner segments
    tic()
    pymesh_trian = pymesh.triangulate(vertices_e,segments,1)                    # This is the PyMesh-way of doing Delaunay-triangulation. Returns vertices and faces
    toc()
    mesh_e       = pymesh.form_mesh(pymesh_trian[0], pymesh_trian[1])           # Form mesh using PyMesh and vertices and faces
    
    # Random
    numver_r     = len(xv_r)
    triang_r     = tri.Triangulation(xv_r, yv_r)                                # Do Delauney-triangulation on the set of x- and y-coordinates
    faces_r      = triang_r.get_masked_triangles()                              # Extracts faces from the Delaunay-triangulation from the matplotlib.tri-package
    vertices_r   = np.array([[xv_r[i], yv_r[i], 0] for i in range(numver_r)])   # Create 2D array containing vertex-coordinates
    mesh_r       = pymesh.form_mesh(vertices_r, faces_r)                        # Form mesh using PyMesh and vertices and faces








