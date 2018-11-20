from foam_import_openmesh import *
from foam_equilibriate import *
from foam_properties import *
from foam_visualise import *
from foam_topology import *
from foam_helper_functions import *

def init_grid_rand(min, max, num_elem, rectangular=True):
    """ Initializes square grid, randomly distributed
        double min: minimum value for both x and y in grid
        double max: maximum value for both x and y in grid
        num_elem int: number of points
        rectangular bool: circular or regtangular
        return: Grid (in 1D) for x and y, respectively
        """

    if num_elem not in locals():
        num_elem = int(np.round(10 * (max - min) ** 2))

    # Sample from rectangular distribution
    if rectangular:
        xv     = np.random.uniform(min, max, size=(1,num_elem))[0]
        yv     = np.random.uniform(min, max, size=(1,num_elem))[0]

    # Sample from circular distribution
    else:
        cntr = (max-min)/2.0                                                    # Center of circle to sample from
        rad  = cntr                                                             # Radius of that circle (equal in this case)
        xv   = np.zeros(num_elem)
        yv   = np.zeros(num_elem)
        for i in range(num_elem):
            u, theta = [(random.random()+random.random())*rad, 2*math.pi*random.random()]
            if u<rad:
                r = u
            else:
                r = 2*rad - u
            xv[i]    = cntr + r * math.cos(theta)
            yv[i]    = cntr + r * math.sin(theta)

    return xv,yv


def create_random_mesh(min,max,rectangular=True):
    """ Create a random mesh in OpenMesh using Matplotlib triangulation
        int min: lower boundary of window
        int max: upper boundary of window
        rectangular bool: circular or regtangular
        return mesh: OpenMesh mesh structure
    """

    xv, yv       = init_grid_rand(min,max,5,rectangular)                        # Initialize random grid
    triang       = tri.Triangulation(xv, yv)                                    # Do Delauney-triangulation on the set of x- and y-coordinates
    faces        = triang.get_masked_triangles()                                # Extracts faces from the Delaunay-triangulation from the matplotlib.tri-package
    num_ver      = len(xv)
    num_fac      = len(faces)
    mesh         = TriMesh()
    vhs          = [mesh.add_vertex(TriMesh.Point(float(xv[i]),float(yv[i]),float(0))) for i in range(num_ver)]

    for i in range(num_fac):
        mesh.add_face(vhs[faces[i][0]], vhs[faces[i][1]], vhs[faces[i][2]])

    # Add mid points to each face
    add_mid_point(mesh)

    return mesh




