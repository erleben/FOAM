from foam_import_openmesh import *
from foam_properties import *


def retriangulate_mesh(mesh):
    """ Retriangulate mesh and assign mid points afterwards
        TriMesh mesh: an existing OpenMesh mesh
    """

    # Delete all faces in mesh
    mesh.request_face_status()
    mesh.request_edge_status()
    mesh.request_halfedge_status()
    for fh in mesh.faces(): mesh.delete_face(fh, False)
    mesh.garbage_collection()

    # Retriangulate
    vh_list      = [vh for vh in mesh.vertices()]                                           # Fetch vertex coordinates
    vertices_tmp = [[mesh.point(vh)[0], mesh.point(vh)[1], 0.0] for vh in vh_list]          # Fetch vertex coordinates from mesh
    vertices     = [list(vh) for vh in zip(*vertices_tmp)]                                  # Tranpose matrix
    triang       = tri.Triangulation(vertices[0], vertices[1])                              # Do Delaunay triangulation
    faces        = triang.get_masked_triangles()                                            # Extract faces
    num_faces    = len(faces)
    for i in range(num_faces): mesh.add_face(vh_list[faces[i][0]], vh_list[faces[i][1]], vh_list[faces[i][2]])

    # Add mid point property
    add_mid_point(mesh)


def average_equilibriate_vertices(mesh, scale=0.1, max_iterations=10, fix_boundaries=True, sliver=True):
    """ Move vertices using average positions of 1-ring neighborhoods to create more realistic mesh
        TriMesh mesh: an existing OpenMesh mesh
        int scale: scale of movement
        int max_iterations: maximum number of iterations
        bool fix_boundaries: fix boundary vertices
    """

    # Add property to each vertex in order to treat each vertex using information from the same mesh
    prop_handle = VPropHandle()
    mesh.add_property(prop_handle, "vprop_float")

    if sliver: remove_sliver_faces(mesh)

    # Move vertices
    for _ in range(max_iterations):
        # Calculate destination
        for vh in mesh.vertices():
            if mesh.is_boundary(vh)==False or fix_boundaries==False:

                mx = 0
                my = 0
                nn = 0
                for n in mesh.vv(vh):                                           # Sum over neighboring coordinates in x- and y-directions
                    mx += mesh.point(n)[0]
                    my += mesh.point(n)[1]
                    nn += 1
                mx /= nn
                my /= nn
                mesh.set_property(prop_handle, vh, TriMesh.Point(mx, my, 0))

        # Move to destination
        for vh in mesh.vertices():
            if mesh.is_boundary(vh)==False or fix_boundaries==False:
                p  = mesh.point(vh)
                m  = mesh.property(prop_handle,vh)
                new_p = p + (m-p)*scale
                mesh.set_point(vh, new_p)

    # Update properties
    update_mid_point(mesh)


def spring_equilibriate_vertices(mesh, k=0.5, l_0=0.36, delta_t=0.05, max_iterations=10):
    """ Move vertices using string approximations between vertices to create more realistic mesh
        TriMesh mesh: an existing OpenMesh mesh
        int k: spring constant
        int l_0: equilibrium length of each spring
        int delta_t: time step size
        int max_iterations: maximum number of iterations
    """

    global spring_handle
    spring_handle = VPropHandle()
    mesh.add_property(spring_handle, "vprop_float")

    for _ in range(max_iterations):
        # Calculate destination
        for vh0 in mesh.vertices():
            v0   = mesh.point(vh0)
            F    = [0, 0]                                                       # Force acting on v0
            for vh in mesh.vv(vh0):
                v     = mesh.point(vh)
                dist  = np.sqrt((v[0]-v0[0])**2+(v[1]-v0[1])**2)                # Distance between vertices
                F_i   = k*(dist-l_0)                                            # Magnitude of force, calculated using Hooke's law
                F_hat = [(v[0]-v0[0])/dist, (v[1]-v0[1])/dist]                  # Unit vector of acting force
                F     = TriMesh.Point(F[0]+F_hat[0]*F_i,F[1]+F_hat[1]*F_i,0)    # Acting force
            mesh.set_property(spring_handle,vh0,F)

        for vh in mesh.vertices():
            F_i   = mesh.property(spring_handle,vh)
            p     = mesh.point(vh)
            new_p = p + delta_t*F_i
            mesh.set_point(vh, new_p)

    # Update properties
    update_mid_point(mesh)


def combine_equilibriate_vertices(mesh, max_iterations=1):
    """ Combine average_equilibriate_vertices() and spring_equilibriate_vertices()
        TriMesh mesh: an existing OpenMesh mesh
        int max_iterations: maximum number of iterations
    """

    remove_sliver_faces(mesh)                                                   # Remove sliver faces on boundary
    remove_shark_fin_faces(mesh)
    for i in range(max_iterations):
        spring_equilibriate_vertices(mesh)                                      # Do spring equilibration
        average_equilibriate_vertices(mesh,sliver=False)                        # Do average position of 1-ring neighborhood equilibration
        retriangulate_mesh(mesh)
        remove_sliver_faces(mesh)                                               # Remove sliver faces on boundary

    remove_shark_fin_faces(mesh)


def remove_shark_fin_faces(mesh):
    """ Removes boundary faces with only one neighboring face
        TriMesh mesh: an existing OpenMesh mesh
    """

    # MEMORY REQUESTS
    mesh.request_face_status()
    mesh.request_vertex_status()
    mesh.request_halfedge_status()
    mesh.request_edge_status()

    # VARIABLES
    preceed      = True                                                             # False, when there are no more shark fins in mesh
    delete_count = 0                                                                # Counter for each iteration through mesh

    while preceed:
        for fh in mesh.faces():                                                     # Iterate over all faces in mesh
            if mesh.is_boundary(fh):
                neigh = 0                                                           # Neighbor counter
                for __ in mesh.ff(fh):
                    neigh += 1

                if neigh < 2:
                    mesh.delete_face(fh, True)
                    delete_count += 1

        mesh.garbage_collection()

        if delete_count == 0:
            preceed = False

        delete_count = 0                                                            # Reset counter


def remove_sliver_faces(mesh, max_angle=2.50):
    """ Removes boundary faces that are very sliver as a result of the Delaunay triangulation
        TriMesh mesh: an existing OpenMesh mesh
        int max_angle: maximum allowed angle on border
    """

    # Request status. Otherwise, Python will crash when calling for a delete operation
    mesh.request_face_status()
    mesh.request_vertex_status()
    mesh.request_halfedge_status()
    mesh.request_edge_status()

    sliver = True

    while sliver:                                                                   # Preceed, while there are still sliver faces on boundary
        delete_counter = 0
        for fh in mesh.faces():                                                     # Iterate over all faces in mesh
            if mesh.is_boundary(fh):
                vertices = [mesh.point(vh) for vh in mesh.fv(fh)]                   # Iterate over the vertices in a face

                # Calculate lengths of sides in triangular face
                a       = np.sqrt((vertices[1][0]-vertices[0][0])**2 + (vertices[1][1]-vertices[0][1])**2)
                b       = np.sqrt((vertices[2][0]-vertices[0][0])**2 + (vertices[2][1]-vertices[0][1])**2)
                c       = np.sqrt((vertices[2][0]-vertices[1][0])**2 + (vertices[2][1]-vertices[1][1])**2)

                theta_a = np.arccos((b**2 + c**2 - a**2) / (2*b*c))
                theta_b = np.arccos((a**2 + c**2 - b**2) / (2*a*c))
                theta_c = np.arccos((a**2 + b**2 - c**2) / (2*a*b))

                if (theta_a > max_angle) | (theta_b > max_angle) | (theta_c > max_angle):
                    mesh.delete_face(fh)
                    delete_counter += 1

        if delete_counter == 0: sliver = False                                      # Terminate, if no faces were deleted in last round

        # Do garbage collection to delete deleted objects from memory
        mesh.garbage_collection()

        # Update pressure values on boundary
        update_boundary_pressure(mesh)


def reduce_boundary_facing_faces(mesh):
    """ Reduce the number of world-facing faces, if more are adjacent to each other
        TriMesh mesh: an existing OpenMesh mesh
    """

    # Variables
    preceed     = True
    mesh.request_edge_status()

    # Preceed, while there are still adjacent world facing faces
    while preceed:

        # Find halfedge on boundary to iterate over
        for heh in mesh.halfedges():
            if mesh.is_boundary(heh):
                heh_original  = heh
                vh_original   = mesh.to_vertex_handle(heh_original)
                coor_original = (mesh.point(vh_original)[0], mesh.point(vh_original)[1])
                break

        # Initialise
        coor_iter          = 0                                                          # Initialise coor_iter as being different from coor_original
        heh_iter           = heh_original                                               # Initialise halfedge iterator
        to_be_flipped_list = []                                                         # Initialise list of edges to be flipped
        iter               = 0

        # Count number of boundary edges
        number_of_edges_on_boundary = np.sum(np.array([1 for eh1 in mesh.edges() if mesh.is_boundary(eh1)]))

        # Do edge flip, if iterator comes across more world facing faces (characterised by having a vertex with only 3 connecting edges)
        while coor_iter != coor_original:
            iter         += 1
            if iter > number_of_edges_on_boundary: break
            heh_iter     = mesh.next_halfedge_handle(heh_iter)
            vh_0         = mesh.to_vertex_handle(heh_iter)
            coor_iter    = (mesh.point(vh_0)[0], mesh.point(vh_0)[1])
            edge_counter = 0

            # Count number of incoming edges
            for __ in mesh.ve(vh_0):
                edge_counter += 1

            # Write to to_be_flipped_list, if there are 3 incoming edges
            if edge_counter == 3:
                for veh in mesh.ve(vh_0):
                    if mesh.is_boundary(veh) == False:
                        to_be_flipped_list.append(veh)

                # Move to next halfedge on boundary
                heh_iter = mesh.next_halfedge_handle(heh_iter)

        # Stop, if no flips were performed in last round
        if len(to_be_flipped_list) == 0:
            preceed = False

        # Perform flips
        for i in range(len(to_be_flipped_list)):
            mesh.flip(to_be_flipped_list[i])
            mesh.garbage_collection()


        # Remove shark fin faces and sliver faces that are results of the edge flips
        remove_shark_fin_faces(mesh)

    # Update properties
    #add_boundary_neighbour_faces(mesh)                                             # Update boundary face neighbour face property
    update_boundary_pressure(mesh)                                                  # Update pressure values on boundary
    #add_mid_point(mesh)                                                             # Update mid point face property









