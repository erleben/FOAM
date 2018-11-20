from foam_import_openmesh import *
import config

# ----------------------------------------------------------------------------------------------------------------------
# Import constants
# ----------------------------------------------------------------------------------------------------------------------
p_imp        = config.initialise['bubble_pressure']
wp_imp       = config.initialise['world_pressure']
gamma_imp    = config.physical['gamma']
kappa_imp    = config.physical['kappa']
delta_t_imp  = config.numerical['delta_t']
h_imp        = config.numerical['h']
draw_dh_imp  = config.test['draw_dh']
idx_draw_imp = config.test['idx_draw']
scale_imp    = config.draw['scale']


# ----------------------------------------------------------------------------------------------------------------------
# Add computational properties
# ----------------------------------------------------------------------------------------------------------------------

def add_mid_point(mesh):
    """ Add mid point property to each face in existing mesh
        TriMesh mesh: an existing OpenMesh mesh
    """

    # Define and make property handle global, so that it is accessible from other functions
    global mid_point_handle
    mid_point_handle = FPropHandle()
    mesh.add_property(mid_point_handle, "fprop_float")

    for fh in mesh.faces():
        cog = TriMesh.Point(0, 0, 0)
        for vh in mesh.fv(fh):
            cog += mesh.point(vh)
        mesh.set_property(mid_point_handle, fh, cog/3)


def update_mid_point(mesh):
    """ Update mid point property to each face in existing mesh
        TriMesh mesh: an existing OpenMesh mesh
    """

    for fh in mesh.faces():
        cog = TriMesh.Point(0, 0, 0)
        for vh in mesh.fv(fh):
            cog += mesh.point(vh)
        mesh.set_property(mid_point_handle, fh, cog/3)


def update_vertex_positions(mesh):
    """ Update vertex positions in simulation, when junctions have been moved. Make the vertex coordinates the centre of the bubbles.
        TriMesh mesh: an existing OpenMesh mesh
    """

    for vh in mesh.vertices():
        if not mesh.is_boundary(vh):
            # Initialise variables
            cog = TriMesh.Point(0, 0, 0)
            neighbours = 0

            for fh0 in mesh.vf(vh):
                if mesh.is_boundary(fh0) == False:
                    neighbours += 1
                    cog += mesh.property(mid_point_handle, fh0)

            mesh.set_point(vh, cog / neighbours)


def add_pressure(mesh, equal = True, approx_pressure = p_imp, world_pressure = wp_imp):
    """ Add pressure property to an existing mesh
        TriMesh mesh: an existing OpenMesh mesh
        bool equal: True if equal pressure value is wanted, False if different pressure values are wanted
        double approx_pressure: value for pressure, if equal pressure across mesh is wanted. If not, it is used as offset for randomisation
        double world_pressure_fraction: fraction, for how much lower the world pressure is than the approximate pressure inside the foam
    """

    # Define and make property handle global, so that it is accessible from other functions
    global pressure_handle
    pressure_handle = VPropHandle()
    mesh.add_property(pressure_handle, "vprop_float")

    for vh in mesh.vertices():
        if equal:                                                               # If equal pressure is desired
            if mesh.is_boundary(vh):
                mesh.set_property(pressure_handle, vh, world_pressure)
            else:
                mesh.set_property(pressure_handle, vh, approx_pressure)
        else:                                                                   # If random pressure is desired
            if mesh.is_boundary(vh):
                mesh.set_property(pressure_handle, vh, world_pressure)
            else:
                neighbors = 0
                for vh1 in mesh.vv(vh):                                         # Iterate over neighboring vertices and count them
                    neighbors += 1
                pressure = approx_pressure + 0.2/neighbors                      # MAGIC FUNCTION!!!!!!!!!!!!!!!!!!!!!!
                mesh.set_property(pressure_handle, vh, pressure)


def update_boundary_pressure(mesh, world_pressure = wp_imp):
    """ Update boundaries existing mesh
        TriMesh mesh: an existing OpenMesh mesh
        bool equal: True if equal pressure value is wanted, False if different pressure values are wanted
        double pressure: value for pressure, if equal pressure across mesh is wanted
    """

    for vh in mesh.vertices():
        if mesh.is_boundary(vh):
            mesh.set_property(pressure_handle, vh, world_pressure)


def add_face_connectivity(mesh):
    """ Add face connectivity. Each junction containing face will have three bubble film connected faces
        TriMesh mesh: an existing OpenMesh mesh
    """

    # Define and make property handle global, so that it is accessible from other functions
    global face_connectivity_handle
    face_connectivity_handle = FPropHandle()
    mesh.add_property(face_connectivity_handle, "fprop_float")

    for fh in mesh.faces():
        if not mesh.is_boundary(fh):
            is_it_boundary = sum([mesh.is_boundary(vh0) for vh0 in mesh.fv(fh)])                                        # Check, if face belongs to junction on boundary of foam

            if is_it_boundary == 0:
                # Create list of edge neighboring faces to fh
                neighbour_faces = np.array([fh0 for fh0 in mesh.ff(fh)])
            elif is_it_boundary == 1:
                neighbour_faces = np.array([])
                for fh0 in mesh.ff(fh):
                    if not mesh.is_boundary(fh0):
                        neighbour_faces = np.append(neighbour_faces, fh0)
                    else:
                        # Visit neighbouring faces to find the other coordinate
                        fh_next = fh0
                        fh_start = fh

                        while mesh.is_boundary(fh_next):
                            for fh_iter in mesh.ff(fh_next):
                                if fh_iter.idx() != fh_start.idx():
                                    fh_start = fh_next                                                                  # Move one step forward
                                    fh_next = fh_iter                                                                   # Move one step forward
                                    break

                        neighbour_faces = np.append(neighbour_faces, fh_next)

            # Properties of current face
            vertices_0 = np.array([[vh0.idx(), vh0, mesh.property(pressure_handle, vh0), mesh.is_boundary(vh0)] for vh0 in mesh.fv(fh)])   # List of vertices and pressures in fh
            x_0        = np.array([mesh.property(mid_point_handle, fh)[0], mesh.property(mid_point_handle, fh)[1]])                        # Junction coordinate using notation from thesis

            # Initialise faces
            fh_ij = neighbour_faces[0]
            fh_ik = neighbour_faces[1]
            fh_jk = neighbour_faces[2]

            # ij. See thesis for an explanation of indices. 0 is the coordinate of fh and i, j, and k describe adjacent bubble cells. Thus, ij refers to the junction on the other end of bubble film wrt 0 between i and j.
            x_ij        = np.array([mesh.property(mid_point_handle, fh_ij)[0], mesh.property(mid_point_handle, fh_ij)[1]])
            vertices_ij = np.array([[vh0.idx(), vh0, mesh.is_boundary(vh0)] for vh0 in mesh.fv(fh_ij)])

            # ik
            x_ik        = np.array([mesh.property(mid_point_handle, fh_ik)[0], mesh.property(mid_point_handle, fh_ik)[1]])
            vertices_ik = np.array([[vh0.idx(), vh0, mesh.is_boundary(vh0)] for vh0 in mesh.fv(fh_ik)])

            # jk
            x_jk        = np.array([mesh.property(mid_point_handle, fh_jk)[0], mesh.property(mid_point_handle, fh_jk)[1]])

            # Indices for vertices in fh in correct order
            if is_it_boundary == 0:
                idx_pi = np.in1d(vertices_0[:, 0], vertices_ij[:, 0]) & np.in1d(vertices_0[:, 0], vertices_ik[:, 0])    # Common vertex of ij and ik, i.e. i
                idx_pj = np.in1d(vertices_0[:, 0], vertices_ij[:, 0]) & np.invert(np.in1d(vertices_0[:, 0], vertices_ik[:, 0]))
                idx_pk = np.in1d(vertices_0[:, 0], vertices_ik[:, 0]) & np.invert(np.in1d(vertices_0[:, 0], vertices_ij[:, 0]))
            elif is_it_boundary:
                is_boundary_ij = np.any(vertices_ij[:, 2])
                is_boundary_ik = np.any(vertices_ik[:, 2])

                # Different boundary scenarios
                if is_boundary_ij and is_boundary_ik:
                    if fh_ij.idx() == fh_ik.idx():                                                                        # Treating special case with isolated bubble on boundary
                        idx_temp = np.in1d(vertices_0[:, 0],vertices_ij[:, 0])                                          # = np.in1d(v_i,v_k), since f_ij=f_ik
                        vertices_temp = vertices_0[idx_temp][:, 1]
                        neigh_ver = np.array([not (sum([not (mesh.is_boundary(vh)) for vh in mesh.vv(vh0)]) > 1) for vh0 in vertices_temp])  # List of: "is vertex a 2-bubble vertex?"
                        idx_i = vertices_temp[neigh_ver][0].idx()                                                       # Vertex index of p_i in vertices_0
                        idx_j = vertices_temp[np.invert(neigh_ver)][0].idx()
                        idx_pi = np.in1d(vertices_0[:, 0], idx_i)
                        idx_pj = np.in1d(vertices_0[:, 0], idx_j)
                        idx_pk = np.invert(idx_pi | idx_pj)
                    else:
                        idx_temp = vertices_0[:, 3]
                        idx_pi = np.array([x for x in idx_temp])  # Weird work around. Otherwise, np.invert does not work
                        idx_pj = np.in1d(vertices_0[:, 0], vertices_ij[:, 0]) & np.invert(idx_pi)
                        idx_pk = np.in1d(vertices_0[:, 0], vertices_ik[:, 0]) & np.invert(idx_pi)
                elif is_boundary_ij and not (is_boundary_ik):
                    idx_pi = np.in1d(vertices_0[:, 0], vertices_ij[:, 0]) & np.in1d(vertices_0[:, 0], vertices_ik[:, 0])
                    idx_pk = np.in1d(vertices_0[:, 0], vertices_ik[:, 0]) & np.invert(np.in1d(vertices_0[:, 0], vertices_ij[:, 0]))
                    idx_pj = np.invert(idx_pi | idx_pk)
                elif is_boundary_ik and not (is_boundary_ij):
                    idx_pi = np.in1d(vertices_0[:, 0], vertices_ij[:, 0]) & np.in1d(vertices_0[:, 0],vertices_ik[:, 0])
                    idx_pj = np.in1d(vertices_0[:, 0], vertices_ij[:, 0]) & np.invert(np.in1d(vertices_0[:, 0], vertices_ik[:, 0]))
                    idx_pk = np.invert(idx_pi | idx_pj)

            # Associated vertex index
            idx_i = vertices_0[idx_pi, 0][0]
            idx_j = vertices_0[idx_pj, 0][0]
            idx_k = vertices_0[idx_pk, 0][0]

            # Vertices
            # Associated vertex index
            vh_i = vertices_0[idx_pi, 1][0]
            vh_j = vertices_0[idx_pj, 1][0]
            vh_k = vertices_0[idx_pk, 1][0]

            # Pressure values
            p_i = vertices_0[idx_pi, 2][0]
            p_j = vertices_0[idx_pj, 2][0]
            p_k = vertices_0[idx_pk, 2][0]

            # Order, so that 'k' is on boundary, if boundary junction
            if is_it_boundary and mesh.is_boundary(vh_k) == False:
                if mesh.is_boundary(vh_i):
                    idx_temp = idx_k
                    idx_k    = idx_i
                    idx_i    = idx_j
                    idx_j    = idx_temp

                    fh_temp = fh_ij
                    fh_ij   = fh_jk
                    fh_jk   = fh_ik
                    fh_ik   = fh_temp

                    vh_temp = vh_k
                    vh_k    = vh_i
                    vh_i    = vh_j
                    vh_j    = vh_temp
                else:
                    idx_temp = idx_k
                    idx_k    = idx_j
                    idx_j    = idx_i
                    idx_i    = idx_temp

                    fh_temp = fh_ij
                    fh_ij   = fh_ik
                    fh_ik   = fh_jk
                    fh_jk   = fh_temp

                    vh_temp = vh_k
                    vh_k    = vh_j
                    vh_j    = vh_i
                    vh_i    = vh_temp


            # Set attribute vector
            # attributes = [idx_i, idx_j, idx_k, x_0ij[0], x_0ij[1], x_0ik[0], x_0ik[1], x_0jk[0], x_0jk[1], p_i, p_j, p_k]
            #attributes = [idx_i, idx_j, idx_k, x_ij[0], x_ij[1], x_ik[0], x_ik[1], x_jk[0], x_jk[1], p_i, p_j, p_k]
            attributes = [idx_i, idx_j, idx_k, fh_ij, fh_ik, fh_jk, vh_i, vh_j, vh_k]

            # Assign as property
            mesh.set_property(face_connectivity_handle, fh, attributes)


# ----------------------------------------------------------------------------------------------------------------------
# Add physical properties
# ----------------------------------------------------------------------------------------------------------------------

def prop_film_lengths(mesh, gamma = gamma_imp):
    """ Calculate film lengths of all films in mesh. Add as property to each film-belonging edge in dual mesh.
        TriMesh mesh: an existing OpenMesh mesh
    """

    # Define and make property handle global, so that it is accessible from other functions
    global film_length_handle                                                           # Define and make property handle global, so that it is accessible from other functions
    film_length_handle = EPropHandle()
    mesh.add_property(film_length_handle, "eprop_float")

    for eh in mesh.edges():
        if mesh.is_boundary(eh) == False:
            f0_h = mesh.face_handle(mesh.halfedge_handle(eh, 0))                        # Face handle for the one neighboring face
            f1_h = mesh.face_handle(mesh.halfedge_handle(eh, 1))                        # Face handle for the other neighboring face
            v0_h = mesh.to_vertex_handle(mesh.halfedge_handle(eh, 0))                   # Vertex handle for the vertex at the end
            v1_h = mesh.to_vertex_handle(mesh.halfedge_handle(eh, 1))                   # Vertex handle for other vertex at the end

            is_double_boundary_film = mesh.is_boundary(f0_h) != mesh.is_boundary(f1_h)
            is_redundant_edge       = mesh.is_boundary(f0_h) and mesh.is_boundary(f1_h)

            if is_redundant_edge==False:                   # Start by treating nice and easy edges inside foam
                # Handle case of boundary film edge in dual mesh
                if is_double_boundary_film:
                    fh_next        = f0_h if mesh.is_boundary(f0_h) else f1_h
                    fh_start       = f0_h if not mesh.is_boundary(f0_h) else f1_h
                    coor_flm_end_0 = np.array([mesh.property(mid_point_handle, fh_start)[0], mesh.property(mid_point_handle, fh_start)[1]])

                    # Visit neighbouring faces to find the other coordinate
                    while mesh.is_boundary(fh_next):
                        for fh in mesh.ff(fh_next):
                            if fh != fh_start:
                                fh_start = fh_next  # Move one step forward
                                fh_next = fh  # Move one step forward
                                break
                    coor_flm_end_1 = np.array([mesh.property(mid_point_handle, fh_next)[0], mesh.property(mid_point_handle, fh_next)[1]])

                # Handle normal case
                else:
                    # Calculate the two edge vectors spanning the part of the bubble film (using the one or the other bubble centrum does not matter
                    coor_flm_end_0 = np.array([mesh.property(mid_point_handle, f0_h)[0], mesh.property(mid_point_handle, f0_h)[1]])
                    coor_flm_end_1 = np.array([mesh.property(mid_point_handle, f1_h)[0], mesh.property(mid_point_handle, f1_h)[1]])

                d              = np.sqrt(np.sum((coor_flm_end_1-coor_flm_end_0)**2))    # Calculate distance between two end points of bubble film
                delta_p        = mesh.property(pressure_handle, v1_h) - mesh.property(pressure_handle, v0_h)
                film_length    = calc_film_length(coor_flm_end_0, coor_flm_end_1, delta_p)

                mesh.set_property(film_length_handle, eh, film_length)
            else:
                pass


def prop_plateau_angles(mesh, gamma = gamma_imp, fh0 = False):
    """ Calculate angles of plateau junctions in mesh. Add as property to each non-boundary face in dual mesh.
        TriMesh mesh: an existing OpenMesh mesh
        TriMesh FaceHandle fh0: Specific face handle, if only the angles of one face are desired
    """

    # HANDLES
    global plateau_angle_handle                                                        # Define and make property handle global, so that it is accessible from other functions
    plateau_angle_handle = FPropHandle()
    mesh.add_property(plateau_angle_handle, "fprop_float")

    for fh in mesh.faces():
        if not mesh.is_boundary(fh):
            thetas = calc_plateau_angles(mesh, fh)

            # Add list of angles as property to current face
            mesh.set_property(plateau_angle_handle, fh, thetas)


def prop_areas(mesh, gamma = gamma_imp, kappa = kappa_imp, delta_t = delta_t_imp):
    """ Calculate areas of bubbles in mesh. Add as property to each bubble-belonging vertex in dual mesh. Do the same with the target area of the bubble
        TriMesh mesh: an existing OpenMesh mesh
        double kappa: constant of permeability
        double delta_t: size of time step in simulation
    """

    # HANDLES
    global area_handle  # Define and make property handle global, so that it is accessible from other functions
    area_handle = VPropHandle()
    mesh.add_property(area_handle, "vprop_float")

    global target_area_handle  # Define and make property handle global, so that it is accessible from other functions
    target_area_handle = VPropHandle()
    mesh.add_property(target_area_handle, "vprop_float")

    # LOOP THROUGH MESH
    for vh in mesh.vertices():
        if not mesh.is_boundary(vh):
            areas       = calc_areas_bubble(mesh, vh, True, gamma_imp, kappa, delta_t)
            area        = areas[0]
            area_target = areas[1]

            mesh.set_property(area_handle, vh, area)
            mesh.set_property(target_area_handle, vh, area_target)


def prop_jacobian(mesh, h = h_imp, draw = False):
    """ Calculate jacobian for junction in mesh. Add as property to each cell vertex in mesh.
        TriMesh mesh: an existing OpenMesh mesh
        double h: describing step distance in approximation of derivative
    """

    # HANDLES
    global jacobian_handle  # Define and make property handle global, so that it is accessible from other functions
    jacobian_handle = FPropHandle()
    mesh.add_property(jacobian_handle, "fprop_float")

    berlin = 0

    for fh in mesh.faces():
        if not mesh.is_boundary(fh):
            ## FETCH CONSTANTS, VARIABLES ETC.
            # Check, if junction is on boundary
            is_it_boundary = sum([mesh.is_boundary(vh0) for vh0 in mesh.fv(fh)])

            # Fetch attribute values
            x_0 = np.array([mesh.property(mid_point_handle, fh)[0], mesh.property(mid_point_handle, fh)[1]])            # Coordinate of junction
            #idx_i, idx_j, idx_k, x_K_x, x_K_y, x_J_x, x_J_y, x_I_x, x_I_y, p_i, p_j, p_k = mesh.property(face_connectivity_handle, fh)
            idx_i, idx_j, idx_k, fh_K, fh_J, fh_I, vh_i, vh_j, vh_k = mesh.property(face_connectivity_handle, fh)


            # Collect vectors in NumPy-array
            # x_K = np.array([x_K_x, x_K_y])                                                                              # Junction oppposite of K
            # x_J = np.array([x_J_x, x_J_y])                                                                              # Junction oppposite of J
            # x_I = np.array([x_I_x, x_I_y])                                                                              # Junction oppposite of I
            x_K = np.array([mesh.property(mid_point_handle, fh_K)[0], mesh.property(mid_point_handle, fh_K)[1]])          # Junction oppposite of K
            x_J = np.array([mesh.property(mid_point_handle, fh_J)[0], mesh.property(mid_point_handle, fh_J)[1]])          # Junction oppposite of J
            x_I = np.array([mesh.property(mid_point_handle, fh_I)[0], mesh.property(mid_point_handle, fh_I)[1]])          # Junction oppposite of I

            # Coordinates of bubble centres
            # c_i = [np.array([mesh.point(vh)[0], mesh.point(vh)[1]]) for vh in mesh.fv(fh) if vh.idx() == idx_i][0]
            # c_j = [np.array([mesh.point(vh)[0], mesh.point(vh)[1]]) for vh in mesh.fv(fh) if vh.idx() == idx_j][0]
            # c_k = [np.array([mesh.point(vh)[0], mesh.point(vh)[1]]) for vh in mesh.fv(fh) if vh.idx() == idx_k][0]
            c_i = np.array([mesh.point(vh_i)[0], mesh.point(vh_i)[1]])
            c_j = np.array([mesh.point(vh_j)[0], mesh.point(vh_j)[1]])
            c_k = np.array([mesh.point(vh_k)[0], mesh.point(vh_k)[1]])

            # Pressure values
            p_i = mesh.property(pressure_handle, vh_i)
            p_j = mesh.property(pressure_handle, vh_j)
            p_k = mesh.property(pressure_handle, vh_k)

            # Create finite difference vectors
            h_x = np.array([h, 0])
            h_y = np.array([0, h])

            # Central difference
            hh = 2 * h

            # Calculate pressure differences
            delta_p_ij = p_j - p_i
            delta_p_ik = p_k - p_i
            delta_p_ji = p_i - p_j
            delta_p_jk = p_k - p_j
            delta_p_ki = p_i - p_k
            delta_p_kj = p_j - p_k

            # if is_it_boundary == False: berlin += 1
            # if is_it_boundary == False and berlin == 1:
            #     berlin += 1
            #     x = np.linspace(0.0001,0.01,100)
            #
            #     # Area
            #     A0 = [(calc_area_slice(x_K, x_0 + np.array([h, 0]), c_i, delta_p_ij) + calc_area_slice(x_0 + np.array([h, 0]), x_J, c_i, delta_p_ik) - calc_area_slice(x_K, x_0 - np.array([h, 0]), c_i, delta_p_ij) - calc_area_slice(x_0 - np.array([h, 0]), x_J, c_i, delta_p_ik)) / (2*h) for h in x]
            #     A1 = [(calc_area_slice(x_K, x_0 + np.array([0, h]), c_i, delta_p_ij) + calc_area_slice(x_0 + np.array([0, h]), x_J, c_i, delta_p_ik) - calc_area_slice(x_K, x_0 - np.array([0, h]), c_i, delta_p_ij) - calc_area_slice(x_0 - np.array([0, h]), x_J, c_i, delta_p_ik)) / (2*h) for h in x]
            #     A2 = [(calc_area_slice(x_K, x_0, c_i, delta_p_ij - h) + calc_area_slice(x_0, x_J, c_i, delta_p_ik - h) - calc_area_slice(x_K, x_0, c_i, delta_p_ij + h) - calc_area_slice(x_0, x_J, c_i, delta_p_ik + h)) / (2*h) for h in x]
            #     A3 = [(calc_area_slice(x_K, x_0, c_i, delta_p_ij + h) + calc_area_slice(x_0, x_J, c_i,  delta_p_ik) - calc_area_slice(x_K, x_0, c_i, delta_p_ij - h) - calc_area_slice(x_0, x_J, c_i, delta_p_ik)) / (2*h) for h in x]
            #     A4 = [(calc_area_slice(x_K, x_0, c_i, delta_p_ij) + calc_area_slice(x_0, x_J, c_i, delta_p_ik + h) - calc_area_slice(x_K, x_0, c_i, delta_p_ij) - calc_area_slice(x_0, x_J, c_i, delta_p_ik - h)) / (2*h) for h in x]
            #
            #     # Angle
            #     T0 = [(calc_angle(x_K, x_J, x_0 + np.array([h, 0]), delta_p_ij, delta_p_ik) - calc_angle(x_K, x_J, x_0 - np.array([h, 0]), delta_p_ij, delta_p_ik)) / (2*h) for h in x]
            #     T1 = [(calc_angle(x_K, x_J, x_0 + np.array([0, h]), delta_p_ij, delta_p_ik) - calc_angle(x_K, x_J, x_0 - np.array([0, h]), delta_p_ij, delta_p_ik)) / (2*h) for h in x]
            #     T2 = [(calc_angle(x_K, x_J, x_0, delta_p_ij - h, delta_p_ik - h) - calc_angle(x_K, x_J, x_0, delta_p_ij + h, delta_p_ik + h)) / (2*h) for h in x]
            #     T3 = [(calc_angle(x_K, x_J, x_0, delta_p_ij + h, delta_p_ik) - calc_angle(x_K, x_J, x_0, delta_p_ij - h, delta_p_ik)) / (2*h) for h in x]
            #     T4 = [(calc_angle(x_K, x_J, x_0, delta_p_ij, delta_p_ik + h) - calc_angle(x_K, x_J, x_0, delta_p_ij, delta_p_ik - h)) / (2*h) for h in x]
            #
            #
            #     plt.figure()
            #     # plt.plot(x, A0, label='dA_i/dx')
            #     # plt.plot(x, A1, label='dA_i/dy')
            #     # plt.plot(x, A2, label='dA_i/dpi')
            #     # plt.plot(x, A3, label='dA_i/dpj')
            #     # plt.plot(x, A4, label='dA_i/dpk')
            #     plt.plot(x, T0, label='d_theta_i/dx')
            #     plt.plot(x, T1, label='d_theta_i/dy')
            #     plt.plot(x, T2, label='d_theta_i/dpi')
            #     plt.plot(x, T3, label='d_theta_i/dpj')
            #     plt.plot(x, T4, label='d_theta_i/dpk')
            #     plt.legend(loc='upper right')
            #     plt.xlim([x[0], x[-1]])
            #     plt.xlabel('h')
            #     plt.ylabel('dA/dv')
            #     plt.show()


            ## CALCULATE ENTRIES IN JACOBIAN
            # First row in Jacobian
            theta_i_x  = (calc_angle(x_K, x_J, x_0 + h_x, delta_p_ij, delta_p_ik) - calc_angle(x_K, x_J, x_0 - h_x, delta_p_ij, delta_p_ik)) / hh
            theta_i_y  = (calc_angle(x_K, x_J, x_0 + h_y, delta_p_ij, delta_p_ik) - calc_angle(x_K, x_J, x_0 - h_y, delta_p_ij, delta_p_ik)) / hh
            theta_i_pi = (calc_angle(x_K, x_J, x_0, delta_p_ij - h, delta_p_ik - h) - calc_angle(x_K, x_J, x_0, delta_p_ij + h, delta_p_ik + h)) / hh
            theta_i_pj = (calc_angle(x_K, x_J, x_0, delta_p_ij + h, delta_p_ik) - calc_angle(x_K, x_J, x_0, delta_p_ij - h, delta_p_ik)) / hh
            theta_i_pk = (calc_angle(x_K, x_J, x_0, delta_p_ij, delta_p_ik + h) - calc_angle(x_K, x_J, x_0, delta_p_ij, delta_p_ik - h)) / hh

            # Second row in Jacobian
            theta_j_x  = (calc_angle(x_I, x_K, x_0 + h_x, delta_p_jk, delta_p_ji) - calc_angle(x_I, x_K, x_0 - h_x, delta_p_jk, delta_p_ji)) / hh
            theta_j_y  = (calc_angle(x_I, x_K, x_0 + h_y, delta_p_jk, delta_p_ji) - calc_angle(x_I, x_K, x_0 - h_y, delta_p_jk, delta_p_ji)) / hh
            theta_j_pi = (calc_angle(x_I, x_K, x_0, delta_p_jk, delta_p_ji + h) - calc_angle(x_I, x_K, x_0, delta_p_jk, delta_p_ji - h)) / hh
            theta_j_pj = (calc_angle(x_I, x_K, x_0, delta_p_jk - h, delta_p_ji) - calc_angle(x_I, x_K, x_0, delta_p_jk + h, delta_p_ji)) / hh
            theta_j_pk = (calc_angle(x_I, x_K, x_0, delta_p_jk + h, delta_p_ji) - calc_angle(x_I, x_K, x_0, delta_p_jk - h, delta_p_ji)) / hh

            # Thrid row in Jacobian
            A_i_x  = (calc_area_slice(x_K, x_0 + h_x, c_i, delta_p_ij) + calc_area_slice(x_0 + h_x, x_J, c_i, delta_p_ik) - calc_area_slice(x_K, x_0 - h_x, c_i, delta_p_ij) - calc_area_slice(x_0 - h_x, x_J, c_i, delta_p_ik)) / hh
            A_i_y  = (calc_area_slice(x_K, x_0 + h_y, c_i, delta_p_ij) + calc_area_slice(x_0 + h_y, x_J, c_i, delta_p_ik) - calc_area_slice(x_K, x_0 - h_y, c_i, delta_p_ij) - calc_area_slice(x_0 - h_y, x_J, c_i, delta_p_ik)) / hh
            A_i_pi = (calc_area_slice(x_K, x_0, c_i, delta_p_ij - h) + calc_area_slice(x_0, x_J, c_i, delta_p_ik - h) - calc_area_slice(x_K, x_0, c_i, delta_p_ij + h) - calc_area_slice(x_0, x_J, c_i, delta_p_ik + h)) / hh
            A_i_pj = (calc_area_slice(x_K, x_0, c_i, delta_p_ij + h) + calc_area_slice(x_0, x_J, c_i, delta_p_ik) - calc_area_slice(x_K, x_0, c_i, delta_p_ij - h) - calc_area_slice(x_0, x_J, c_i, delta_p_ik)) / hh
            A_i_pk = (calc_area_slice(x_K, x_0, c_i, delta_p_ij) + calc_area_slice(x_0, x_J, c_i, delta_p_ik + h) - calc_area_slice(x_K, x_0, c_i, delta_p_ij) - calc_area_slice(x_0, x_J, c_i, delta_p_ik - h)) / hh

            # Fourth row in Jacobian
            A_j_x  = (calc_area_slice(x_I, x_0 + h_x, c_j, delta_p_jk) + calc_area_slice(x_0 + h_x, x_K, c_j, delta_p_ji) - calc_area_slice(x_I, x_0 - h_x, c_j, delta_p_jk) - calc_area_slice(x_0 - h_x, x_K, c_j, delta_p_ji)) / hh
            A_j_y  = (calc_area_slice(x_I, x_0 + h_y, c_j, delta_p_jk) + calc_area_slice(x_0 + h_y, x_K, c_j, delta_p_ji) - calc_area_slice(x_I, x_0 - h_y, c_j, delta_p_jk) - calc_area_slice(x_0 - h_y, x_K, c_j, delta_p_ji)) / hh
            A_j_pi = (calc_area_slice(x_I, x_0, c_j, delta_p_jk) + calc_area_slice(x_0, x_K, c_j, delta_p_ji + h) - calc_area_slice(x_I, x_0, c_j, delta_p_jk) - calc_area_slice(x_0, x_K, c_j, delta_p_ji - h)) / hh
            A_j_pj = (calc_area_slice(x_I, x_0, c_j, delta_p_jk - h) + calc_area_slice(x_0, x_K, c_j, delta_p_ji - h) - calc_area_slice(x_I, x_0, c_j, delta_p_jk + h) - calc_area_slice(x_0, x_K, c_j, delta_p_ji + h)) / hh
            A_j_pk = (calc_area_slice(x_I, x_0, c_j, delta_p_jk + h) + calc_area_slice(x_0, x_K, c_j, delta_p_ji) - calc_area_slice(x_I, x_0, c_j, delta_p_jk - h) - calc_area_slice(x_0, x_K, c_j, delta_p_ji)) / hh

            # Fifth row in Jacobian
            A_k_x  = (calc_area_slice(x_J, x_0 + h_x, c_k, delta_p_ki) + calc_area_slice(x_0 + h_x, x_I, c_k, delta_p_kj) - calc_area_slice(x_J, x_0 - h_x, c_k, delta_p_ki) - calc_area_slice(x_0 - h_x, x_I, c_k, delta_p_kj)) / hh
            A_k_y  = (calc_area_slice(x_J, x_0 + h_y, c_k, delta_p_ki) + calc_area_slice(x_0 + h_y, x_I, c_k, delta_p_kj) - calc_area_slice(x_J, x_0 - h_y, c_k, delta_p_ki) - calc_area_slice(x_0 - h_y, x_I, c_k, delta_p_kj)) / hh
            A_k_pi = (calc_area_slice(x_J, x_0, c_k, delta_p_ki + h) + calc_area_slice(x_0, x_I, c_k, delta_p_kj) - calc_area_slice(x_J, x_0, c_k, delta_p_ki - h) - calc_area_slice(x_0, x_I, c_k, delta_p_kj)) / hh
            A_k_pj = (calc_area_slice(x_J, x_0, c_k, delta_p_ki) + calc_area_slice(x_0, x_I, c_k, delta_p_kj + h) - calc_area_slice(x_J, x_0, c_k, delta_p_ki) - calc_area_slice(x_0, x_I, c_k, delta_p_kj - h)) / hh
            A_k_pk = (calc_area_slice(x_J, x_0, c_k, delta_p_ki - h) + calc_area_slice(x_0, x_I, c_k, delta_p_kj - h) - calc_area_slice(x_J, x_0, c_k, delta_p_ki + h) - calc_area_slice(x_0, x_I, c_k, delta_p_kj + h)) / hh

            ## ASSEMBLE JACOBIAN
            if is_it_boundary == False:
                jacobian = np.array([[theta_i_x, theta_i_y, theta_i_pi, theta_i_pj, theta_i_pk],
                                     [theta_j_x, theta_j_y, theta_j_pi, theta_j_pj, theta_j_pk],
                                     [A_i_x, A_i_y, A_i_pi, A_i_pj, A_i_pk],
                                     [A_j_x, A_j_y, A_j_pi, A_j_pj, A_j_pk],
                                     [A_k_x, A_k_y, A_k_pi, A_k_pj, A_k_pk]])
            else:
                jacobian = np.array([[theta_i_x, theta_i_y, theta_i_pi, theta_i_pj, theta_i_pk],
                                     [theta_j_x, theta_j_y, theta_j_pi, theta_j_pj, theta_j_pk],
                                     [A_i_x, A_i_y, A_i_pi, A_i_pj, A_i_pk],
                                     [A_j_x, A_j_y, A_j_pi, A_j_pj, A_j_pk]])

            ## SET PROPERTY
            mesh.set_property(jacobian_handle, fh, jacobian)


def prop_is_flipped(mesh):
    # HANDLES
    global is_flipped_handle                                                                                            # Define and make property handle global, so that it is accessible from other functions
    is_flipped_handle = FPropHandle()
    mesh.add_property(is_flipped_handle, "fprop_float")

    for fh in mesh.faces():
        mesh.set_property(is_flipped_handle , fh, False)



# ----------------------------------------------------------------------------------------------------------------------
# Helper functions to property addition functions
# ----------------------------------------------------------------------------------------------------------------------


def test_finite_difference(mesh, fh):
    """ Draw a h + dF/dh plot to test, how large h can be before influencing the derivative
        OpenMesh TriMesh mesh: dual foam mesh
        OpenMesh FaceHandle fh: junction to investigate
    """
    x = np.linspace(0.0001,0.01,100)

    # Area
    A0 = [(calc_area_slice(x_K, x_0 + np.array([h, 0]), c_i, delta_p_ij) + calc_area_slice(x_0 + np.array([h, 0]), x_J, c_i, delta_p_ik) - calc_area_slice(x_K, x_0 - np.array([h, 0]), c_i, delta_p_ij) - calc_area_slice(x_0 - np.array([h, 0]), x_J, c_i, delta_p_ik)) / (2*h) for h in x]
    A1 = [(calc_area_slice(x_K, x_0 + np.array([0, h]), c_i, delta_p_ij) + calc_area_slice(x_0 + np.array([0, h]), x_J, c_i, delta_p_ik) - calc_area_slice(x_K, x_0 - np.array([0, h]), c_i, delta_p_ij) - calc_area_slice(x_0 - np.array([0, h]), x_J, c_i, delta_p_ik)) / (2*h) for h in x]
    A2 = [(calc_area_slice(x_K, x_0, c_i, delta_p_ij - h) + calc_area_slice(x_0, x_J, c_i, delta_p_ik - h) - calc_area_slice(x_K, x_0, c_i, delta_p_ij + h) - calc_area_slice(x_0, x_J, c_i, delta_p_ik + h)) / (2*h) for h in x]
    A3 = [(calc_area_slice(x_K, x_0, c_i, delta_p_ij + h) + calc_area_slice(x_0, x_J, c_i,  delta_p_ik) - calc_area_slice(x_K, x_0, c_i, delta_p_ij - h) - calc_area_slice(x_0, x_J, c_i, delta_p_ik)) / (2*h) for h in x]
    A4 = [(calc_area_slice(x_K, x_0, c_i, delta_p_ij) + calc_area_slice(x_0, x_J, c_i, delta_p_ik + h) - calc_area_slice(x_K, x_0, c_i, delta_p_ij) - calc_area_slice(x_0, x_J, c_i, delta_p_ik - h)) / (2*h) for h in x]

    # Angle
    T0 = [(calc_angle(x_K, x_J, x_0 + np.array([h, 0]), delta_p_ij, delta_p_ik) - calc_angle(x_K, x_J, x_0 - np.array([h, 0]), delta_p_ij, delta_p_ik)) / (2*h) for h in x]
    T1 = [(calc_angle(x_K, x_J, x_0 + np.array([0, h]), delta_p_ij, delta_p_ik) - calc_angle(x_K, x_J, x_0 - np.array([0, h]), delta_p_ij, delta_p_ik)) / (2*h) for h in x]
    T2 = [(calc_angle(x_K, x_J, x_0, delta_p_ij - h, delta_p_ik - h) - calc_angle(x_K, x_J, x_0, delta_p_ij + h, delta_p_ik + h)) / (2*h) for h in x]
    T3 = [(calc_angle(x_K, x_J, x_0, delta_p_ij + h, delta_p_ik) - calc_angle(x_K, x_J, x_0, delta_p_ij - h, delta_p_ik)) / (2*h) for h in x]
    T4 = [(calc_angle(x_K, x_J, x_0, delta_p_ij, delta_p_ik + h) - calc_angle(x_K, x_J, x_0, delta_p_ij, delta_p_ik - h)) / (2*h) for h in x]

    plt.figure()
    plt.plot(x, A0, label='dA_i/dx')
    plt.plot(x, A1, label='dA_i/dy')
    plt.plot(x, A2, label='dA_i/dpi')
    plt.plot(x, A3, label='dA_i/dpj')
    plt.plot(x, A4, label='dA_i/dpk')
    plt.legend(loc='upper right')
    plt.xlim([x[0], x[-1]])
    plt.xlabel('h')
    plt.ylabel('dA/dv')
    plt.show()

    plt.figure()
    plt.plot(x, T0, label='d_theta_i/dx')
    plt.plot(x, T1, label='d_theta_i/dy')
    plt.plot(x, T2, label='d_theta_i/dpi')
    plt.plot(x, T3, label='d_theta_i/dpj')
    plt.plot(x, T4, label='d_theta_i/dpk')
    plt.xlim([x[0], x[-1]])
    plt.xlabel('h')
    plt.ylabel('d_theta/dv')
    plt.show()


def calc_angle(x_K, x_J, x_0, delta_pj, delta_pk, gamma=gamma_imp):
    """ Calculate an angle given three points and two pressure differences
        np.array x_K, x_J: 2-dimensional junction vectors opposite of bubble k and j, i.e. adjacent to bubble i
        double delta_pj, delta_pk: pressure differences between cell i with j and k, respectively
        return double theta: angle
    """

    x_ij = x_K - x_0
    x_ik = x_J - x_0

    delta_pj = -delta_pj
    delta_pk = -delta_pk

    d_j = np.sqrt(np.sum((x_ij) ** 2))
    d_k = np.sqrt(np.sum((x_ik) ** 2))
    R_j = 2 * gamma / abs(delta_pj) if delta_pj != 0 else np.inf
    R_k = 2 * gamma / abs(delta_pk) if delta_pk != 0 else np.inf

    # Treat cases, where the radius of curvature is smaller than the radius of actual radius of curvature of the arc.
    if 0.5 * d_j > R_j: R_j = 0.5 * d_j
    if 0.5 * d_k > R_k: R_k = 0.5 * d_k

    #if np.dot(x_ij, x_ik) / (d_j * d_k)>1: print repr(np.dot(x_ij, x_ik) / (d_j * d_k)), repr(d_j / (2 * R_j))

    # Handle case, when fraction is just above 1 due to numerical limits (i.e. 1.00000000000000002)
    work_around = False
    if abs(np.dot(x_ij, x_ik) / (d_j * d_k)) > 1:
        work_around = True
        repl        = np.sign(np.dot(x_ij, x_ik) / (d_j * d_k))                                                         # Replacement

    if np.linalg.det([x_ij, x_ik]) > 0:
        if not work_around:
            theta = np.arccos(np.dot(x_ij, x_ik) / (d_j * d_k)) + (np.sign(delta_pj) * np.arcsin(d_j / (2 * R_j)) + np.sign(delta_pk) * np.arcsin(d_k / (2 * R_k)))
        else:
            theta = np.arccos(repl) + (np.sign(delta_pj) * np.arcsin(d_j / (2 * R_j)) + np.sign(delta_pk) * np.arcsin(d_k / (2 * R_k)))
    elif np.linalg.det([x_ij, x_ik]) < 0:
        if not work_around:
            theta = (2 * np.pi - np.arccos(np.dot(x_ij, x_ik) / (d_j * d_k))) + (np.sign(delta_pj) * np.arcsin(d_j / (2 * R_j)) + np.sign(delta_pk) * np.arcsin(d_k / (2 * R_k)))
        else:
            theta = (2 * np.pi - np.arccos(repl)) + (np.sign(delta_pj) * np.arcsin(d_j / (2 * R_j)) + np.sign(delta_pk) * np.arcsin(d_k / (2 * R_k)))
    else:
        theta = (np.sign(delta_pj) * np.arcsin(d_j / (2 * R_j)) + np.sign(delta_pk) * np.arcsin(d_k / (2 * R_k)))

    return theta


def calc_plateau_angles(mesh, fh, delta_alpha = [], gamma = gamma_imp, yes=False):
    """ Calculate the three angles in a specific Plateau junction
        TriMesh mesh: an existing OpenMesh mesh
        TriMesh FaceHandle fh: Specific face handle, if only the angles of one face are desired
        NumPy array delta_alpha: If a test of a different set of parameters is desired, this value is a 5x1 NymPy. Otherwise, it is an empty array
    """

    # Initialise
    thetas     = np.empty([3,2])

    # Fetch variables for face
    x_0 = np.array([mesh.property(mid_point_handle, fh)[0], mesh.property(mid_point_handle, fh)[1]])
    idx_i, idx_j, idx_k, fh_K, fh_J, fh_I, vh_i, vh_j, vh_k = mesh.property(face_connectivity_handle, fh)

    # Collect vectors in NumPy-array
    x_K = np.array([mesh.property(mid_point_handle, fh_K)[0], mesh.property(mid_point_handle, fh_K)[1]])                # Junction oppposite of K
    x_J = np.array([mesh.property(mid_point_handle, fh_J)[0], mesh.property(mid_point_handle, fh_J)[1]])                # Junction oppposite of J
    x_I = np.array([mesh.property(mid_point_handle, fh_I)[0], mesh.property(mid_point_handle, fh_I)[1]])                # Junction oppposite of I

    # Pressure values
    p_i = mesh.property(pressure_handle, vh_i)
    p_j = mesh.property(pressure_handle, vh_j)
    p_k = mesh.property(pressure_handle, vh_k)

    # Change variables, if desired
    if len(delta_alpha) != 0:
        x_0 = x_0 + np.array([delta_alpha[0][0], delta_alpha[1][0]])
        p_i = p_i + delta_alpha[2]
        p_j = p_j + delta_alpha[3]
        p_k = p_k + delta_alpha[4]

    # Calculate pressure differences
    delta_p_ij = p_j - p_i
    delta_p_ik = p_k - p_i
    delta_p_ji = p_i - p_j
    delta_p_jk = p_k - p_j
    delta_p_ki = p_i - p_k
    delta_p_kj = p_j - p_k

    # Calculate angles
    theta_i = calc_angle(x_K, x_J, x_0, delta_p_ij, delta_p_ik)
    theta_j = calc_angle(x_I, x_K, x_0, delta_p_jk, delta_p_ji)
    theta_k = calc_angle(x_J, x_I, x_0, delta_p_ki, delta_p_kj)

    # print calc_angle(x_0jk, x_0ij, delta_p_jk, delta_p_ji), calc_angle(x_0ij, x_0jk, delta_p_ij, delta_p_jk)

    # Add to list
    thetas[0, :] = [theta_i, idx_i]
    thetas[1, :] = [theta_j, idx_j]
    thetas[2, :] = [theta_k, idx_k]

    return thetas


def calc_area_slice(x_J, x_K, x_i, delta_p, gamma=gamma_imp):
    """ Calculate area segment between three adjacent points and bubble centre, if not specified
        np.array x_j, x_k: 2-dimensional 2d points on each end of bubble film
        double delta_p: pressure difference between current bubble and bubble on other side of bubble film
        return double area: area of slice
    """

    d         = np.sqrt(np.sum((x_K - x_J) ** 2))
    area_tria = 0.5 * abs(np.cross(x_J - x_i, x_K - x_i))

    if delta_p != 0:
        r = 2.0 * gamma / abs(delta_p)                                                                                  # Curvature radius
        #print d,r
        if 0.5 * d > abs(r): r = 0.5 * d                                                                                # Treat case, where the radius of curvature is smaller than the radius of actual radius of curvature of the arc.
        theta     = np.arccos((2 * abs(r)**2 - d**2) / (2 * abs(r)**2))
        area_temp = 0.5 * r ** 2 * (theta - np.sin(theta))
        area_segm = area_temp if delta_p > 0 else -area_temp
    else:
        area_segm = 0

    area = area_tria - area_segm

    return area


def calc_areas_bubble(mesh, vh, idx_change = np.array([-1, -1, -1, -1]), delta_alpha = np.array([0,0,0,0,0]), calc_target = True, gamma = gamma_imp, kappa = kappa_imp, delta_t = delta_t_imp):
    """ Calculate area and target area of bubble, if specified
        OpenMesh TriMesh mesh: dual foam mesh
        OpenMesh vertex handle vh: vertex handle of bubble cell
        int fh_change_idx: index of face handle indicating junction with modified position (in case of calculating a finite difference)
        double xy_change: Change in junction position (in case of calculating a finite difference)
        int vh_change_idx: index of vertex handle indicating bubble cell with modified pressure (in case of calculating a finite difference)
        double p_change: Change in pressure (in case of calculating a finite difference)
        bool calc_target: specification of whether target area is to be calculated or not
        int gamma: surface tension
        NumPy array delta_alpha: vector of variables changes in order, [x, y, p_bubble, p_neighbour_0, p_neighbour_1], where (x, y) are the changes to the junction to neighbour 0 and 1
        NumPy array indices: If change in one junction and pressures - indices identifying junction in question and adjacent neighbour bubbles
        return np.array areas: list of area and target area, if specified
    """

    if mesh.is_boundary(vh): return 0

    # Unpack indices for changes in position and pressure
    fh_idx   = idx_change[0]
    vh_i_idx = idx_change[1]
    vh_j_idx = idx_change[2]
    vh_k_idx = idx_change[3]

    # Find coordinate of this vertex
    coor_bubble = np.array([mesh.point(vh)[0], mesh.point(vh)[1]])
    idx_bubble  = vh.idx()
    if idx_bubble == vh_i_idx: p_bubble = mesh.property(pressure_handle, vh) + delta_alpha[2][0]
    elif idx_bubble == vh_j_idx: p_bubble = mesh.property(pressure_handle, vh) + delta_alpha[3][0]
    elif idx_bubble == vh_k_idx: p_bubble = mesh.property(pressure_handle, vh) + delta_alpha[4][0]
    else: p_bubble = mesh.property(pressure_handle, vh)

    # Initialise area counter
    area = 0
    if calc_target: area_target_temp = 0

    # Build faces list
    faces = [fh for fh in mesh.vf(vh) if mesh.is_boundary(fh) == False]
    faces.append(faces[0])

    for i in range(len(faces) - 1):
        # Find vertex belonging to two last faces and to adjacent bubble
        vertices_0 = np.array([[vh0.idx(), vh0] for vh0 in mesh.fv(faces[i]) if vh0.idx() != idx_bubble])
        vertices_1 = np.array([[vh0.idx(), vh0] for vh0 in mesh.fv(faces[i + 1]) if vh0.idx() != idx_bubble])
        idx        = np.in1d(vertices_0[:, 0], vertices_1[:, 0])

        # Find common vertex != vh_bubble_idx, i.e. vertex of neighbour cell in question
        if sum(idx) == 0:                                                                                               # Handle case of boundary film, where two or more world bubble vertices are present
            common_vh = [vh0 for vh0 in mesh.fv(faces[i]) if mesh.is_boundary(vh0)][0]
        else:
            common_vh = vertices_0[idx, 1][0]

        coor_0 = np.array([mesh.property(mid_point_handle, faces[i])[0], mesh.property(mid_point_handle, faces[i])[1]])
        coor_1 = np.array([mesh.property(mid_point_handle, faces[i + 1])[0], mesh.property(mid_point_handle, faces[i + 1])[1]])

        # Update junction coordinates, if calculating finite difference
        delta_alpha_xy = np.array([np.double(delta_alpha[0]), np.double(delta_alpha[1])])
        if faces[i].idx() == fh_idx: coor_0 += delta_alpha_xy
        elif faces[i + 1].idx() == fh_idx: coor_1 += delta_alpha_xy

        # Update pressure, if calculating finite difference
        if common_vh.idx() == vh_i_idx: p_neighbour = mesh.property(pressure_handle, common_vh) + delta_alpha[2][0]
        elif common_vh.idx() == vh_j_idx: p_neighbour = mesh.property(pressure_handle, common_vh) + delta_alpha[3][0]
        elif common_vh.idx() == vh_k_idx: p_neighbour = mesh.property(pressure_handle, common_vh) + delta_alpha[4][0]
        else: p_neighbour = mesh.property(pressure_handle, common_vh)

        delta_p = p_neighbour - p_bubble

        # Calculate area of bubble
        area += calc_area_slice(coor_0, coor_1, coor_bubble, delta_p, gamma)

        # Calculate target area of bubble
        if calc_target: area_target_temp += delta_p * calc_film_length(coor_0, coor_1, delta_p, gamma)


    if calc_target: area_target = area + delta_t * kappa * area_target_temp

    #if calc_target: print area, area_target, area_target_temp, delta_p, len(faces)

    # if vh.idx()==105 and calc_target: print p_bubble, area, area_target, area_target-area
    # if vh.idx() == 105: print p_bubble

    if calc_target: return np.array([area, area_target])
    else: return area


def calc_film_length(x_j, x_k, delta_p, gamma=gamma_imp):
    """ Calculate film_length
        np.array x_j, x_k: 2-dimensional 2d points on each end of bubble film
        double delta_p: pressure difference between current bubble and bubble on other side of bubble film
        return double area: area of slice
    """

    d       = np.sqrt(np.sum((x_k - x_j) ** 2))                                                                         # Calculate distance between two end points of bubble film

    if delta_p != 0:
        r = 2.0 * gamma / abs(delta_p)  # Curvature radius
        if 0.5 * d > r: r = 0.5 * d                                                                                     # Treat case, where the radius of curvature is smaller than the radius of actual radius of curvature of the arc.
        theta = np.arccos((2 * abs(r)**2 - d**2) / (2 * abs(r)**2))
        film_length = theta * r
    else:
        film_length = d

    return film_length


def calc_jacobian(mesh, fh, h = h_imp, draw_dh = draw_dh_imp, idx_draw = idx_draw_imp, delta_alpha = []):
    is_it_boundary = sum([mesh.is_boundary(vh0) for vh0 in mesh.fv(fh)])

    # Fetch attribute values
    x_0 = np.array([mesh.property(mid_point_handle, fh)[0], mesh.property(mid_point_handle, fh)[1]])  # Coordinate of junction
    idx_i, idx_j, idx_k, fh_K, fh_J, fh_I, vh_i, vh_j, vh_k = mesh.property(face_connectivity_handle, fh)

    # Collect vectors in NumPy-array
    x_K = np.array([mesh.property(mid_point_handle, fh_K)[0], mesh.property(mid_point_handle, fh_K)[1]])  # Junction oppposite of K
    x_J = np.array([mesh.property(mid_point_handle, fh_J)[0], mesh.property(mid_point_handle, fh_J)[1]])  # Junction oppposite of J
    x_I = np.array([mesh.property(mid_point_handle, fh_I)[0], mesh.property(mid_point_handle, fh_I)[1]])  # Junction oppposite of I

    # Coordinates of bubble centres
    c_i = np.array([mesh.point(vh_i)[0], mesh.point(vh_i)[1]])
    c_j = np.array([mesh.point(vh_j)[0], mesh.point(vh_j)[1]])
    c_k = np.array([mesh.point(vh_k)[0], mesh.point(vh_k)[1]])

    # Pressure values
    p_i = mesh.property(pressure_handle, vh_i)
    p_j = mesh.property(pressure_handle, vh_j)
    p_k = mesh.property(pressure_handle, vh_k)

    # Change variables, if desired
    if len(delta_alpha) != 0:
        x_0 = x_0 + np.array([delta_alpha[0][0], delta_alpha[1][0]])
        p_i = p_i + delta_alpha[2][0]
        p_j = p_j + delta_alpha[3][0]
        p_k = p_k + delta_alpha[4][0]

    # Create finite difference vectors
    h_x = np.array([h, 0])
    h_y = np.array([0, h])

    # Central difference
    hh = 2 * h

    # Calculate pressure differences
    delta_p_ij = p_j - p_i
    delta_p_ik = p_k - p_i
    delta_p_ji = p_i - p_j
    delta_p_jk = p_k - p_j
    delta_p_ki = p_i - p_k
    delta_p_kj = p_j - p_k

    ## CALCULATE ENTRIES IN JACOBIAN
    # First row in Jacobian
    theta_i_x  = (calc_angle(x_K, x_J, x_0 + h_x, delta_p_ij, delta_p_ik) - calc_angle(x_K, x_J, x_0 - h_x, delta_p_ij, delta_p_ik)) / hh
    theta_i_y  = (calc_angle(x_K, x_J, x_0 + h_y, delta_p_ij, delta_p_ik) - calc_angle(x_K, x_J, x_0 - h_y, delta_p_ij, delta_p_ik)) / hh
    theta_i_pi = (calc_angle(x_K, x_J, x_0, delta_p_ij - h, delta_p_ik - h) - calc_angle(x_K, x_J, x_0, delta_p_ij + h, delta_p_ik + h)) / hh
    theta_i_pj = (calc_angle(x_K, x_J, x_0, delta_p_ij + h, delta_p_ik) - calc_angle(x_K, x_J, x_0, delta_p_ij - h, delta_p_ik)) / hh
    theta_i_pk = (calc_angle(x_K, x_J, x_0, delta_p_ij, delta_p_ik + h) - calc_angle(x_K, x_J, x_0, delta_p_ij, delta_p_ik - h)) / hh


    # Second row in Jacobian
    theta_j_x  = (calc_angle(x_I, x_K, x_0 + h_x, delta_p_jk, delta_p_ji) - calc_angle(x_I, x_K, x_0 - h_x, delta_p_jk, delta_p_ji)) / hh
    theta_j_y  = (calc_angle(x_I, x_K, x_0 + h_y, delta_p_jk, delta_p_ji) - calc_angle(x_I, x_K, x_0 - h_y, delta_p_jk, delta_p_ji)) / hh
    theta_j_pi = (calc_angle(x_I, x_K, x_0, delta_p_jk, delta_p_ji + h) - calc_angle(x_I, x_K, x_0, delta_p_jk, delta_p_ji - h)) / hh
    theta_j_pj = (calc_angle(x_I, x_K, x_0, delta_p_jk - h, delta_p_ji - h) - calc_angle(x_I, x_K, x_0, delta_p_jk + h, delta_p_ji + h)) / hh
    theta_j_pk = (calc_angle(x_I, x_K, x_0, delta_p_jk + h, delta_p_ji) - calc_angle(x_I, x_K, x_0, delta_p_jk - h, delta_p_ji)) / hh


    # Thrid row in Jacobian
    A_i_x  = (calc_area_slice(x_K, x_0 + h_x, c_i, delta_p_ij) + calc_area_slice(x_0 + h_x, x_J, c_i, delta_p_ik) - calc_area_slice(x_K, x_0 - h_x, c_i, delta_p_ij) - calc_area_slice(x_0 - h_x, x_J, c_i, delta_p_ik)) / hh
    A_i_y  = (calc_area_slice(x_K, x_0 + h_y, c_i, delta_p_ij) + calc_area_slice(x_0 + h_y, x_J, c_i, delta_p_ik) - calc_area_slice(x_K, x_0 - h_y, c_i, delta_p_ij) - calc_area_slice(x_0 - h_y, x_J, c_i, delta_p_ik)) / hh
    #A_i_pi = (calc_area_slice(x_K, x_0, c_i, delta_p_ij - h) + calc_area_slice(x_0, x_J, c_i, delta_p_ik - h) - calc_area_slice(x_K, x_0, c_i, delta_p_ij + h) - calc_area_slice(x_0, x_J, c_i, delta_p_ik + h)) / hh
    A_i_pi = (calc_areas_bubble(mesh, vh_i, np.array([-1, idx_i, -1, -1]), np.array([[0],[0],[h],[0],[0]]), calc_target=False) - calc_areas_bubble(mesh, vh_i, np.array([-1, idx_i, -1, -1]), np.array([[0],[0],[-h],[0],[0]]), calc_target=False)) / hh
    A_i_pj = (calc_area_slice(x_K, x_0, c_i, delta_p_ij + h) - calc_area_slice(x_K, x_0, c_i, delta_p_ij - h)) / hh
    A_i_pk = (calc_area_slice(x_0, x_J, c_i, delta_p_ik + h) - calc_area_slice(x_0, x_J, c_i, delta_p_ik - h)) / hh

    # Fourth row in Jacobian
    A_j_x  = (calc_area_slice(x_I, x_0 + h_x, c_j, delta_p_jk) + calc_area_slice(x_0 + h_x, x_K, c_j, delta_p_ji) - calc_area_slice(x_I, x_0 - h_x, c_j, delta_p_jk) - calc_area_slice(x_0 - h_x, x_K, c_j, delta_p_ji)) / hh
    A_j_y  = (calc_area_slice(x_I, x_0 + h_y, c_j, delta_p_jk) + calc_area_slice(x_0 + h_y, x_K, c_j, delta_p_ji) - calc_area_slice(x_I, x_0 - h_y, c_j, delta_p_jk) - calc_area_slice(x_0 - h_y, x_K, c_j, delta_p_ji)) / hh
    A_j_pi = (calc_area_slice(x_0, x_K, c_j, delta_p_ji + h) - calc_area_slice(x_0, x_K, c_j, delta_p_ji - h)) / hh
    #A_j_pj = (calc_area_slice(x_I, x_0, c_j, delta_p_jk - h) + calc_area_slice(x_0, x_K, c_j, delta_p_ji - h) - calc_area_slice(x_I, x_0, c_j, delta_p_jk + h) - calc_area_slice(x_0, x_K, c_j, delta_p_ji + h)) / hh
    A_j_pj = (calc_areas_bubble(mesh, vh_j, np.array([-1, -1, idx_j, -1]), np.array([[0],[0],[0],[h],[0]]), calc_target=False) - calc_areas_bubble(mesh, vh_j, np.array([-1, -1, vh_j, -1]), np.array([[0],[0],[0],[-h],[0]]), calc_target=False)) / hh
    A_j_pk = (calc_area_slice(x_I, x_0, c_j, delta_p_jk + h) - calc_area_slice(x_I, x_0, c_j, delta_p_jk - h)) / hh

    # Fifth row in Jacobian
    A_k_x  = (calc_area_slice(x_J, x_0 + h_x, c_k, delta_p_ki) + calc_area_slice(x_0 + h_x, x_I, c_k, delta_p_kj) - calc_area_slice(x_J, x_0 - h_x, c_k, delta_p_ki) - calc_area_slice(x_0 - h_x, x_I, c_k, delta_p_kj)) / hh
    A_k_y  = (calc_area_slice(x_J, x_0 + h_y, c_k, delta_p_ki) + calc_area_slice(x_0 + h_y, x_I, c_k, delta_p_kj) - calc_area_slice(x_J, x_0 - h_y, c_k, delta_p_ki) - calc_area_slice(x_0 - h_y, x_I, c_k, delta_p_kj)) / hh
    A_k_pi = (calc_area_slice(x_J, x_0, c_k, delta_p_ki + h) - calc_area_slice(x_J, x_0, c_k, delta_p_ki - h)) / hh
    A_k_pj = (calc_area_slice(x_0, x_I, c_k, delta_p_kj + h) - calc_area_slice(x_0, x_I, c_k, delta_p_kj - h)) / hh
    #A_k_pk = (calc_area_slice(x_J, x_0, c_k, delta_p_ki - h) + calc_area_slice(x_0, x_I, c_k, delta_p_kj - h) - calc_area_slice(x_J, x_0, c_k, delta_p_ki + h) - calc_area_slice(x_0, x_I, c_k, delta_p_kj + h)) / hh
    A_k_pk = (calc_areas_bubble(mesh, vh_k, np.array([-1, -1, -1, idx_k]), np.array([[0],[0],[0],[0],[h]]), calc_target=False) - calc_areas_bubble(mesh, vh_k, np.array([-1, -1, -1, idx_k]), np.array([[0],[0],[0],[0],[-h]]), calc_target=False)) / hh

    ## ASSEMBLE JACOBIAN
    if is_it_boundary == False:
        jacobian = np.array([[theta_i_x, theta_i_y, theta_i_pi, theta_i_pj, theta_i_pk],
                             [theta_j_x, theta_j_y, theta_j_pi, theta_j_pj, theta_j_pk],
                             [A_i_x, A_i_y, A_i_pi, A_i_pj, A_i_pk],
                             [A_j_x, A_j_y, A_j_pi, A_j_pj, A_j_pk],
                             [A_k_x, A_k_y, A_k_pi, A_k_pj, A_k_pk]])
    else:
        jacobian = np.array([[theta_i_x, theta_i_y, theta_i_pi, theta_i_pj],
                             [theta_j_x, theta_j_y, theta_j_pi, theta_j_pj],
                             [A_i_x, A_i_y, A_i_pi, A_i_pj],
                             [A_j_x, A_j_y, A_j_pi, A_j_pj]])


    # Draw 2 plots of the finite differences of an A and a theta wrt h, if desired
    if draw_dh and fh.idx() == idx_draw: test_finite_difference(mesh, fh)

    return jacobian


def calc_global_model(mesh, h = h_imp):
    # Count number of junctions
    number_junctions = 0
    for fh0 in mesh.faces():
        if not mesh.is_boundary(fh0):
            number_junctions += 1

    # Count number of bubble cell
    number_cells = 0
    for vh0 in mesh.vertices():
        if not mesh.is_boundary(vh0):
            number_cells += 1

    # Initialise Jacobian and minus_F
    rows = number_junctions * 2 + number_cells
    cols = rows
    jacobian = np.zeros([rows, cols])
    minus_F  = np.zeros([cols, 1])

    # Initialise iteration counter
    face_counter = 0

    # Create mapping lists
    vertex_mapping = np.array([vh0.idx() for vh0 in mesh.vertices() if mesh.is_boundary(vh0) == False])
    face_mapping   = np.array([fh0.idx() for fh0 in mesh.faces() if mesh.is_boundary(fh0) == False])

    # Build angle part of Jacobian
    for fh in mesh.faces():
        if not mesh.is_boundary(fh):
            # Fetch attributes
            idx_i, idx_j, idx_k, fh_K, fh_J, fh_I, vh_i, vh_j, vh_k = mesh.property(face_connectivity_handle, fh)
            fh_idx = fh.idx()

            # Collect vectors in NumPy-array
            x_0 = np.array([mesh.property(mid_point_handle, fh)[0], mesh.property(mid_point_handle, fh)[1]])
            x_K = np.array([mesh.property(mid_point_handle, fh_K)[0], mesh.property(mid_point_handle, fh_K)[1]])        # Junction oppposite of K
            x_J = np.array([mesh.property(mid_point_handle, fh_J)[0], mesh.property(mid_point_handle, fh_J)[1]])        # Junction oppposite of J
            x_I = np.array([mesh.property(mid_point_handle, fh_I)[0], mesh.property(mid_point_handle, fh_I)[1]])        # Junction oppposite of I

            # Pressure values
            p_i = mesh.property(pressure_handle, vh_i)
            p_j = mesh.property(pressure_handle, vh_j)
            p_k = mesh.property(pressure_handle, vh_k)

            # Calculate pressure differences
            delta_p_ij = p_j - p_i
            delta_p_ik = p_k - p_i
            delta_p_ji = p_i - p_j
            delta_p_jk = p_k - p_j
            delta_p_ki = p_i - p_k
            delta_p_kj = p_j - p_k

            # Central difference
            hh = 2 * h

            # Create finite difference vectors
            h_x = np.array([h, 0])
            h_y = np.array([0, h])

            # Calculate Jacobian entries
            theta_i_x0 = (calc_angle(x_K, x_J, x_0 + h_x, delta_p_ij, delta_p_ik) - calc_angle(x_K, x_J, x_0 - h_x, delta_p_ij, delta_p_ik)) / hh
            theta_i_y0 = (calc_angle(x_K, x_J, x_0 + h_y, delta_p_ij, delta_p_ik) - calc_angle(x_K, x_J, x_0 - h_y, delta_p_ij, delta_p_ik)) / hh
            theta_i_xK = (calc_angle(x_K + h_x, x_J, x_0, delta_p_ij, delta_p_ik) - calc_angle(x_K - h_x, x_J, x_0, delta_p_ij, delta_p_ik)) / hh
            theta_i_yK = (calc_angle(x_K + h_y, x_J, x_0, delta_p_ij, delta_p_ik) - calc_angle(x_K - h_y, x_J, x_0, delta_p_ij, delta_p_ik)) / hh
            theta_i_xJ = (calc_angle(x_K, x_J + h_x, x_0, delta_p_ij, delta_p_ik) - calc_angle(x_K, x_J - h_x, x_0, delta_p_ij, delta_p_ik)) / hh
            theta_i_yJ = (calc_angle(x_K, x_J + h_y, x_0, delta_p_ij, delta_p_ik) - calc_angle(x_K, x_J - h_y, x_0, delta_p_ij, delta_p_ik)) / hh
            theta_i_pi = (calc_angle(x_K, x_J, x_0, delta_p_ij - h, delta_p_ik - h) - calc_angle(x_K, x_J, x_0, delta_p_ij + h, delta_p_ik + h)) / hh
            theta_i_pj = (calc_angle(x_K, x_J, x_0, delta_p_ij + h, delta_p_ik) - calc_angle(x_K, x_J, x_0, delta_p_ij - h, delta_p_ik)) / hh
            if not mesh.is_boundary(vh_k): theta_i_pk = (calc_angle(x_K, x_J, x_0, delta_p_ij, delta_p_ik + h) - calc_angle(x_K, x_J, x_0, delta_p_ij, delta_p_ik - h)) / hh
            theta_j_x0 = (calc_angle(x_I, x_K, x_0 + h_x, delta_p_jk, delta_p_ji) - calc_angle(x_I, x_K, x_0 - h_x, delta_p_jk, delta_p_ji)) / hh
            theta_j_y0 = (calc_angle(x_I, x_K, x_0 + h_y, delta_p_jk, delta_p_ji) - calc_angle(x_I, x_K, x_0 - h_y, delta_p_jk, delta_p_ji)) / hh
            theta_j_xK = (calc_angle(x_I, x_K + h_x, x_0, delta_p_jk, delta_p_ji) - calc_angle(x_I, x_K - h_x, x_0, delta_p_jk, delta_p_ji)) / hh
            theta_j_yK = (calc_angle(x_I, x_K + h_y, x_0, delta_p_jk, delta_p_ji) - calc_angle(x_I, x_K - h_y, x_0, delta_p_jk, delta_p_ji)) / hh
            theta_j_xI = (calc_angle(x_I + h_x, x_K, x_0, delta_p_jk, delta_p_ji) - calc_angle(x_I - h_x, x_K, x_0, delta_p_jk, delta_p_ji)) / hh
            theta_j_yI = (calc_angle(x_I + h_y, x_K, x_0, delta_p_jk, delta_p_ji) - calc_angle(x_I - h_y, x_K, x_0, delta_p_jk, delta_p_ji)) / hh
            theta_j_pi = (calc_angle(x_I, x_K, x_0, delta_p_jk, delta_p_ji + h) - calc_angle(x_I, x_K, x_0, delta_p_jk, delta_p_ji - h)) / hh
            theta_j_pj = (calc_angle(x_I, x_K, x_0, delta_p_jk - h, delta_p_ji - h) - calc_angle(x_I, x_K, x_0, delta_p_jk + h, delta_p_ji + h)) / hh
            if not mesh.is_boundary(vh_k): theta_j_pk = (calc_angle(x_I, x_K, x_0, delta_p_jk + h, delta_p_ji) - calc_angle(x_I, x_K, x_0, delta_p_jk - h, delta_p_ji)) / hh

            #print theta_j_yI
            #print theta_i_y0


            # Calculate minus_F entries
            theta_delta_i = 2.0 * np.pi / 3.0 - calc_angle(x_K, x_J, x_0, delta_p_ij, delta_p_ik)
            theta_delta_j = 2.0 * np.pi / 3.0 - calc_angle(x_I, x_K, x_0, delta_p_jk, delta_p_ji)

            # Find indices for angles
            wrt_x0 = np.where(face_mapping == fh_idx)[0][0]                                                             # Column numbers
            wrt_y0 = number_junctions + wrt_x0
            wrt_xK = np.where(face_mapping == fh_K.idx())[0][0]
            wrt_yK = number_junctions + wrt_xK
            wrt_xJ = np.where(face_mapping == fh_J.idx())[0][0]
            wrt_yJ = number_junctions + wrt_xJ
            wrt_xI = np.where(face_mapping == fh_I.idx())[0][0]
            wrt_yI = number_junctions + wrt_xI
            wrt_pi = number_junctions * 2 + np.where(vertex_mapping == idx_i)[0][0]
            wrt_pj = number_junctions * 2 + np.where(vertex_mapping == idx_j)[0][0]
            if not mesh.is_boundary(vh_k): wrt_pk = number_junctions * 2 + np.where(vertex_mapping == idx_k)[0][0]
            cell_i = np.where(face_mapping == fh_idx)[0][0] * 2                                                                                   # Row numbers
            cell_j = cell_i + 1

            # Update entry values in Jacobian for these three cells
            jacobian[cell_i, wrt_x0] = theta_i_x0
            jacobian[cell_i, wrt_y0] = theta_i_y0

            jacobian[cell_i, wrt_xK] = theta_i_xK
            jacobian[cell_i, wrt_yK] = theta_i_yK
            jacobian[cell_i, wrt_xJ] = theta_i_xJ
            jacobian[cell_i, wrt_yJ] = theta_i_yJ

            jacobian[cell_i, wrt_pi] = theta_i_pi
            jacobian[cell_i, wrt_pj] = theta_i_pj
            if not mesh.is_boundary(vh_k): jacobian[cell_i, wrt_pk] = theta_i_pk
            jacobian[cell_j, wrt_x0] = theta_j_x0
            jacobian[cell_j, wrt_y0] = theta_j_y0

            jacobian[cell_j, wrt_xK] = theta_j_xK
            jacobian[cell_j, wrt_yK] = theta_j_yK
            jacobian[cell_j, wrt_xI] = theta_j_xI
            jacobian[cell_j, wrt_yI] = theta_j_yI

            jacobian[cell_j, wrt_pi] = theta_j_pi
            jacobian[cell_j, wrt_pj] = theta_j_pj
            if not mesh.is_boundary(vh_k): jacobian[cell_j, wrt_pk] = theta_j_pk

            # Update entries in -F(alpha)
            minus_F[cell_i] = theta_delta_i
            minus_F[cell_j] = theta_delta_j

            is_it_boundary = sum([mesh.is_boundary(vh0) for vh0 in mesh.fv(fh)])

            # Increment counter
            face_counter += 1

    # Initialise iteration counter
    vertex_counter = face_counter * 2

    # Build area part of Jacobian
    for vh in mesh.vertices():
        if mesh.is_boundary(vh) == False:
            # Fetch index for bubble
            cell = np.where(vertex_mapping == vh.idx())[0][0] + number_junctions * 2

            is_it_boundary = sum([mesh.is_boundary(vh0) for vh0 in mesh.vv(vh)])

            # Handle dAdx and dAdy-cases
            for fh0 in mesh.vf(vh):
                if mesh.is_boundary(fh0) == False:
                    A_plus_h  = calc_areas_bubble(mesh, vh, idx_change=np.array([fh0.idx(), -1, -1, -1]), delta_alpha=np.array([[h], [0], [0], [0], [0]]), calc_target=False)
                    A_minus_h = calc_areas_bubble(mesh, vh, idx_change=np.array([fh0.idx(), -1, -1, -1]), delta_alpha=np.array([[-h], [0], [0], [0], [0]]), calc_target=False)
                    dAdx = (A_plus_h - A_minus_h) / hh

                    A_plus_h  = calc_areas_bubble(mesh, vh, idx_change=np.array([fh0.idx(), -1, -1, -1]), delta_alpha=np.array([[0], [h], [0], [0], [0]]), calc_target=False)
                    A_minus_h = calc_areas_bubble(mesh, vh, idx_change=np.array([fh0.idx(), -1, -1, -1]), delta_alpha=np.array([[0], [-h], [0], [0], [0]]), calc_target=False)
                    dAdy = (A_plus_h - A_minus_h) / hh

                    # Find indices
                    wrt_x = np.where(face_mapping == fh0.idx())[0][0]
                    wrt_y = number_junctions + wrt_x

                    # Assign value in Jacobian
                    jacobian[cell, wrt_x] = dAdx
                    jacobian[cell, wrt_y] = dAdy

            # Handle dAdp-cases
            for vh0 in mesh.vv(vh):
                if mesh.is_boundary(vh0) == False:
                    # Calculate finite difference
                    A_plus_h  = calc_areas_bubble(mesh, vh, idx_change=np.array([-1, vh0.idx(), -1, -1]), delta_alpha=np.array([[0], [0], [h], [0], [0]]), calc_target=False)
                    A_minus_h = calc_areas_bubble(mesh, vh, idx_change=np.array([-1, vh0.idx(), -1, -1]), delta_alpha=np.array([[0], [0], [-h], [0], [0]]), calc_target=False)
                    dAdp = (A_plus_h - A_minus_h) / hh

                    # Find indices
                    wrt_p = number_junctions * 2 + np.where(vertex_mapping == vh0.idx())[0][0]

                    # Assign value in Jacobian
                    jacobian[cell, wrt_p] = dAdp

            # Handle dA0/dp0
            A_plus_h  = calc_areas_bubble(mesh, vh, idx_change=np.array([-1, vh.idx(), -1, -1]), delta_alpha=np.array([[0], [0], [h], [0], [0]]), calc_target=False)
            A_minus_h = calc_areas_bubble(mesh, vh, idx_change=np.array([-1, vh.idx(), -1, -1]), delta_alpha=np.array([[0], [0], [-h], [0], [0]]), calc_target=False)
            dAdp      = (A_plus_h - A_minus_h) / hh
            wrt_p     = number_junctions * 2 + np.where(vertex_mapping == vh.idx())[0][0]
            jacobian[cell, wrt_p] = dAdp

            # Calculate and update entries in -F(alpha)
            areas         = calc_areas_bubble(mesh, vh)
            delta_areas   = areas[1] - areas[0]
            minus_F[cell] = delta_areas

            # Increment counter
            vertex_counter += 1


    return jacobian, minus_F, vertex_mapping, face_mapping



# ----------------------------------------------------------------------------------------------------------------------
# Save and load mesh along with properties
# ----------------------------------------------------------------------------------------------------------------------

def save_foam_mesh(mesh, iteration, method_name):
    """ Save mesh at any time in iteration along with its attributes. Load again with load_foam_mesh
        int iteration: iteration number, for how far the simulation has come
        string method_name: name of method ('local', 'quasi_global' or 'global')
    """

    # Save positions in array
    fh_idxs = []
    xs = []
    ys = []
    for fh in mesh.faces():
        if not mesh.is_boundary(fh):
            # Fetch coordinates
            xy = mesh.property(mid_point_handle, fh)
            x = xy[0]
            y = xy[1]

            # Store coordinates in array
            fh_idxs.append(fh.idx())
            xs.append(x)
            ys.append(y)

    # Save pressures in array
    vh_idxs = []
    ps = []
    for vh in mesh.vertices():
        # Fetch pressure value
        p = mesh.property(pressure_handle, vh)

        # Store pressure in array
        vh_idxs.append(vh.idx())
        ps.append(p)

    # Collect in NumPy-array
    coordinates = np.array([fh_idxs, xs, ys])
    pressures = np.array([vh_idxs, ps])

    # Save mesh and attribute vectors
    np.save('saved/saved_mesh_' + method_name + '_' + str(iteration) + '_iterations_xy', coordinates)
    np.save('saved/saved_mesh_' + method_name + '_' + str(iteration) + '_iterations_p', pressures)
    write_mesh(mesh, 'saved/saved_mesh_' + method_name + '_' + str(iteration) + '.obj')


def load_foam_mesh(method_name, iteration):
    """ Load mesh and assign pressures and
        int iteration: iteration number, for how far the simulation has come
        string method_name: name of method ('local', 'quasi_global' or 'global')
    """

    # Load mesh (vertices and faces)
    mesh = TriMesh()
    read_mesh(mesh, 'saved/initial_mesh/saved_mesh_' + method_name + '_' + str(iteration) + '.obj')
    #read_mesh(mesh, 'saved/initial_mesh_same_pressure/saved_mesh_' + method_name + '_' + str(iteration) + '.obj')
    #read_mesh(mesh, 'saved/large_mesh/saved_mesh_' + method_name + '_' + str(iteration) + '.obj')

    # Initialise properties
    global pressure_handle
    global mid_point_handle
    pressure_handle = VPropHandle()
    mid_point_handle = FPropHandle()
    mesh.add_property(pressure_handle, "vprop_float")
    mesh.add_property(mid_point_handle, "fprop_float")

    # Add junction coordinates
    coordinates = np.load('saved/initial_mesh/saved_mesh_' + method_name + '_' + str(iteration) + '_iterations_xy.npy')
    #coordinates = np.load('saved/initial_mesh_same_pressure/saved_mesh_' + method_name + '_' + str(iteration) + '_iterations_xy.npy')
    #coordinates = np.load('saved/large_mesh/saved_mesh_' + method_name + '_' + str(iteration) + '_iterations_xy.npy')
    for fh in mesh.faces():
        if not mesh.is_boundary(fh):
            idx  = np.where(coordinates[0, :] == fh.idx())[0][0]
            coor = TriMesh.Point(coordinates[1,idx], coordinates[2,idx], 0)
            mesh.set_property(mid_point_handle, fh, coor)

    # Save pressures in array
    pressures = np.load('saved/initial_mesh/saved_mesh_' + method_name + '_' + str(iteration) + '_iterations_p.npy')
    #pressures = np.load('saved/initial_mesh_same_pressure/saved_mesh_' + method_name + '_' + str(iteration) + '_iterations_p.npy')
    #pressures = np.load('saved/large_mesh/saved_mesh_' + method_name + '_' + str(iteration) + '_iterations_p.npy')
    for vh in mesh.vertices():
        # Fetch pressure value
        idx = np.where(pressures[0, :] == vh.idx())[0][0]
        p   = pressures[1,idx]
        mesh.set_property(pressure_handle, vh, p)

    add_face_connectivity(mesh)

    return mesh




