from foam_import_openmesh import *
import foam_properties as fp
import foam_visualise as fv
import foam_topology as ft
import matplotlib.pyplot as plt
import matplotlib.path as mplPath
import config


# ----------------------------------------------------------------------------------------------------------------------
# Import constants
# ----------------------------------------------------------------------------------------------------------------------
beta_imp       = config.backtrack['beta']
diff_restr_imp = config.backtrack['diff_restr']
c_imp          = config.backtrack['c_suff_decr']
k_max_imp      = config.backtrack['k_max']
merit_tau_imp  = config.test['merit_tau']
idx_draw_imp   = config.test['idx_draw']
gamma_imp      = config.physical['gamma']


# ----------------------------------------------------------------------------------------------------------------------
# Simulation
# ----------------------------------------------------------------------------------------------------------------------

def solve_local_model(mesh=0, max_iterations = 1, method = 1, merit_tau = merit_tau_imp, pressure_colour = False, colour_interval = 1, save = False, load = False, iterations_start=1, save_interval=10, k_max=k_max_imp):
    """ Solve foam system locally (either local or quasi-global) for a given number of iterations
        OpenMesh TriMesh mesh: dual foam mesh with property attributes
        int max_iterations: number of iterations
        int method: 1 for local, 2 for quasi_global
        int merit_tau: face index for junction - draw merit as function of tau from 0 to 1
        bool pressure_colour: draw additional foam with colourmap indicating pressures in cells
        int colour_interval: if pressure_colour, draw for every colour_interval iterations
        bool save: save mesh at end of run
        bool load: load mesh instead of using initial mesh
        int iterations_start: iteration counter corresponding to the iteration, the saved mesh was saved at, if load mesh
        int save_interval: interval, for which mesh is saved, if save == True
        int k_max: maximum iterations for inner loop for junction variable updates
    """

    # Progress bar
    print("At iteration: "),

    # Initialise index counter, if a test plot of the search direction is desired
    if merit_tau: fh_idx = 0

    # Identify folder name
    if method == 1:
        folder_name = 'local'
    elif method == 2:
        folder_name = 'quasi_global'

    # Load saved mesh or use mesh parsed to function
    if load:
        mesh = fp.load_foam_mesh(method_name=folder_name, iteration=iterations_start)
    elif not load and mesh == 0:
        raise ValueError('Mesh not specified. Please specify iterations_start according to file to be loaded or pass existing mesh to function')

    # Add face connectivity property
    fp.add_face_connectivity(mesh)

    # Initialise k-counter list
    k_counter_list = []

    # Simulate
    for i in range(max_iterations):
        if load: i = i + iterations_start

        # Add needed attribute properties to all junctions for making topology changes and update face connectivity property
        fp.update_vertex_positions(mesh)
        fp.add_face_connectivity(mesh)

        # Save mesh
        if save and i% save_interval == 0: fp.save_foam_mesh(mesh, iteration=i, method_name=folder_name)

        # Draw mesh
        if i%10 == 0: fv.draw_primary_mesh(mesh, title=folder_name + '/svg/simulation_iteration_' + str(i).zfill(3), pressure_colour = False, draw_dual=False)
        if pressure_colour and i%colour_interval==0: fv.draw_primary_mesh(mesh, title=folder_name + '/svg/simulation_iteration_pressure_' + str(i).zfill(3), pressure_colour = True)

        # Add to progress bar
        print(str(i+1) + " "),

        # Perform topology changes
        ft.perform_topology_changes(mesh)
        fp.update_boundary_pressure(mesh)
        fp.add_face_connectivity(mesh)

        # Create update variable property, if method 2
        if method == 2:
            global update_variables_handle                                                                              # Define and make property handle global, so that it is accessible from other functions
            update_variables_handle = FPropHandle()
            mesh.add_property(update_variables_handle, "fprop_float")

        # Initialise k-counter and face counter
        k_counter    = 0
        face_counter = 0

        # Iterate over faces
        for fh in mesh.faces():
            if not mesh.is_boundary(fh):
                ## VARIABLES
                # Attributes
                attributes = mesh.property(fp.face_connectivity_handle, fh)
                idx_i, idx_j, idx_k, fh_K, fh_J, fh_I, vh_i, vh_j, vh_k = attributes
                thetas = fp.calc_plateau_angles(mesh, fh)

                # Pressure values
                p_i = mesh.property(fp.pressure_handle, vh_i)
                p_j = mesh.property(fp.pressure_handle, vh_j)
                p_k = mesh.property(fp.pressure_handle, vh_k)

                # Define targets and current values
                theta_i_T  = 2.0 * np.pi / 3.0
                theta_j_T  = 2.0 * np.pi / 3.0
                theta_i    = thetas[0, 0]
                theta_j    = thetas[1, 0]
                areas_temp = fp.calc_areas_bubble(mesh, vh_i)
                area_i_T   = areas_temp[1]
                area_i     = areas_temp[0]
                areas_temp = fp.calc_areas_bubble(mesh, vh_j)
                area_j_T   = areas_temp[1]
                area_j     = areas_temp[0]
                if not mesh.is_boundary(vh_k):
                    areas_temp = fp.calc_areas_bubble(mesh, vh_k)
                    area_k_T   = areas_temp[1]
                    area_k     = areas_temp[0]

                # Create target vector and current values
                if not mesh.is_boundary(vh_k):
                    objective = np.array([[theta_i_T], [theta_j_T], [area_i_T], [area_j_T], [area_k_T]])
                    current   = np.array([[theta_i], [theta_j], [area_i], [area_j], [area_k]])
                else:
                    objective = np.array([[theta_i_T], [theta_j_T], [area_i_T], [area_j_T]])
                    current   = np.array([[theta_i], [theta_j], [area_i], [area_j]])

                minus_F   = objective - current

                # # Make area differences 0, if they are too small to prevent overflow ehen calculating delta_alpha
                # if abs(areas_i) < 1e-10: areas_i = 0
                # if abs(areas_j) < 1e-10: areas_j = 0
                # if abs(areas_k) < 1e-10: areas_k = 0

                # Count number of neighboring bubbles (also counting world bubble if on boundary)
                neighbours_i = np.sum(np.array([1 for vh0 in mesh.vv(vh_i) if mesh.is_boundary(vh0) == False]))
                neighbours_j = np.sum(np.array([1 for vh0 in mesh.vv(vh_j) if mesh.is_boundary(vh0) == False]))
                if not mesh.is_boundary(vh_k): neighbours_k = np.sum(np.array([1 for vh0 in mesh.vv(vh_k) if mesh.is_boundary(vh0) == False]))
                is_boundary_i = True if np.sum(np.array([1 for vh0 in mesh.vv(vh_i) if mesh.is_boundary(vh0)]))>0 else False
                is_boundary_j = True if np.sum(np.array([1 for vh0 in mesh.vv(vh_j) if mesh.is_boundary(vh0)]))>0 else False
                if not mesh.is_boundary(vh_k): is_boundary_k = True if np.sum(np.array([1 for vh0 in mesh.vv(vh_k) if mesh.is_boundary(vh0)])) > 0 else False
                if is_boundary_i: neighbours_i += 1
                if is_boundary_j: neighbours_j += 1
                if not mesh.is_boundary(vh_k) and is_boundary_k: neighbours_k += 1

                # Initialise delta_alpha
                delta_alpha_total = np.array([[0], [0], [0], [0], [0]])


                # Initialise iterator
                j = 0

                ## SOLVE MODEL AND UPDATE VALUES
                # Initiate k loop
                while 0.5 * np.dot(np.transpose(minus_F),minus_F)[0][0]>1e-6:
                    # Update iteration counter
                    j += 1

                    # Solve system of linear equations nabla_F * delta_alpha = -F(alpha)
                    if mesh.is_boundary(vh_k):                                                                          # If on boundary
                        nabla_F = fp.calc_jacobian(mesh, fh, delta_alpha=delta_alpha_total)
                        #minus_F = np.array([theta_i, theta_j, area_i, area_j])
                        #minus_F.shape = (len(minus_F), 1)
                        delta_alpha_start = np.dot(np.linalg.pinv(nabla_F), minus_F)
                        delta_alpha_start = np.append(delta_alpha_start, [0.0])
                        delta_alpha_start.shape = (5, 1)
                    else:                                                                                               # If inside foam
                        nabla_F = fp.calc_jacobian(mesh, fh, delta_alpha=delta_alpha_total)
                        #minus_F       = np.array([theta_i, theta_j, area_i, area_j, area_k])
                        #minus_F.shape = (len(minus_F), 1)
                        delta_alpha_start  =  np.dot(np.linalg.pinv(nabla_F), minus_F)

                    # Backtrack to find tau. OBS: WE USE FIXED TAU SO FAR, JUST IN THE TESTING PHASE
                    tau = backtrack_local(mesh, fh, delta_alpha_start, minus_F, nabla_F, attributes, beta_imp, update_variables=delta_alpha_total)

                    # Update delta_alpha
                    delta_alpha        = tau * delta_alpha_start
                    delta_alpha_total  = delta_alpha_total + delta_alpha

                    # Update variables temporarily
                    thetas  = fp.calc_plateau_angles(mesh, fh, delta_alpha=delta_alpha_total)
                    theta_i = thetas[0, 0]
                    theta_j = thetas[1, 0]
                    area_i  = fp.calc_areas_bubble(mesh, vh_i, delta_alpha=delta_alpha_total)[0]
                    area_j  = fp.calc_areas_bubble(mesh, vh_j, delta_alpha=delta_alpha_total)[0]
                    if not mesh.is_boundary(vh_k): area_k  = fp.calc_areas_bubble(mesh, vh_k, delta_alpha=delta_alpha_total)[0]

                    if not mesh.is_boundary(vh_k):
                        current = np.array([[theta_i], [theta_j], [area_i], [area_j], [area_k]])
                    else:
                        current = np.array([[theta_i], [theta_j], [area_i], [area_j]])

                    minus_F   = objective - current

                    if tau < 1e-10 or j > k_max: break


                # Update counters
                k_counter    += i
                face_counter += 1

                # Update x, y, p_i, p_j, p_k
                xy       = mesh.property(fp.mid_point_handle, fh)
                delta_xy = TriMesh.Point(delta_alpha_total[0][0], delta_alpha_total[1][0], 0)

                # Plot tau and search direction
                if mesh.is_boundary(vh_k) and merit_tau: fh_idx += 1
                if mesh.is_boundary(vh_k) and merit_tau and fh_idx == idx_draw_imp:
                    draw_merit(mesh, fh, delta_alpha_total)

                # Update variables differently, depending on whether method 1 or 2
                if method == 1:
                    mesh.set_property(fp.mid_point_handle, fh, xy + delta_xy)

                    if abs(delta_alpha_total[2][0]) > 1e-10:
                        mesh.set_property(fp.pressure_handle, vh_i, p_i + delta_alpha_total[2][0]/neighbours_i)
                    else: pass

                    if abs(delta_alpha_total[3][0]) > 1e-10:
                        mesh.set_property(fp.pressure_handle, vh_j, p_j + delta_alpha_total[3][0]/neighbours_j)
                    else: pass

                    if abs(delta_alpha_total[4][0]) > 1e-10 and mesh.is_boundary(vh_k) == False:
                        mesh.set_property(fp.pressure_handle, vh_k, p_k + delta_alpha_total[4][0]/neighbours_k)
                    else: pass

                    # Update vertex positions to use in coor_bubble
                    fp.update_vertex_positions(mesh)

                elif method == 2:
                    if abs(delta_alpha_total[2][0]) > 1e-10: delta_pi = delta_alpha_total[2][0]/neighbours_i
                    else: delta_pi = 0

                    if abs(delta_alpha_total[3][0]) > 1e-10: delta_pj = delta_alpha_total[3][0]/neighbours_j
                    else: delta_pj = 0

                    if abs(delta_alpha_total[4][0]) > 1e-10 and mesh.is_boundary(vh_k) == False: delta_pk = delta_alpha_total[4][0]/neighbours_k
                    else: delta_pk = 0

                    # Save delta_alpha as property
                    mesh.set_property(update_variables_handle, fh, [delta_xy, delta_pi, delta_pj, delta_pk])

        # If method 2, then update all variables
        if method == 2:
            for fh in mesh.faces():
                if not mesh.is_boundary(fh):
                    # Fetch alpha
                    _, _, _, _, _, _, vh_i, vh_j, vh_k = mesh.property(fp.face_connectivity_handle, fh)
                    xy  = mesh.property(fp.mid_point_handle, fh)
                    p_i = mesh.property(fp.pressure_handle, vh_i)
                    p_j = mesh.property(fp.pressure_handle, vh_j)
                    p_k = mesh.property(fp.pressure_handle, vh_k)

                    # Fetch delta_alpha
                    delta_xy, delta_pi, delta_pj, delta_pk = mesh.property(update_variables_handle, fh)

                    # Update variables
                    mesh.set_property(fp.mid_point_handle, fh, xy + delta_xy)
                    mesh.set_property(fp.pressure_handle, vh_i, p_i + delta_pi)
                    mesh.set_property(fp.pressure_handle, vh_j, p_j + delta_pj)
                    mesh.set_property(fp.pressure_handle, vh_k, p_k + delta_pk)

        k_counter_list.append(k_counter/face_counter)

    # Do last visualisations
    fv.draw_primary_mesh(mesh, title = folder_name + '/svg/simulation_iteration_' + str(i + 1).zfill(3), pressure_colour=False)
    if pressure_colour and (i+1)%colour_interval==0: fv.draw_primary_mesh(mesh, title=folder_name + '/svg/simulation_iteration_pressure_' + str(i + 1).zfill(3), pressure_colour=True)
    if save and (i+1)%save_interval==0: fp.save_foam_mesh(mesh, iteration=i+1, method_name=folder_name)

    # Save average k-counter
    #np.save('saved/k_counter_'+folder_name, k_counter_list)


def solve_global_model(mesh, max_iterations = 1, pressure_colour = False, colour_interval = 1, save = False, load = False, iterations_start=1, save_interval=10):
    """ Solve global system for a given number of iterations
        OpenMesh TriMesh mesh: dual foam mesh with property attributes
        int max_iterations: number of iterations
        bool pressure_colour: draw additional foam with colourmap indicating pressures in cells
        int colour_interval: if pressure_colour, draw for every colour_interval iterations
    """

    # Progress bar
    print("At iteration: "),

    # Load saved mesh or use mesh parsed to function
    if load:
        mesh = fp.load_foam_mesh(method_name='global', iteration=iterations_start)
    elif not load and mesh == 0:
        raise ValueError('Mesh not specified. Please specify iterations_start according to file to be loaded or pass existing mesh to function')

    # Simulate
    for i in range(max_iterations):
        if load: i = i + iterations_start

        # Add needed attribute properties to all junctions for making topology changes
        fp.add_face_connectivity(mesh)
        fp.update_vertex_positions(mesh)

        # Draw mesh
        fv.draw_primary_mesh(mesh, title='global/svg/simulation_iteration_' + str(i).zfill(3), pressure_colour = False)
        if pressure_colour and i%colour_interval==0: fv.draw_primary_mesh(mesh, title='global/svg/simulation_iteration_pressure_' + str(i).zfill(3), pressure_colour = True)

        # Save mesh
        if save and i% save_interval == 0: fp.save_foam_mesh(mesh, iteration=i, method_name='global')

        # Add to progress bar
        print(str(i+1) + " "),

        # Perform topology changes
        ft.perform_topology_changes(mesh)
        fp.add_face_connectivity(mesh)


        ## UPDATE VARIABLES
        # Update positions
        number_junctions = 0
        for fh0 in mesh.faces():
            if not mesh.is_boundary(fh0):
                number_junctions += 1

        # Solve system
        nabla_F, minus_F, vertex_mapping, face_mapping = fp.calc_global_model(mesh)
        #print np.linalg.matrix_rank(nabla_F), len(minus_F)

        #for x in minus_F: print x
        delta_alpha = np.dot(np.linalg.pinv(nabla_F), minus_F)

        # Backtrack to find tau
        taus = backtrack_global(mesh, delta_alpha, face_mapping, vertex_mapping, nabla_F, minus_F)
        delta_alpha_new = taus * delta_alpha

        for fh in mesh.faces():
            if mesh.is_boundary(fh) == False:
                x        = np.where(face_mapping == fh.idx())[0][0]
                y        = x + number_junctions
                xy       = mesh.property(fp.mid_point_handle, fh)
                dx       = delta_alpha_new[x][0] if abs(delta_alpha_new[x][0]) > 1e-10 else 0.0
                dy       = delta_alpha_new[y][0] if abs(delta_alpha_new[y][0]) > 1e-10 else 0.0
                delta_xy = TriMesh.Point(dx, dy, 0)
                mesh.set_property(fp.mid_point_handle, fh, xy + delta_xy)

        for vh in mesh.vertices():
            if mesh.is_boundary(vh) == False:
                p       = mesh.property(fp.pressure_handle, vh)
                p_idx   = number_junctions * 2 + np.where(vertex_mapping == vh.idx())[0][0]
                delta_p = delta_alpha_new[p_idx][0]

                if abs(delta_p) > 1e-10:
                    mesh.set_property(fp.pressure_handle, vh, p + delta_p)
                else:
                    pass

    fv.draw_primary_mesh(mesh, title = 'global/svg/simulation_iteration_' + str(i + 1).zfill(3), pressure_colour=False)
    if pressure_colour and (i+1)%colour_interval==0: fv.draw_primary_mesh(mesh, title='global/svg/simulation_iteration_pressure_' + str(i + 1).zfill(3), pressure_colour=True)
    if save and (i+1)%save_interval==0: fp.save_foam_mesh(mesh, iteration=i+1, method_name='global')


def backtrack_local(mesh, fh, delta_alpha, minus_F, jacobian, attributes, beta = beta_imp, diff_restr = diff_restr_imp, update_variables=[]):
    """ Backtrack to find a tau, making tau * delta_alpha a valid step
        OpenMesh TriMesh mesh: dual foam mesh with property attributes
        TriMesh FaceHandle fh: for junction in question
        NumPy array delta_alpha: search direction
        NumPy array attributes: attributes for fh
        double beta: update factor for tau
        double diff_restr: restriction for how much the pressure can change through diffusion in this iteration.
        NumPy array update_variables: If x, y and p's need to be updated
    """


    # INITIALISE CONSTANTS
    tau             = 1
    idx_i, idx_j, idx_k, fh_K, fh_J, fh_I, vh_i, vh_j, vh_k = attributes
    x_0             = np.array([mesh.property(fp.mid_point_handle, fh)[0], mesh.property(fp.mid_point_handle, fh)[1]])
    x_K             = np.array([mesh.property(fp.mid_point_handle, fh_K)[0], mesh.property(fp.mid_point_handle, fh_K)[1]])  # Junction oppposite of K
    x_J             = np.array([mesh.property(fp.mid_point_handle, fh_J)[0], mesh.property(fp.mid_point_handle, fh_J)[1]])  # Junction oppposite of J
    x_I             = np.array([mesh.property(fp.mid_point_handle, fh_I)[0], mesh.property(fp.mid_point_handle, fh_I)[1]])  # Junction oppposite of I
    p_i             = mesh.property(fp.pressure_handle, vh_i)
    p_j             = mesh.property(fp.pressure_handle, vh_j)
    p_k             = mesh.property(fp.pressure_handle, vh_k)

    # Update variable positions, if desired
    if len(update_variables) != 0:
        x_0 = x_0 + np.array([update_variables[0][0], update_variables[1][0]])
        p_i = p_i + update_variables[2][0]
        p_j = p_j + update_variables[3][0]
        if not mesh.is_boundary(vh_k): p_k = p_k + update_variables[4][0]

    # Calculate initial surface energy
    energy_start    = np.sum(np.array([fp.calc_film_length(x_0, x_K, p_i - p_j),
                                       fp.calc_film_length(x_0, x_J, p_i - p_k),
                                       fp.calc_film_length(x_0, x_I, p_j - p_k)]))                                      # Potential energy (proportional to film length)



    ## GEOMETRIC CONDITIONS
    # Initialise variables
    x_0_next = x_0 + tau * np.array([delta_alpha[0][0], delta_alpha[1][0]])

    # Backtrack
    while not is_inside_bubble(mesh, x_0_next, vh_i, vh_j, vh_k) and mesh.is_boundary(vh_k) == False:
        tau *= beta
        x_0_next = x_0 + tau * np.array([delta_alpha[0][0], delta_alpha[1][0]])


    ## DECREASE IN POTENTIAL ENERGY
    # Initialise variables
    x_0_next    = x_0 + tau * np.array([delta_alpha[0][0], delta_alpha[1][0]])
    delta_p_ij  = p_i - p_j
    delta_p_ik  = p_i - p_k
    delta_p_jk  = p_j - p_k
    energy_next = np.sum(np.array([fp.calc_film_length(x_0_next, x_K, delta_p_ij),
                                   fp.calc_film_length(x_0_next, x_J, delta_p_ik),
                                   fp.calc_film_length(x_0_next, x_I, delta_p_jk)]))

    # Backtrack
    while energy_next > energy_start:
        tau *= beta
        if tau < 1e-10: return 0

        x_0_next    = x_0 + tau * np.array([delta_alpha[0][0], delta_alpha[1][0]])
        energy_next = np.sum(np.array([fp.calc_film_length(x_0_next, x_K, delta_p_ij),
                                       fp.calc_film_length(x_0_next, x_J, delta_p_ik),
                                       fp.calc_film_length(x_0_next, x_I, delta_p_jk)]))

    ## SUFFICIENT DECREASE
    # backtrack
    while merit(mesh, fh, tau * delta_alpha, update_variables) > merit_gradient(delta_alpha, jacobian, minus_F, tau, c_imp):
        #print 'jegerher'
        tau *= beta
        if tau < 1e-10: return 0

    ## RESTRICT DIFFUSION
    while abs(tau * delta_alpha[2][0]/p_i) > diff_restr or abs(tau * delta_alpha[3][0]/p_j) > diff_restr or abs(tau * delta_alpha[4][0]/p_k) > diff_restr:
        tau *= beta
        if tau < 1e-10: return 0

    ## CHECK, IF FILMS HAVE CROSSED EACH OTHER DUE TO TOO LARGE PRESSURE VALUES
    angles = fp.calc_plateau_angles(mesh,fh,tau*delta_alpha)[:,0]
    while np.any(angles<0):
        tau *= beta
        if tau < 1e-10: return 0

        angles = fp.calc_plateau_angles(mesh, fh, tau * delta_alpha)[:, 0]

    return tau


def backtrack_global(mesh, delta_alpha, face_mapping, vertex_mapping, jacobian, minus_F):
    """ Perform backtrack on global delta_alpha
        OpenMesh TriMesh mesh: dual foam mesh
        NumPy array delta_alpha: search direction
        NumPy array face_mapping: face mapping (iterator vs. index)
        NumPy array vertex_mapping: vertex mapping (iterator vs. index)
        bool global_backtrack: one tau for whole vector, or individual taus
        double beta: update fraction
        double diff_restr: limit on diffusion
        double gamma: surface tension
        return tau/taus: damping of search direction - one or individual
    """
    ## CONSTANTS
    # Count number of junctions and number of cells
    number_junctions = 0
    for fh0 in mesh.faces():
        if not mesh.is_boundary(fh0):
            number_junctions += 1

    # INITIALISE VARIABLES
    taus = np.zeros([len(delta_alpha), 1])

    # RUN VERIFICATION LOOP
    # Check, if inside cell and that curvature radii do not get smaller than half the radius between the junctions (junction position part)
    for fh in mesh.faces():
        if not mesh.is_boundary(fh):
            attributes = mesh.property(fp.face_connectivity_handle, fh)
            idx_i, idx_j, idx_k, fh_K, fh_J, fh_I, vh_i, vh_j, vh_k = attributes

            # Count number of neighbours
            neighbours_i = np.sum(np.array([1 for vh0 in mesh.vv(vh_i) if mesh.is_boundary(vh0) == False]))
            neighbours_j = np.sum(np.array([1 for vh0 in mesh.vv(vh_j) if mesh.is_boundary(vh0) == False]))
            if not mesh.is_boundary(vh_k): neighbours_k = np.sum(np.array([1 for vh0 in mesh.vv(vh_k) if mesh.is_boundary(vh0) == False]))

            # Find indices
            x_idx    = np.where(face_mapping == fh.idx())[0][0]                                                         # x in2dex in delta_alpha
            y_idx    = x_idx + number_junctions
            pi_idx   = number_junctions * 2 + np.where(vertex_mapping == idx_i)[0][0]
            pj_idx   = number_junctions * 2 + np.where(vertex_mapping == idx_j)[0][0]
            if not mesh.is_boundary(vh_k): pk_idx = number_junctions * 2 + np.where(vertex_mapping == idx_k)[0][0]
            theta_i   = np.where(face_mapping == fh.idx())[0][0] * 2
            theta_j   = theta_i + 1
            area_i    = np.where(vertex_mapping == idx_i)[0][0] + number_junctions * 2
            area_j    = np.where(vertex_mapping == idx_j)[0][0] + number_junctions * 2
            if not mesh.is_boundary(vh_k): area_k = np.where(vertex_mapping == idx_k)[0][0] + number_junctions * 2

            # Delta alpha
            if not mesh.is_boundary(vh_k):
                delta_alpha_junction = np.array([[delta_alpha[x_idx][0]], [delta_alpha[y_idx][0]], [delta_alpha[pi_idx][0]], [delta_alpha[pj_idx][0]], [delta_alpha[pk_idx][0]]])
                minus_F_junction     = np.array([[minus_F[theta_i][0]], [minus_F[theta_j][0]], [minus_F[area_i][0]], [minus_F[area_j][0]], [minus_F[area_k][0]]])
                jacobian_junction    = np.array([[jacobian[theta_i][x_idx], jacobian[theta_i][y_idx], jacobian[theta_i][pi_idx], jacobian[theta_i][pj_idx], jacobian[theta_i][pk_idx]],
                                                 [jacobian[theta_j][x_idx], jacobian[theta_j][y_idx], jacobian[theta_j][pi_idx], jacobian[theta_j][pj_idx], jacobian[theta_j][pk_idx]],
                                                 [jacobian[area_i][x_idx], jacobian[area_i][y_idx], jacobian[area_i][pi_idx], jacobian[area_i][pj_idx], jacobian[area_i][pk_idx]],
                                                 [jacobian[area_j][x_idx], jacobian[area_j][y_idx], jacobian[area_j][pi_idx], jacobian[area_j][pj_idx], jacobian[area_j][pk_idx]],
                                                 [jacobian[area_k][x_idx], jacobian[area_k][y_idx], jacobian[area_k][pi_idx], jacobian[area_k][pj_idx], jacobian[area_k][pk_idx]]])
            else:
                delta_alpha_junction = np.array([[delta_alpha[x_idx][0]], [delta_alpha[y_idx][0]], [delta_alpha[pi_idx][0]], [delta_alpha[pj_idx][0]], [0]])
                minus_F_junction     = np.array([[minus_F[theta_i][0]], [minus_F[theta_j][0]], [minus_F[area_i][0]], [minus_F[area_j][0]]])
                jacobian_junction    = np.array([[jacobian[theta_i][x_idx], jacobian[theta_i][y_idx], jacobian[theta_i][pi_idx], jacobian[theta_i][pj_idx]],
                                                 [jacobian[theta_j][x_idx], jacobian[theta_j][y_idx], jacobian[theta_j][pi_idx], jacobian[theta_j][pj_idx]],
                                                 [jacobian[area_i][x_idx], jacobian[area_i][y_idx], jacobian[area_i][pi_idx], jacobian[area_i][pj_idx]],
                                                 [jacobian[area_j][x_idx], jacobian[area_j][y_idx], jacobian[area_j][pi_idx], jacobian[area_j][pj_idx]]])


            # Perform local backtrack on this particular delta_alpha
            tau = backtrack_local(mesh, fh, delta_alpha_junction, minus_F_junction, jacobian_junction, attributes, beta=beta_imp, diff_restr=diff_restr_imp)

            # Update tau values in large tau vector
            taus[x_idx]  = tau
            taus[y_idx]  = tau
            taus[pi_idx] += tau / neighbours_i
            taus[pj_idx] += tau / neighbours_j
            if not mesh.is_boundary(vh_k): taus[pk_idx] += tau / neighbours_k

    return taus


def is_inside_bubble(mesh, x_0, vh_i, vh_j, vh_k):
    """ Determine, whether a certain junction point is inside or the three adjacent bubble cells
        NumPy array x_0: junction coordinate
        OpenMesh TriMesh mesh: dual foam mesh
        OpenMesh vertex handles v_i, v_j, v_k: vertex handles of face containing junction of interest
        return bool: True, if inside, False, if outside
    """

    # Build lists of faces containing the junctions that border the three cells
    faces_i = np.array([fh for fh in mesh.vf(vh_i) if mesh.is_boundary(fh) == False])
    faces_j = np.array([fh for fh in mesh.vf(vh_j) if mesh.is_boundary(fh) == False])
    faces_k = np.array([fh for fh in mesh.vf(vh_k) if mesh.is_boundary(fh) == False])

    # Create a path object representing the border of straight lines for each cell
    cell_i    = mplPath.Path(np.array([[mesh.property(fp.mid_point_handle, fh0)[0], mesh.property(fp.mid_point_handle, fh0)[1]] for fh0 in faces_i]))
    cell_j    = mplPath.Path(np.array([[mesh.property(fp.mid_point_handle, fh0)[0], mesh.property(fp.mid_point_handle, fh0)[1]] for fh0 in faces_j]))
    cell_k    = mplPath.Path(np.array([[mesh.property(fp.mid_point_handle, fh0)[0], mesh.property(fp.mid_point_handle, fh0)[1]] for fh0 in faces_k]))
    cell_dual = mplPath.Path(np.array([[mesh.point(vh_i)[0], mesh.point(vh_i)[1]], [mesh.point(vh_j)[0], mesh.point(vh_j)[1]], [mesh.point(vh_k)[0], mesh.point(vh_k)[1]]]))

    # Check, if x_0 is inside any of the cells
    inside_i            = cell_i.contains_point(x_0)
    inside_j            = cell_j.contains_point(x_0)
    inside_k            = cell_k.contains_point(x_0)
    inside_world_bubble = cell_dual.contains_point(x_0)                                                                 # If junction is on boundary, do work-around to see, if junction is moving towards world bubble. Can be improved. Now just checking, if inside dual face

    if inside_i or inside_j or inside_k or inside_world_bubble: return True
    else: return False


def curvature_radius_too_small(mesh, x_0, vh_i, vh_j, vh_k):
    """ Determine, whether a certain junction point is inside or the three adjacent bubble cells
        NumPy array x_0: junction coordinate
        OpenMesh TriMesh mesh: dual foam mesh
        OpenMesh vertex handles v_i, v_j, v_k: vertex handles of face containing junction of interest
        return bool: True, if inside, False, if outside
    """

    # Build lists of faces containing the junctions that border the three cells
    faces_i = np.array([fh for fh in mesh.vf(vh_i) if mesh.is_boundary(fh) == False])
    faces_j = np.array([fh for fh in mesh.vf(vh_j) if mesh.is_boundary(fh) == False])
    faces_k = np.array([fh for fh in mesh.vf(vh_k) if mesh.is_boundary(fh) == False])

    # Create a path object representing the border of straight lines for each cell
    cell_i    = mplPath.Path(np.array([[mesh.property(fp.mid_point_handle, fh0)[0], mesh.property(fp.mid_point_handle, fh0)[1]] for fh0 in faces_i]))
    cell_j    = mplPath.Path(np.array([[mesh.property(fp.mid_point_handle, fh0)[0], mesh.property(fp.mid_point_handle, fh0)[1]] for fh0 in faces_j]))
    cell_k    = mplPath.Path(np.array([[mesh.property(fp.mid_point_handle, fh0)[0], mesh.property(fp.mid_point_handle, fh0)[1]] for fh0 in faces_k]))
    cell_dual = mplPath.Path(np.array([[mesh.point(vh_i)[0], mesh.point(vh_i)[1]], [mesh.point(vh_j)[0], mesh.point(vh_j)[1]], [mesh.point(vh_k)[0], mesh.point(vh_k)[1]]]))

    # Check, if x_0 is inside any of the cells
    inside_i            = cell_i.contains_point(x_0)
    inside_j            = cell_j.contains_point(x_0)
    inside_k            = cell_k.contains_point(x_0)
    inside_world_bubble = cell_dual.contains_point(x_0)                                                                 # If junction is on boundary, do work-around to see, if junction is moving towards world bubble. Can be improved. Now just checking, if inside dual face

    if inside_i or inside_j or inside_k or inside_world_bubble: return True
    else: return False


def merit_gradient(delta_alpha, jacobian, minus_F, tau, c):
    """ Calculate upper boundary of slope line of gradient of merit function at tau=0
        NumPy array delta_alpha: search direction
        NumPy array jacobian: jacobian matrix of F
        NumPy array minus_F: jacobian matrix of F
        double tau: step length
        double c: damping factor of slope
        return dtheta_dtau: merit of straight line at tau*delta_alpha of line with slope dtheta_dtau_tau=0 * c
    """

    F = -minus_F

    if len(F)==4:
        delta_alpha = delta_alpha[0:4]

    theta_0 = 0.5 * np.dot(np.transpose(F), F)[0][0]

    dtheta_dtau = theta_0 + c * tau * np.dot(np.dot(np.transpose(F), jacobian), delta_alpha)[0][0]

    return dtheta_dtau


def merit(mesh, fh, delta_alpha, update_variables=[]):
    """ Calculate merit of F
        OpenMesh TriMesh mesh: dual foam mesh
        OpenMesh FaceHandle fh: face of junction
        NumPy array delta_alpha: search direction
        return theta: merit function at delta_alpha
    """

    # Prepare variables
    idx_i, idx_j, idx_k, fh_K, fh_J, fh_I, vh_i, vh_j, vh_k = mesh.property(fp.face_connectivity_handle, fh)
    x_0      = np.array([mesh.property(fp.mid_point_handle, fh)[0], mesh.property(fp.mid_point_handle, fh)[1]])
    x_K      = np.array([mesh.property(fp.mid_point_handle, fh_K)[0], mesh.property(fp.mid_point_handle, fh_K)[1]])     # Junction oppposite of K
    x_J      = np.array([mesh.property(fp.mid_point_handle, fh_J)[0], mesh.property(fp.mid_point_handle, fh_J)[1]])     # Junction oppposite of J
    x_I      = np.array([mesh.property(fp.mid_point_handle, fh_I)[0], mesh.property(fp.mid_point_handle, fh_I)[1]])     # Junction oppposite of I
    p_i      = mesh.property(fp.pressure_handle, vh_i)
    p_j      = mesh.property(fp.pressure_handle, vh_j)
    p_k      = mesh.property(fp.pressure_handle, vh_k)

    # Update variable positions, if desired
    if len(update_variables) != 0:
        x_0 = x_0 + np.array([update_variables[0][0], update_variables[1][0]])
        p_i = p_i + update_variables[2][0]
        p_j = p_j + update_variables[3][0]
        if not mesh.is_boundary(vh_k): p_k = p_k + update_variables[4][0]

    # Search direction
    delta_xy = np.array([delta_alpha[0][0], delta_alpha[1][0]])
    x_0_new  = x_0 + delta_xy
    p_i_new  = p_i + delta_alpha[2][0]
    p_j_new  = p_j + delta_alpha[3][0]
    p_k_new  = p_k + delta_alpha[4][0]

    # Prepare areas
    area_i_target = fp.calc_areas_bubble(mesh, vh_i)[1]
    area_j_target = fp.calc_areas_bubble(mesh, vh_j)[1]
    if not mesh.is_boundary(vh_k): area_k_target = fp.calc_areas_bubble(mesh, vh_k)[1]
    else: area_k_target = 0

    area_i_new      = fp.calc_areas_bubble(mesh, vh_i, np.array([fh.idx(), idx_i, idx_j, idx_k]), delta_alpha, calc_target = False)
    area_j_new      = fp.calc_areas_bubble(mesh, vh_j)[0]#, np.array([fh.idx(), idx_i, idx_j, idx_k]), delta_alpha, calc_target = False)
    if mesh.is_boundary(vh_k) == False: area_k_new = fp.calc_areas_bubble(mesh, vh_k, np.array([fh.idx(), idx_i, idx_j, idx_k]), delta_alpha, calc_target = False)
    else: area_k_new = 0

    # Calculate entries in F(alpha + delta_alpha)
    F_theta_i  = 2 * np.pi / 3 - fp.calc_angle(x_K, x_J, x_0_new, p_j_new - p_i_new, p_k_new - p_i_new)
    F_theta_j  = 2 * np.pi / 3 - fp.calc_angle(x_I, x_K, x_0_new, p_k_new - p_j_new, p_i_new - p_j_new)
    F_area_i   = area_i_target - area_i_new
    F_area_j   = area_j_target - area_j_new
    F_area_k   = area_k_target - area_k_new

    F = np.array([F_theta_i, F_theta_j, F_area_i, F_area_j, F_area_k])

    theta = 0.5 * np.dot(np.transpose(F), F)

    return theta


def draw_merit(mesh, fh, delta_alpha):
    """ Make plot of merit function along delta alpha with tau ranging from 0 to 1
        OpenMesh TriMesh mesh: dual foam mesh
        OpenMesh FaceHandle fh: junction to investigate
    """

    # Prepare variables
    x   =  np.linspace(0, 1, 50)
    F   = [merit(mesh, fh, tau * delta_alpha) for tau in x]

    # Draw slope line
    s = [-0.6*tau+F[0] for tau in x]

    # Plot
    plt.figure()
    plt.plot(x, F, label='Merit function')
    plt.plot(x, s, 'r',label='Sufficient decrease')
    plt.xlabel(r'$\tau$')
    plt.ylabel('$1/2F^TF$')
    plt.legend()
    plt.savefig('images/search_direction.eps', format='eps', dpi=1000)



# ----------------------------------------------------------------------------------------------------------------------
# Experiments
# ----------------------------------------------------------------------------------------------------------------------

def draw_convergence_plot(method, num_meshes = 1000):
    """ Draw convergence plot of data generated with certain method
    """

    if method == 1:
        method_name = 'local'
    elif method == 2:
        method_name = 'quasi_global'
    elif method == 3:
        method_name = 'global'

    # Angle part
    merit_list = []
    merit_list_norm = []                                                                                                # Normalised with number of cells
    for i in range(num_meshes):
        merit     = 0
        num_junc  = 0
        mesh      = fp.load_foam_mesh(method_name, iteration=i)
        fp.add_face_connectivity(mesh)

        for fh in mesh.faces():
            if not mesh.is_boundary(fh):
                # Update junction counter
                num_junc += 1

                # Fetch attrisbutes
                _, _, _, fh_K, fh_J, fh_I, vh_i, vh_j, vh_k = mesh.property(fp.face_connectivity_handle, fh)
                x_0 = np.array([mesh.property(fp.mid_point_handle, fh)[0], mesh.property(fp.mid_point_handle, fh)[1]])
                x_K = np.array([mesh.property(fp.mid_point_handle, fh_K)[0], mesh.property(fp.mid_point_handle, fh_K)[1]])        # Junction oppposite of K
                x_J = np.array([mesh.property(fp.mid_point_handle, fh_J)[0], mesh.property(fp.mid_point_handle, fh_J)[1]])        # Junction oppposite of J
                x_I = np.array([mesh.property(fp.mid_point_handle, fh_I)[0], mesh.property(fp.mid_point_handle, fh_I)[1]])        # Junction oppposite of I
                p_i = mesh.property(fp.pressure_handle, vh_i)
                p_j = mesh.property(fp.pressure_handle, vh_j)
                p_k = mesh.property(fp.pressure_handle, vh_k)

                # Calculate pressure differences
                delta_p_ij = p_j - p_i
                delta_p_ik = p_k - p_i
                delta_p_ji = p_i - p_j
                delta_p_jk = p_k - p_j

                # Calculate angles
                theta_delta_i = 2.0 * np.pi / 3.0 - fp.calc_angle(x_K, x_J, x_0, delta_p_ij, delta_p_ik)
                theta_delta_j = 2.0 * np.pi / 3.0 - fp.calc_angle(x_I, x_K, x_0, delta_p_jk, delta_p_ji)

                # Difference from equilibrium
                merit += theta_delta_i**2 + theta_delta_j**2

        # Area part
        for vh in mesh.vertices():
            if not mesh.is_boundary(vh):
                # Calculate target area and current area
                areas = fp.calc_areas_bubble(mesh, vh)
                delta_areas = areas[1] - areas[0]
                merit += delta_areas**2

        # Append this iteration to lists
        merit_list.append(0.5*merit)
        merit_list_norm.append(0.5 * merit/num_junc)

        print i

    # Save data
    np.save('saved/experiments/merit_' + method_name + '_' + str(i) + '_energy.npy', merit_list)
    np.save('saved/experiments/merit_norm_' + method_name + '_' + str(i) + '_energy.npy', merit_list_norm)

    plt.figure()
    plt.plot(merit_list, '.')
    plt.xlabel('Iteration', fontsize=18)
    plt.ylabel('$1/2F^TF$', fontsize=18)
    plt.savefig('images/merit_' + method_name +'.eps', format='eps', dpi=1000)
    plt.show()

    plt.figure()
    plt.plot(merit_list_norm, '.')
    plt.xlabel('Iteration', fontsize=18)
    plt.ylabel('$1/2F^TF / num_{cells}$', fontsize=18)
    plt.savefig('images/merit_' + method_name +'.eps', format='eps', dpi=1000)
    plt.show()


def draw_film_lengths(method, num_meshes = 1000):
    """ Draw and save total film length in foam and save the data to .npy file
    """

    if method == 1:
        method_name = 'local'
    elif method == 2:
        method_name = 'quasi_global'
    elif method == 3:
        method_name = 'global'

    # Initialise energy list
    energy_foam_list = []
    for i in range(num_meshes):
        print i

        # Initialise and load
        energy_foam = 0
        mesh = fp.load_foam_mesh(method_name, iteration=i)
        fp.add_face_connectivity(mesh)

        # Iterate over faces and count film lengths (count half film lengths, the other end is calculated when reaching the face on the other side)
        for fh in mesh.faces():
            if not mesh.is_boundary(fh):
                # Fetch attributes
                _, _, _, fh_K, fh_J, fh_I, vh_i, vh_j, vh_k = mesh.property(fp.face_connectivity_handle, fh)
                x_0 = np.array([mesh.property(fp.mid_point_handle, fh)[0], mesh.property(fp.mid_point_handle, fh)[1]])
                x_K = np.array([mesh.property(fp.mid_point_handle, fh_K)[0], mesh.property(fp.mid_point_handle, fh_K)[1]])        # Junction oppposite of K
                x_J = np.array([mesh.property(fp.mid_point_handle, fh_J)[0], mesh.property(fp.mid_point_handle, fh_J)[1]])        # Junction oppposite of J
                x_I = np.array([mesh.property(fp.mid_point_handle, fh_I)[0], mesh.property(fp.mid_point_handle, fh_I)[1]])        # Junction oppposite of I
                p_i = mesh.property(fp.pressure_handle, vh_i)
                p_j = mesh.property(fp.pressure_handle, vh_j)
                p_k = mesh.property(fp.pressure_handle, vh_k)

                # Calculate surface energy for junctions
                energy_junction = 0.5 * np.sum(np.array([fp.calc_film_length(x_0, x_K, p_i - p_j),
                                                         fp.calc_film_length(x_0, x_J, p_i - p_k),
                                                         fp.calc_film_length(x_0, x_I, p_j - p_k)]))

                energy_foam += energy_junction


        # Add to list
        energy_foam_list.append(energy_foam)

    # Save data
    np.save('saved/experiments/' + method_name + '_' + str(num_meshes) + '_energy', energy_foam_list)

    # Plot
    plt.figure()
    plt.plot(energy_foam_list, 'k')
    plt.xlabel('Iterations', fontsize=18)
    plt.ylabel('Total film length', fontsize=18)
    plt.savefig('images/energy_'+str(method_name)+'.eps', format='eps', dpi=1000)
    plt.show()

    plt.figure()
    plt.plot(energy_foam_list, 'k')
    plt.xlabel('Iterations', fontsize=18)
    plt.ylabel('Total film length', fontsize=18)
    plt.yscale('log')
    plt.savefig('images/energy_logy_'+str(method_name)+'.eps', format='eps', dpi=1000)
    plt.show()


def draw_mean_cell_area(method, num_meshes = 1000):
    """ Draw and save average cell area in foam and save the data to .npy file
    """

    if method == 1:
        method_name = 'local'
    elif method == 2:
        method_name = 'quasi_global'
    elif method == 3:
        method_name = 'global'

    # Initialise energy list
    average_area_list        = []
    average_area_err_list    = []
    # average_area_nb_list     = []                                                                                       # List without boundary cells
    # average_area_nb_err_list = []

    for i in range(num_meshes):
        print i

        # Initialise and load
        area_foam_list    = []
        # area_foam_nb_list = []
        mesh = fp.load_foam_mesh(method_name, iteration=i)
        fp.add_face_connectivity(mesh)

        # Iterate over faces and count film lengths (count half film lengths, the other end is calculated when reaching the face on the other side)
        for vh in mesh.vertices():
            if not mesh.is_boundary(vh):
                # is_boundary_cell = np.sum(np.array([1 for vh0 in mesh.vv(vh) if mesh.is_boundary(vh0)]))

                # Append to list
                area = fp.calc_areas_bubble(mesh, vh, calc_target = False)
                area_foam_list.append(area)
                # if is_boundary_cell: area_foam_nb_list.append(area)

        # Add to list
        N = len(area_foam_list)
        average_area_list.append(np.sum(area_foam_list)/N)                                                              # Mean
        average_area_err_list.append(np.std(area_foam_list)/np.sqrt(N))                                                 # Uncertainty on mean

        # N_nb = len(area_foam_nb_list)
        # average_area_nb_list.append(np.sum(area_foam_nb_list) / N_nb)                                                   # Mean
        # average_area_nb_err_list.append(np.std(area_foam_nb_list) / np.sqrt(N_nb))                                      # Uncertainty on mean



    # Save data
    np.save('saved/experiments/' + method_name + '_' + str(i) + '_average_area_equal_boundary', average_area_list)
    np.save('saved/experiments/' + method_name + '_' + str(i) + '_average_area_equal_boundary_error', average_area_err_list)
    # np.save('saved/experiments/' + method_name + '_' + str(i) + '_average_area_nb_large', average_area_nb_list)
    # np.save('saved/experiments/' + method_name + '_' + str(i) + '_average_area_nb_large_error', average_area_nb_err_list)

    # # Plot
    # plt.figure()
    # plt.plot(average_area_list)
    # plt.xlabel('Iterations', fontsize=18)
    # plt.ylabel('Mean cell area', fontsize=18)
    # plt.savefig('images/average_area_'+str(method_name)+'.eps', format='eps', dpi=1000)
    # plt.show()


def aboav_weaire(method, frame_number):
    # Define method
    if method == 1:
        method_name = 'local'
    elif method == 2:
        method_name = 'quasi_global'
    elif method == 3:
        method_name = 'global'

    # Load mesh
    mesh = fp.load_foam_mesh(method_name, iteration=frame_number)
    fp.add_face_connectivity(mesh)

    # Initialise lists
    neighbour_mean_list   = np.zeros((8,1))
    neighbour_number_list = np.zeros((8,1))
    err_3  = []
    err_4  = []
    err_5  = []
    err_6  = []
    err_7  = []
    err_8  = []
    err_9  = []
    err_10 = []

    # Iterate over in in-foam cells
    for vh in mesh.vertices():
        is_boundary = np.sum(np.array([1 for vh0 in mesh.vv(vh) if mesh.is_boundary(vh0)])) > 0

        # Iterate over cell neighbours, if cell is not boundary cell
        if not is_boundary:
            neighbours      = 0
            neighbour_sides = 0

            for vh0 in mesh.vv(vh):
                neighbours += 1

                # Count numbers of neighbours of that cell
                is_boundary_1        = np.sum(np.array([1 for vh1 in mesh.vv(vh0) if mesh.is_boundary(vh1)])) > 0
                neighbour_face_count = np.sum(np.array([1 for vh1 in mesh.vv(vh0) if not mesh.is_boundary(vh1)]))
                if is_boundary_1: neighbour_face_count += 1

                # Update total neighbour cell neighbour sides
                neighbour_sides += neighbour_face_count

            # Find mean number of sides of neighbour bubbles
            mean_neighbour_sides = neighbour_sides/neighbours

            if neighbours<11:
                # Update lists
                neighbour_mean_list[neighbours-3][0]   = neighbour_mean_list[neighbours-3][0] + mean_neighbour_sides
                neighbour_number_list[neighbours-3][0] = neighbour_number_list[neighbours-3][0] + 1

                # Add to python lists for error calculations
                if neighbours == 3:
                    err_3.append(mean_neighbour_sides)
                elif neighbours == 4:
                    err_4.append(mean_neighbour_sides)
                elif neighbours == 5:
                    err_5.append(mean_neighbour_sides)
                elif neighbours == 6:
                    err_6.append(mean_neighbour_sides)
                elif neighbours == 7:
                    err_7.append(mean_neighbour_sides)
                elif neighbours == 8:
                    err_8.append(mean_neighbour_sides)
                elif neighbours == 9:
                    err_9.append(mean_neighbour_sides)
                elif neighbours == 10:
                    err_10.append(mean_neighbour_sides)

    # Find means
    for i in range(8):
        if neighbour_number_list[i][0] > 0:
            neighbour_mean_list[i][0] = neighbour_mean_list[i][0] / neighbour_number_list[i][0] * (i+3)

    # Create list of errors on mean
    error_list = [np.std(err_3) / np.sqrt(len(err_3)),
                  np.std(err_4) / np.sqrt(len(err_4)),
                  np.std(err_5) / np.sqrt(len(err_5)),
                  np.std(err_6) / np.sqrt(len(err_6)),
                  np.std(err_7) / np.sqrt(len(err_7)),
                  np.std(err_8) / np.sqrt(len(err_8)),
                  np.std(err_9) / np.sqrt(len(err_9)),
                  np.std(err_10) / np.sqrt(len(err_10))]

    # Save
    np.save('saved/experiments/' + method_name + '_' + str(frame_number) + '_aboav_weaire.npy', neighbour_mean_list)
    np.save('saved/experiments/' + method_name + '_' + str(frame_number) + '_aboav_weaire_error.npy', error_list)



def second_moment(method, frame_number):
    # Define method
    if method == 1:
        method_name = 'local'
    elif method == 2:
        method_name = 'quasi_global'
    elif method == 3:
        method_name = 'global'

    # Load mesh
    mesh = fp.load_foam_mesh(method_name, iteration=frame_number)
    fp.add_face_connectivity(mesh)

    # Initialise lists
    neighbour_face_count_list = []

    # Iterate over in-foam cells
    for vh in mesh.vertices():
        is_boundary = np.sum(np.array([1 for vh0 in mesh.vv(vh) if mesh.is_boundary(vh0)])) > 0

        # Iterate over cell neighbours, if cell is not boundary cell
        if not is_boundary:
            neighbour_face_count = np.sum(np.array([1 for vh0 in mesh.vv(vh) if not mesh.is_boundary(vh0)]))            # number of neighbours
            neighbour_face_count_list.append(neighbour_face_count)


    n_mean     = np.mean(neighbour_face_count_list)
    n_mean_err = np.std(neighbour_face_count_list)/np.sqrt(len(neighbour_face_count_list))
    hist   = np.histogram(neighbour_face_count_list,bins=max(neighbour_face_count_list)-min(neighbour_face_count_list)+1)[0]
    hist_float = []
    for x in hist: hist_float.append(float(x))
    p_n    = hist_float/np.sum(hist)



    print p_n

    # Find mu_2
    mu_2 = 0
    for i in range(max(neighbour_face_count_list)-min(neighbour_face_count_list)+1):
        idx = i + min(neighbour_face_count_list)
        mu_2 += (idx - n_mean)**2 * p_n[i]

    # Do error propagation
    mu_2_err = 0
    for i in range(max(neighbour_face_count_list)-min(neighbour_face_count_list)+1):
        idx = i + min(neighbour_face_count_list)
        mu_2_err += (2 * n_mean - 2 * idx)**2 * p_n[i]**2 * n_mean_err**2
    mu_2_err = np.sqrt(mu_2_err)


    print mu_2, mu_2_err






