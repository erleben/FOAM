import config
import matplotlib.pyplot as plt
import foam_equilibriate as fe
import foam_properties as fp
import foam_visualise as fv
import foam_simulation as fs
from foam_initialise import *



if __name__ == "__main__":
    # Read variables from config-file (Initially for random grid)
    min   = config.length['min']
    max   = config.length['max']
    equal = config.initialise['equal']
    read  = config.open['read']                                                                                         # True if read, False if write
    param = config.test['param_test']


    #-------------------------------------------------------------------------------------------------------------------
    # Initialise mesh: Save or read mesh?
    #-------------------------------------------------------------------------------------------------------------------
    if read:
        foam_mesh = TriMesh()
        read_mesh(foam_mesh, "saved_foam_mesh.obj")
        fp.add_pressure(foam_mesh, equal)
        fp.add_mid_point(foam_mesh)
        fp.update_vertex_positions(foam_mesh)
        fp.add_face_connectivity(foam_mesh)
    else:
        foam_mesh = create_random_mesh(min, max)
        fp.add_pressure(foam_mesh, equal)
        fe.combine_equilibriate_vertices(foam_mesh)
        fp.update_vertex_positions(foam_mesh)
        write_mesh(foam_mesh, "saved_foam_mesh.obj")
        fp.add_pressure(foam_mesh, equal)



    #-------------------------------------------------------------------------------------------------------------------
    # Do parameter test (to be moved to separate file)
    #-------------------------------------------------------------------------------------------------------------------
    if param: comprehensive_parameter_test(foam_mesh)


    #-------------------------------------------------------------------------------------------------------------------
    # Solve numerical model
    #-------------------------------------------------------------------------------------------------------------------

    plt.close()

    # QUASI-GLOBAL AND LOCAL

    #fs.solve_local_model(foam_mesh, max_iterations = 1700, method = 2, merit_tau = False, pressure_colour = False, colour_interval=100, save=True, save_interval=1)
    #fs.solve_local_model(max_iterations=700, method=2, merit_tau=False, pressure_colour=True, colour_interval=100, save=True, load=True, iterations_start=1018, save_interval=1)


    # # GLOBAL
    #
    # foam_mesh = TriMesh()
    # read_mesh(foam_mesh, "saved_foam_mesh.obj")
    # fp.add_pressure(foam_mesh, equal)
    # fp.add_mid_point(foam_mesh)
    # fp.update_vertex_positions(foam_mesh)
    # fp.add_face_connectivity(foam_mesh)
    #
    #fs.solve_global_model(foam_mesh, max_iterations=50, pressure_colour=True, colour_interval=10, load=True, save=True, save_interval=1, iterations_start=100)


    #-------------------------------------------------------------------------------------------------------------------
    # MERIT FUNCTION
    #-------------------------------------------------------------------------------------------------------------------
    # for fh in foam_mesh.faces():
    #     if not foam_mesh.is_boundary(fh):
    #         is_boundary = np.sum(np.array([1 for vh in foam_mesh.fv(fh) if foam_mesh.is_boundary(vh)]))
    #
    #         if not is_boundary:
    #             hs = np.linspace(0.00001,1,500)
    #             T_0 = []
    #             T_1 = []
    #             T_2 = []
    #             T_3 = []
    #             T_4 = []
    #             A_0 = []
    #             A_1 = []
    #             A_2 = []
    #             A_3 = []
    #             A_4 = []
    #
    #             for hx in hs:
    #                 j = fp.calc_jacobian(foam_mesh, fh, h = hx)
    #                 T_0.append(j[0][0])
    #                 T_1.append(j[0][1])
    #                 T_2.append(j[0][2])
    #                 T_3.append(j[0][3])
    #                 T_4.append(j[0][4])
    #                 A_0.append(j[2][0])
    #                 A_1.append(j[2][1])
    #                 A_2.append(j[2][2])
    #                 A_3.append(j[2][3])
    #                 A_4.append(j[2][4])
    #             break
    #
    # plt.figure()
    # plt.plot(hs, T_0, label='$/x$')
    # plt.plot(hs, T_1, label='$/y$')
    # plt.plot(hs, T_2, label='$/p_i$')
    # plt.plot(hs, T_3, label='$/p_j$')
    # plt.plot(hs, T_4, label='$/p_k$')
    # plt.xlabel('$h$', fontsize=18)
    # plt.ylabel('$\delta \Theta / \delta h$', fontsize=18)
    # plt.legend()
    # plt.xscale('log')
    # plt.savefig('images/hT_.eps', format='eps', dpi=1000)
    # plt.show()
    #
    # plt.figure()
    # plt.plot(hs, A_0, label='$/x$')
    # plt.plot(hs, A_1, label='$/y$')
    # plt.plot(hs, A_2, label='$/p_i$')
    # plt.plot(hs, A_3, label='$/p_j$')
    # plt.plot(hs, A_4, label='$/p_k$')
    # plt.xlabel('$h$', fontsize=18)
    # plt.ylabel('$\delta A / \delta h$', fontsize=18)
    # plt.legend()
    # plt.xscale('log')
    # plt.savefig('images/hA_.eps', format='eps', dpi=1000)
    # plt.show()


    #-------------------------------------------------------------------------------------------------------------------
    # PRESSURE
    #-------------------------------------------------------------------------------------------------------------------
    # # LOCAL
    # mesh = fp.load_foam_mesh(method_name='local', iteration=30)
    # pressure_list_local = []
    # for vh in mesh.vertices():
    #     if not mesh.is_boundary(vh):
    #         pressure_list_local.append(mesh.property(fp.pressure_handle, vh))
    #
    #
    ## QUASI
    # mesh = fp.load_foam_mesh(method_name='quasi_global', iteration=30)
    # pressure_list_quasi = []
    # for vh in mesh.vertices():
    #     if not mesh.is_boundary(vh):
    #         pressure_list_quasi.append(mesh.property(fp.pressure_handle, vh))
    #
    #
    # pressure_list_np = np.array(pressure_list_quasi)
    # np.save('saved/experiments/pressure_list_large_quasi.npy', pressure_list_np)
    #
    #
    # # GLOBAL
    # mesh = fp.load_foam_mesh(method_name='global', iteration=30)
    # pressure_list_global = []
    # for vh in mesh.vertices():
    #     if not mesh.is_boundary(vh):
    #         pressure_list_global.append(mesh.property(fp.pressure_handle, vh))
    #         print mesh.property(fp.pressure_handle, vh)
    #
    #
    # np.save('saved/experiments/pressure_list_local.npy', pressure_list_local)
    # np.save('saved/experiments/pressure_list_quasi.npy', pressure_list_quasi)
    # np.save('saved/experiments/pressure_list_global.npy', pressure_list_global)


    #plt.xscale('log')
    #plt.yscale('log')
    #plt.xlabel('Iteration', fontsize=18)
    #plt.ylabel('$1/2F^TF$', fontsize=18)
    #plt.savefig('../thesis/images/total_merit_nobacktrack.eps', format='eps', dpi=1000)


    #-------------------------------------------------------------------------------------------------------------------
    # NUMBER OF CELLS
    #-------------------------------------------------------------------------------------------------------------------


    # ## QUASI-GLOBAL -- LARGE
    # cell_count_list = []
    # for i in range(1700):
    #     print i
    #     mesh = fp.load_foam_mesh(method_name='quasi_global', iteration=i)
    #     cell_counter = 0
    #     for vh in mesh.vertices():
    #         if not mesh.is_boundary(vh):
    #             cell_counter += 1
    #
    #     cell_count_list.append(cell_counter)
    #
    #
    # cell_count_list_np = np.array(cell_count_list)
    # np.save('saved/experiments/cell_counter_list_quasi_large.npy', cell_count_list_np)


    # print cell_count_list[10], cell_count_list[100], cell_count_list[500], cell_count_list[1000]


    #fs.second_moment(2, 10)

    # mesh = fp.load_foam_mesh(method_name='quasi_global', iteration=1700)
    # cell_counter = 0
    # for vh in mesh.vertices():
    #     if not mesh.is_boundary(vh):
    #         cell_counter += 1
    #
    # print cell_counter




    #-------------------------------------------------------------------------------------------------------------------
    # SURFACE ENERGY
    #-------------------------------------------------------------------------------------------------------------------

    # ex1 = np.load('saved/experiments/quasi_global_1699_energy.npy')
    # plt.figure()
    # plt.plot(ex1, label='Local')
    # plt.xlabel('Iteration', fontsize=18)
    # plt.ylabel('Total film length', fontsize=18)
    # plt.savefig('images/film_length_local.eps', format='eps', dpi=1000)
    # plt.show()
    #
    # plt.figure()
    # plt.plot(ex1, label='Quasi-global')
    # plt.xlabel('Iteration', fontsize=18)
    # plt.ylabel('Total film length', fontsize=18)
    # plt.savefig('images/film_length_quasi.eps', format='eps', dpi=1000)
    # plt.show()


    # pressures = np.load('saved/saved_mesh_' + 'quasi_global' + '_' + str(50) + '_iterations_p.npy')
    # plt.figure(1)
    # plt.hist(pressures[1,:],bins=15)
    # plt.xlabel('Pressure')
    # plt.ylabel('Counts')
    # plt.title('Pressure histogram, quasi-global, after 50 iterations')
    # plt.show()
    #
    # pressures = np.load('saved/saved_mesh_' + 'quasi_global' + '_' + str(150) + '_iterations_p.npy')
    # plt.figure(2)
    # plt.hist(pressures[1,:],bins=15)
    # plt.xlabel('Pressure')
    # plt.ylabel('Counts')
    # plt.title('Pressure histogram, quasi-global, after 150 iterations')
    # plt.show()


    # mesh = fp.load_foam_mesh(method_name='local', iteration=869)
    # fv.draw_primary_mesh(mesh, title='local' + '/simulation_iteration_dual' + str(869).zfill(3), pressure_colour=False, draw_dual = True)
    # mesh = fp.load_foam_mesh(method_name='local', iteration=870)
    # fv.draw_primary_mesh(mesh, title='local' + '/simulation_iteration_dual' + str(870).zfill(3), pressure_colour=False, draw_dual = True)

    # #
    # mesh = fp.load_foam_mesh(method_name='quasi_global', iteration=1000)
    # fv.draw_primary_mesh(mesh, title='quasi_global' + '/simulation_iteration_test_2_' + str(1000).zfill(3), pressure_colour=False, draw_dual = True)

    #fs.draw_film_lengths(2, num_meshes=1700)


    #-------------------------------------------------------------------------------------------------------------------
    # Experiments
    #-------------------------------------------------------------------------------------------------------------------
    #fs.draw_film_lengths(1, num_meshes=1700)
    #fs.draw_film_lengths(2, num_meshes=1700)
    #fs.draw_mean_cell_area(1, num_meshes=1700)
    #fs.draw_mean_cell_area(2, num_meshes=1700)
    #fs.draw_mean_cell_area(2, num_meshes=1000)

    #fs.draw_convergence_plot(1, num_meshes=1700)
    #fs.draw_convergence_plot(2, num_meshes=1700)
    #fs.draw_convergence_plot(3, num_meshes=100)

    #fs.aboav_weaire(2, 1000)


    #fs.second_moment(2, 1000)



