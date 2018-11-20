from foam_import_openmesh import *
import foam_properties as fp
import foam_equilibriate as fe
import config
import numpy as np

# ----------------------------------------------------------------------------------------------------------------------
# Import constants
# ----------------------------------------------------------------------------------------------------------------------
wp_imp       = config.initialise['world_pressure']

def T1_process(mesh, eh):
    """ Do a T1-process, i.e. perform an edge-flip on a specified edge
        TriMesh mesh: an existing OpenMesh mesh
        openmesh.EdgeHandle edge_handle: a handle to an edge in the trimesh
    """

    # Fetch original junction coordinates
    fh_0 = mesh.face_handle(mesh.halfedge_handle(eh, 0))
    fh_1 = mesh.face_handle(mesh.halfedge_handle(eh, 1))
    x_I  = mesh.property(fp.mid_point_handle, fh_0)
    x_J  = mesh.property(fp.mid_point_handle, fh_1)

    # Create new point
    x_JI    = x_I - x_J
    x_perp  = TriMesh.Point(-x_JI[1], x_JI[0], 0)                                                                       # Perpendicular vector
    x_I_new = (x_I + x_J) / 2.0 + x_perp / 2.0
    x_J_new = (x_I + x_J) / 2.0 - x_perp / 2.0

    # Flip edge
    mesh.flip(eh)
    mesh.set_property(fp.mid_point_handle, fh_0, x_I_new)
    mesh.set_property(fp.mid_point_handle, fh_1, x_J_new)

    # Update "is flipped"-handle
    mesh.set_property(fp.is_flipped_handle, fh_0, True)
    mesh.set_property(fp.is_flipped_handle, fh_1, True)


def T2_process(mesh, vh, cell):
    """ Do a T2-process, i.e. perform a face collapse of a face in the primary mesh
        TriMesh mesh: an existing OpenMesh mesh
        openmesh.VertexHandle face_handle: a handle to an vertex in the dual trimesh
    """

    # Request status. Otherwise, a breakdown will occur, when trying to delete a mesh element
    mesh.request_face_status()
    mesh.request_vertex_status()
    mesh.request_halfedge_status()
    mesh.request_edge_status()

    if cell == 3:
        # Identify vertices that will eventually span the resulting face of the T2-operation
        vertices = np.array([[mesh.to_vertex_handle(heh), mesh.to_vertex_handle(heh).idx()] for heh in mesh.voh(vh)])
        mid_point_new = TriMesh.Point(0, 0, 0)
        for fh0 in mesh.vf(vh):
            if mesh.is_boundary(fh0) == False:
                    mid_point_new += mesh.property(fp.mid_point_handle, fh0) / 3.0

        # Check, if were dealing with a boundary cell with faces facing the world (some boundary cells are built around non-boundary faces
        boundary_cell = True if np.sum(np.array([1 for fh0 in mesh.vf(vh) if mesh.is_boundary(fh0)])) > 0 else False

        if boundary_cell == False:
            # Find edge to collapse, ultimately resulting in a face collapse at the same time
            for heh in mesh.voh(vh):
                if mesh.is_collapse_ok(heh):
                    mesh.collapse(heh)
                    break

            # Find resulting face
            for fh0 in mesh.vf(vertices[0][0]):
                indices = np.array([vh0.idx() for vh0 in mesh.fv(fh0)])
                if np.sum(np.in1d(vertices[:, 1], indices)) == 3:
                    fh_new = fh0
                    break

        else:
            # Find resulting face
            for fh0 in mesh.vf(vh):
                half_boundary_face = np.sum(np.array([1 for vh0 in mesh.fv(fh0) if mesh.is_boundary(vh0)]))
                if half_boundary_face == 0:                                                                                 # Identify resulting face
                    fh_new = fh0
                    break

            # delete_list = []
            # for fh0 in mesh.vf(vh):
            #     if mesh.is_boundary(fh0):
            #         delete_list.append(fh0)
            delete_list = [fh0 for fh0 in mesh.vf(vh) if mesh.is_boundary(fh0)]
            for fh0 in delete_list: mesh.delete_face(fh0, True)

            mesh.set_property(fp.pressure_handle, vh, wp_imp)

        mesh.set_property(fp.mid_point_handle, fh_new, mid_point_new)

    else:                                                                                                               # If cell == 2
        delete_list = [fh0 for fh0 in mesh.vf(vh) if mesh.is_boundary(fh0)]
        for fh0 in delete_list: mesh.delete_face(fh0, True)
        mesh.set_property(fp.pressure_handle, vh, wp_imp)

    mesh.garbage_collection()



def perform_topology_changes(mesh):
    # T2
    # cell_areas = np.array([mesh.property(fp.area_handle, vh0) for vh0 in mesh.vertices() if mesh.is_boundary(vh0) == False])
    # min_area   = np.mean(cell_areas) - 0.5 * np.std(cell_areas)
    #print np.mean(cell_areas), min_area
    boundary_T12_performed = False                                                                                      # Initialise flag to check, if T1 or T2 is performed on boundary
    min_area   = 0.0639229521832

    for vh in mesh.vertices():
        if mesh.is_boundary(vh) == False:
            neighbours    = np.sum(np.array([1 for vh0 in mesh.vv(vh) if mesh.is_boundary(vh0)==False]))
            boundary_cell = True if np.sum(np.array([1 for vh0 in mesh.vv(vh) if mesh.is_boundary(vh0)])) > 0 else False
            if boundary_cell: neighbours += 1

            if neighbours == 3:
                if fp.calc_areas_bubble(mesh, vh, calc_target = False) < min_area:
                    T2_process(mesh, vh, cell=3)

                    if boundary_cell: boundary_T12_performed = True
            elif neighbours == 2:
                if fp.calc_areas_bubble(mesh, vh, calc_target=False) < min_area:
                    T2_process(mesh, vh, cell=2)

                    if boundary_cell: boundary_T12_performed = True

    mesh.garbage_collection()

    # T1
    # edge_lengths = np.array([mesh.property(fp.film_length_handle, eh) for eh in mesh.edges() if mesh.is_boundary(eh)==False])
    # edge_lengths = edge_lengths[edge_lengths != np.array(None)]                                                         # Remove 'None's from list, which result from redundant edge's missing film length handles
    # min_length   = np.mean(edge_lengths)-2 * np.std(edge_lengths)
    # print np.mean(edge_lengths), min_length

    fp.update_vertex_positions(mesh)
    fp.prop_is_flipped(mesh)

    min_length   = 0.0638045508824

    for eh in mesh.edges():
        if mesh.is_boundary(eh) == False:
            heh_0 = mesh.halfedge_handle(eh, 0)
            heh_1 = mesh.halfedge_handle(eh, 1)
            vh_0 = mesh.to_vertex_handle(heh_0)
            vh_1 = mesh.to_vertex_handle(heh_1)
            fh_0 = mesh.face_handle(mesh.halfedge_handle(eh, 0))  # Face handle for the one neighboring face
            fh_1 = mesh.face_handle(mesh.halfedge_handle(eh, 1))  # Face handle for the other neighboring face
            boundary_cell = True if mesh.is_boundary(fh_0) or mesh.is_boundary(fh_1) else False

            if not mesh.property(fp.is_flipped_handle, fh_0) and not mesh.property(fp.is_flipped_handle, fh_1) and not boundary_cell and mesh.is_flip_ok(eh):
                coor_0  = np.array([mesh.property(fp.mid_point_handle, fh_0)[0], mesh.property(fp.mid_point_handle, fh_0)[1]])
                coor_1  = np.array([mesh.property(fp.mid_point_handle, fh_1)[0], mesh.property(fp.mid_point_handle, fh_1)[1]])
                delta_p = mesh.property(fp.pressure_handle, vh_1) - mesh.property(fp.pressure_handle, vh_0)
                film_length = fp.calc_film_length(coor_0, coor_1, delta_p)

                if film_length < min_length:
                    T1_process(mesh, eh)
                    if boundary_cell: boundary_T12_performed = True

    # Reduce boundary facing vertices
    if boundary_T12_performed: fe.reduce_boundary_facing_faces(mesh)

    # Update properties
    fp.update_vertex_positions(mesh)
    fp.add_face_connectivity(mesh)

    mesh.garbage_collection()


