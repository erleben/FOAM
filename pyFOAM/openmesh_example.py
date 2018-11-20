# -------------------------------------------------------------------------------
# Imports
# -------------------------------------------------------------------------------
import math
import sys
import numpy as np
import random
import collections
import copy
import svgwrite
import matplotlib.tri as tri
import svgutils.transform as st
import matplotlib.pyplot as plt
from reportlab.graphics import renderPDF
from svglib.svglib import svg2rlg
from operator import add


# -------------------------------------------------------------------------------
# Additional import management
# -------------------------------------------------------------------------------
## Manage paths
sys.path.append("/Users/jtbondorf/bin/OpenMesh/build/Build/python")             # Add OpenMesh library to Python's search path
sys.path.append("/Users/jtbondorf/Documents/Universitet/Speciale/2016_JOHAN_TEICHERT_BONDORF/src")

## Import other files
import config
from parameter_test import *

## Import OpenMesh
from openmesh import *


# ----------------------------------------------------------------------------------------------------------------------
# Create mesh
# ----------------------------------------------------------------------------------------------------------------------

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


# ----------------------------------------------------------------------------------------------------------------------
# Add or update properties
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


def add_pressure(mesh, equal=True, approx_pressure=1, world_pressure_fraction=0.9):
    """ Add pressure property to an existing mesh
        TriMesh mesh: an existing OpenMesh mesh
        bool equal: True if equal pressure value is wanted, False if different pressure values are wanted
        double approx_pressure: value for pressure, if equal pressure across mesh is wanted. If not, it is used as offset for randomisation
        double world_pressure_fraction: fraction, for how much lower the world pressure is than the approximate pressure inside the foam
    """


    # Define lower world pressure
    global world_pressure
    if equal:
        world_pressure = approx_pressure
    else:
        world_pressure = world_pressure_fraction * approx_pressure

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


def update_boundary_pressure(mesh):
    """ Update boundaries existing mesh
        TriMesh mesh: an existing OpenMesh mesh
        bool equal: True if equal pressure value is wanted, False if different pressure values are wanted
        double pressure: value for pressure, if equal pressure across mesh is wanted
    """

    for vh in mesh.vertices():
        if mesh.is_boundary(vh):
            mesh.set_property(pressure_handle, vh, world_pressure)


def add_boundary_neighbour_faces(mesh):
    """ Add neighbour face property to each boundary face existing mesh
        TriMesh mesh: an existing OpenMesh mesh
    """

    # Define and make property handle global, so that it is accessible from other functions
    global boundary_neighbour_faces_handle
    boundary_neighbour_faces_handle = FPropHandle()
    mesh.add_property(boundary_neighbour_faces_handle, "fprop_float")

    for fh0 in mesh.faces():
        if mesh.is_boundary(fh0):
            for eh in mesh.fe(fh0):
                right_or_left = mesh.is_boundary(eh)
                break

            neighbour_faces = [mesh.property(mid_point_handle, fh) for fh in mesh.ff(fh0)]
            if right_or_left == False: neighbour_faces = [neighbour_faces[1], neighbour_faces[0]]

            mesh.set_property(boundary_neighbour_faces_handle, fh0,
                              neighbour_faces)  # For boundary faces, add face handles for neighbour faces
        else:
            mesh.set_property(boundary_neighbour_faces_handle, fh0, [])


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
            is_it_boundary = sum([mesh.is_boundary(vh0) for vh0 in mesh.fv(fh)])  # Check, if face bleongs to junction on boundary of foam

            if is_it_boundary == 0:
                # Create list of edge neighboring faces to fh
                neighbour_faces = np.array([fh0 for fh0 in mesh.ff(fh)])
                neighbour_faces = np.append(neighbour_faces, neighbour_faces[0])
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
                                    fh_start = fh_next  # Move one step forward
                                    fh_next = fh_iter  # Move one step forward
                                    break

                        neighbour_faces = np.append(neighbour_faces, fh_next)

                neighbour_faces = np.append(neighbour_faces, neighbour_faces[0])

            # Properties of current face
            vertices_i = np.array([[vh0.idx(), vh0, mesh.property(pressure_handle, vh0), mesh.is_boundary(vh0)] for vh0 in mesh.fv(fh)])  # List of vertices and pressures in fh
            x_i = np.array([mesh.property(mid_point_handle, fh)[0], mesh.property(mid_point_handle, fh)[1]])  # Junction coordinate using notation from thesis

            # Initialise list containing of 3 angles
            variables = np.empty([3, 8])

            for i in range(len(neighbour_faces) - 1):
                # Initialise faces
                f_j = neighbour_faces[i]
                f_k = neighbour_faces[i + 1]

                # j
                x_j        = np.array([mesh.property(mid_point_handle, f_j)[0], mesh.property(mid_point_handle, f_j)[1]])
                x_ij       = x_j - x_i
                vertices_j = np.array([[vh0.idx(), vh0, mesh.is_boundary(vh0)] for vh0 in mesh.fv(f_j)])

                # k
                x_k        = np.array([mesh.property(mid_point_handle, f_k)[0], mesh.property(mid_point_handle, f_k)[1]])
                x_ik       = x_k - x_i
                vertices_k = np.array([[vh0.idx(), vh0, mesh.is_boundary(vh0)] for vh0 in mesh.fv(f_k)])

                # Indices for vertices in fh in correct order
                if is_it_boundary == 0:
                    idx_pi = np.in1d(vertices_i[:, 0], vertices_j[:, 0]) & np.in1d(vertices_i[:, 0], vertices_k[:, 0])  # Common vertex of j and k
                    idx_pj = np.in1d(vertices_i[:, 0], vertices_j[:, 0]) & np.invert(np.in1d(vertices_i[:, 0], vertices_k[:, 0]))
                    idx_pk = np.in1d(vertices_i[:, 0], vertices_k[:, 0]) & np.invert(np.in1d(vertices_i[:, 0], vertices_j[:, 0]))
                elif is_it_boundary:
                    is_boundary_j = np.any(vertices_j[:, 2])
                    is_boundary_k = np.any(vertices_k[:, 2])

                    # Different boundary scenarios
                    if is_boundary_j and is_boundary_k:
                        if f_j.idx() == f_k.idx():  # Treating special cae with isolated bubble on boundary
                            idx_temp = np.in1d(vertices_i[:, 0],vertices_j[:, 0])  # = np.in1d(v_i,v_k), since f_j=f_k
                            vertices_temp = vertices_i[idx_temp][:, 1]
                            neigh_ver = np.array([not (sum([not (mesh.is_boundary(vh)) for vh in mesh.vv(vh0)]) > 1) for vh0 in vertices_temp])  # List of: "is vertex a 2-bubble vertex?"
                            idx_i = vertices_temp[neigh_ver][0].idx()  # Vertex index of p_i in vertices_i
                            idx_j = vertices_temp[np.invert(neigh_ver)][0].idx()
                            idx_pi = np.in1d(vertices_i[:, 0], idx_i)
                            idx_pj = np.in1d(vertices_i[:, 0], idx_j)
                            idx_pk = np.invert(idx_pi | idx_pj)
                        else:
                            idx_temp = vertices_i[:, 3]
                            idx_pi = np.array([x for x in idx_temp])  # Weird work around. Otherwise, np.invert does not work
                            idx_pj = np.in1d(vertices_i[:, 0], vertices_j[:, 0]) & np.invert(idx_pi)
                            idx_pk = np.in1d(vertices_i[:, 0], vertices_k[:, 0]) & np.invert(idx_pi)
                    elif is_boundary_j and not (is_boundary_k):
                        idx_pi = np.in1d(vertices_i[:, 0], vertices_j[:, 0]) & np.in1d(vertices_i[:, 0], vertices_k[:, 0])
                        idx_pk = np.in1d(vertices_i[:, 0], vertices_k[:, 0]) & np.invert(np.in1d(vertices_i[:, 0], vertices_j[:, 0]))
                        idx_pj = np.invert(idx_pi | idx_pk)
                    elif is_boundary_k and not (is_boundary_j):
                        idx_pi = np.in1d(vertices_i[:, 0], vertices_j[:, 0]) & np.in1d(vertices_i[:, 0],vertices_k[:, 0])
                        idx_pj = np.in1d(vertices_i[:, 0], vertices_j[:, 0]) & np.invert(np.in1d(vertices_i[:, 0], vertices_k[:, 0]))
                        idx_pk = np.invert(idx_pi | idx_pj)

                # Associated vertex index
                idx_i    = vertices_i[idx_pi, 0][0]

                # Pressure values
                p_i = vertices_i[idx_pi, 2][0]
                p_j = vertices_i[idx_pj, 2][0]
                p_k = vertices_i[idx_pk, 2][0]

                # Set row in variable matrix
                variables[i, :] = [idx_i, x_ij[0], x_ij[1], x_ik[0], x_ik[1], p_i, p_j, p_k]

            # Assign as property
            mesh.set_property(face_connectivity_handle, fh, variables)


def calc_film_lengths(mesh, gamma=73.0):
    """ Calculate film lengths of all films in mesh. Add as property to each film-belonging edge in dual mesh.
        TriMesh mesh: an existing OpenMesh mesh
    """

    # Define and make property handle global, so that it is accessible from other functions
    global film_length_handle                                               # Define and make property handle global, so that it is accessible from other functions
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
                    fh_next = f0_h if mesh.is_boundary(f0_h) else f1_h
                    fh_start = f0_h if not mesh.is_boundary(f0_h) else f1_h
                    coor_flm_end_0 = np.array(
                        [mesh.property(mid_point_handle, fh_start)[0], mesh.property(mid_point_handle, fh_start)[1]])

                    # Visit neighbouring faces to find the other coordinate
                    while mesh.is_boundary(fh_next):
                        for fh in mesh.ff(fh_next):
                            if fh != fh_start:
                                fh_start = fh_next  # Move one step forward
                                fh_next = fh  # Move one step forward
                                break
                    coor_flm_end_1 = np.array(
                        [mesh.property(mid_point_handle, fh_next)[0], mesh.property(mid_point_handle, fh_next)[1]])

                # Handle normal case
                else:
                    # Calculate the two edge vectors spanning the part of the bubble film (using the one or the other bubble centrum does not matter
                    coor_flm_end_0 = np.array([mesh.property(mid_point_handle, f0_h)[0], mesh.property(mid_point_handle, f0_h)[1]])
                    coor_flm_end_1 = np.array([mesh.property(mid_point_handle, f1_h)[0], mesh.property(mid_point_handle, f1_h)[1]])

                d              = np.sqrt(np.sum((coor_flm_end_1-coor_flm_end_0)**2))    # Calculate distance between two end points of bubble film
                delta_p        = mesh.property(pressure_handle, v1_h) - mesh.property(pressure_handle, v0_h)

                if delta_p != 0:
                    r              = 2.0 * gamma / abs(delta_p)                          # Curvature radius
                    theta          = np.arcsin(d/r)
                    film_length    = 2*theta*r
                else:
                    film_length    = d

                mesh.set_property(film_length_handle, eh, film_length)

            else:
                pass


def calc_plateau_angles_original(mesh, gamma=73.0):
    """ Calculate angles of plateau junctions in mesh. Add as property to each non-boundary face in dual mesh.
        TriMesh mesh: an existing OpenMesh mesh
    """

    # HANDLES
    global plateau_angle_handle                                                        # Define and make property handle global, so that it is accessible from other functions
    plateau_angle_handle = FPropHandle()
    mesh.add_property(plateau_angle_handle, "fprop_float")

    for fh in mesh.faces():
        if not mesh.is_boundary(fh):
            is_it_boundary = sum([mesh.is_boundary(vh0) for vh0 in mesh.fv(fh)])       # Check, if face bleongs to junction on boundary of foam

            if is_it_boundary==0:
                # Create list of edge neighboring faces to fh
                neighbour_faces = np.array([fh0 for fh0 in mesh.ff(fh)])
                neighbour_faces = np.append(neighbour_faces,neighbour_faces[0])
            elif is_it_boundary==1:
                neighbour_faces = np.array([])
                for fh0 in mesh.ff(fh):
                    if not mesh.is_boundary(fh0):
                        neighbour_faces = np.append(neighbour_faces,fh0)
                    else:
                        # Visit neighbouring faces to find the other coordinate
                        fh_next  = fh0
                        fh_start = fh

                        while mesh.is_boundary(fh_next):
                            for fh_iter in mesh.ff(fh_next):
                                if fh_iter.idx() != fh_start.idx():
                                    fh_start = fh_next                                  # Move one step forward
                                    fh_next  = fh_iter                                  # Move one step forward
                                    break

                        neighbour_faces = np.append(neighbour_faces,fh_next)

                neighbour_faces = np.append(neighbour_faces,neighbour_faces[0])

            # Properties of current face
            vertices_i = np.array([[vh0.idx(), vh0, mesh.property(pressure_handle,vh0), mesh.is_boundary(vh0)] for vh0 in mesh.fv(fh)])   # List of vertices and pressures in fh
            x_i        = np.array([mesh.property(mid_point_handle, fh)[0],mesh.property(mid_point_handle, fh)[1]]) # Junction coordinate using notation from thesis

            # Initialise list containing of 3 angles
            thetas = np.empty([3,2])

            for i in range(len(neighbour_faces)-1):
                # Initialise faces
                f_j = neighbour_faces[i]
                f_k = neighbour_faces[i+1]

                # j
                x_j        = np.array([mesh.property(mid_point_handle, f_j)[0],mesh.property(mid_point_handle, f_j)[1]])
                x_ij       = x_j - x_i
                d_j        = np.sqrt(np.sum((x_ij) ** 2))
                vertices_j = np.array([[vh0.idx(), vh0, mesh.is_boundary(vh0)] for vh0 in mesh.fv(f_j)])

                # k
                x_k        = np.array([mesh.property(mid_point_handle, f_k)[0],mesh.property(mid_point_handle, f_k)[1]])
                x_ik       = x_k - x_i
                d_k        = np.sqrt(np.sum((x_ik) ** 2))
                vertices_k = np.array([[vh0.idx(), vh0, mesh.is_boundary(vh0)] for vh0 in mesh.fv(f_k)])


                # Indices for vertices in fh in correct order
                if is_it_boundary==0:
                    idx_pi = np.in1d(vertices_i[:,0], vertices_j[:,0]) & np.in1d(vertices_i[:,0], vertices_k[:,0])    # Common vertex of j and k
                    idx_pj = np.in1d(vertices_i[:,0], vertices_j[:,0]) & np.invert(np.in1d(vertices_i[:,0], vertices_k[:,0]))
                    idx_pk = np.in1d(vertices_i[:,0], vertices_k[:,0]) & np.invert(np.in1d(vertices_i[:,0], vertices_j[:,0]))
                elif is_it_boundary:
                    is_boundary_j = np.any(vertices_j[:,2])
                    is_boundary_k = np.any(vertices_k[:,2])

                    # Different boundary scenarios
                    if is_boundary_j and is_boundary_k:
                        if f_j.idx() == f_k.idx():                                                 # Treating special case with isolated bubble on boundary
                            ex = 0
                            idx_temp      = np.in1d(vertices_i[:,0],vertices_j[:,0])               # = np.in1d(v_i,v_k), since f_j=f_k
                            vertices_temp = vertices_i[idx_temp][:,1]
                            neigh_ver     = np.array([not(sum([not(mesh.is_boundary(vh)) for vh in mesh.vv(vh0)])>1) for vh0 in vertices_temp]) # List of: "is vertex a 2-bubble vertex?"
                            idx_i         = vertices_temp[neigh_ver][0].idx()                      # Vertex index of p_i in vertices_i
                            idx_j         = vertices_temp[np.invert(neigh_ver)][0].idx()
                            idx_pi        = np.in1d(vertices_i[:,0], idx_i)
                            idx_pj        = np.in1d(vertices_i[:,0], idx_j)
                            idx_pk        = np.invert(idx_pi | idx_pj)
                        else:
                            ex = 1
                            idx_temp = vertices_i[:,3]
                            idx_pi = np.array([x for x in idx_temp])                                   # Weird work around. Otherwise, np.invert does not work
                            idx_pj = np.in1d(vertices_i[:,0], vertices_j[:,0]) & np.invert(idx_pi)
                            idx_pk = np.in1d(vertices_i[:,0], vertices_k[:,0]) & np.invert(idx_pi)
                    elif is_boundary_j and not(is_boundary_k):
                        ex = 2
                        idx_pi = np.in1d(vertices_i[:,0], vertices_j[:,0]) & np.in1d(vertices_i[:,0], vertices_k[:,0])
                        idx_pk = np.in1d(vertices_i[:,0], vertices_k[:,0]) & np.invert(np.in1d(vertices_i[:,0], vertices_j[:,0]))
                        idx_pj = np.invert(idx_pi | idx_pk)
                    elif is_boundary_k and not(is_boundary_j):
                        ex = 3
                        idx_pi = np.in1d(vertices_i[:,0], vertices_j[:,0]) & np.in1d(vertices_i[:,0], vertices_k[:,0])
                        idx_pj = np.in1d(vertices_i[:,0], vertices_j[:,0]) & np.invert(np.in1d(vertices_i[:,0], vertices_k[:,0]))
                        idx_pk = np.invert(idx_pi | idx_pj)

                # Pressure values
                p_i = vertices_i[idx_pi, 2][0]
                p_j = vertices_i[idx_pj, 2][0]
                p_k = vertices_i[idx_pk, 2][0]

                # Pressure difference values
                delta_pj = p_j - p_i
                delta_pk = p_k - p_i

                # Curvature radii
                R_j = 2 * gamma / abs(delta_pj) if delta_pj != 0 else np.inf
                R_k = 2 * gamma / abs(delta_pk) if delta_pk != 0 else np.inf

                # Angle, notation used: see thesis
                # if np.linalg.det([x_ij,x_ik])>0:
                #     theta = np.arccos(np.dot(x_ij,x_ik)/(d_j*d_k)) + 0.5*(np.sign(delta_pj)*np.arcsin(d_j/(2*R_j)) + np.sign(delta_pk)*np.arcsin(d_k/(2*R_k)))
                # else:
                #     theta = (2 * np.pi - np.arccos(np.dot(x_ij, x_ik) / (d_j * d_k))) + 0.5 * (np.sign(delta_pj) * np.arcsin(d_j / (2 * R_j)) + np.sign(delta_pk) * np.arcsin(d_k / (2 * R_k)))

                theta = calc_angle(x_ij, x_ik, delta_pj, delta_pk)

                # Add to list
                idx_i       = vertices_i[idx_pi,0][0]
                thetas[i,:] = [theta,idx_i]

            # print is_it_boundary, thetas

            # Add list of angles as property to current face
            mesh.set_property(plateau_angle_handle, fh, thetas)


def calc_plateau_angles_new(mesh, gamma=73.0):
    """ Calculate angles of plateau junctions in mesh. Add as property to each non-boundary face in dual mesh.
        TriMesh mesh: an existing OpenMesh mesh
    """

    # HANDLES
    global plateau_angle_handle                                                        # Define and make property handle global, so that it is accessible from other functions
    plateau_angle_handle = FPropHandle()
    mesh.add_property(plateau_angle_handle, "fprop_float")

    # Add information about connectivity of faces
    add_face_connectivity(foam_mesh)

    for fh in mesh.faces():
        if not mesh.is_boundary(fh):
            thetas    = np.empty([3,2])
            variables = mesh.property(face_connectivity_handle, fh)                    # i, x_ij[0], x_ij[1], x_ik[0], x_ik[1], p_i, p_j, p_k x 3

            for i in range(3):
                # Fetch variables for one vertex
                [idx_i, x_ij_x, x_ij_y, x_ik_x, x_ik_y, p_i, p_j, p_k] = variables[i]
                x_ij = np.array([x_ij_x, x_ij_y])
                x_ik = np.array([x_ik_x, x_ik_y])

                # Calculate pressure differences
                delta_pj = p_j - p_i
                delta_pk = p_k - p_i

                # Calculate angle
                theta = calc_angle(x_ij, x_ik, delta_pj, delta_pk)

                # Add to list
                thetas[i,:] = [theta, idx_i]


            # Add list of angles as property to current face
            mesh.set_property(plateau_angle_handle, fh, thetas)


def calc_bubble_areas(mesh, gamma=73.0):
    """ Calculate areas of bubbles in mesh. Add as property to each bubble-belonging vertex in dual mesh.
        TriMesh mesh: an existing OpenMesh mesh
    """

    # HANDLES
    global bubble_area_handle  # Define and make property handle global, so that it is accessible from other functions
    bubble_area_handle = VPropHandle()
    mesh.add_property(bubble_area_handle, "vprop_float")

    # VARIABLES

    for vh in mesh.vertices():
        if not mesh.is_boundary(vh):
            # Find coordinate of this vertex
            coor_bubble = np.array([mesh.point(vh)[0], mesh.point(vh)[1]])
            idx_bubble  = vh.idx()
            p_bubble    = mesh.property(pressure_handle, vh)

            # Initialise area counter
            area_bubble = 0

            # Treat segments in bubble
            faces = []

            # Build faces list
            for fh in mesh.vf(vh):
                if mesh.is_boundary(fh) == False: faces.append(fh)
            faces.append(faces[0])

            for i in range(len(faces) - 1):
                coor_0 = np.array([mesh.property(mid_point_handle, faces[i])[0],mesh.property(mid_point_handle, faces[i])[1]])
                coor_1 = np.array([mesh.property(mid_point_handle, faces[i + 1])[0],mesh.property(mid_point_handle, faces[i + 1])[1]])

                # Find vertex belonging to two last faces and to adjacent bubble
                vertices_0 = np.array([[vh0.idx(), vh0] for vh0 in mesh.fv(faces[i]) if vh0.idx() != idx_bubble])
                vertices_1 = np.array([[vh0.idx(), vh0] for vh0 in mesh.fv(faces[i + 1]) if vh0.idx() != idx_bubble])
                idx        = np.in1d(vertices_0[:, 0], vertices_1[:, 0])

                if sum(idx) == 0:  # Handle case of boundary film, where two or more world bubble vertices are present
                    common_vh = [vh0 for vh0 in mesh.fv(faces[i]) if mesh.is_boundary(vh0)][0]
                else:
                    common_vh = vertices_0[idx, 1][0]
                delta_p = p_bubble - mesh.property(pressure_handle, common_vh)
                d       = np.sqrt(np.sum((coor_1 - coor_0) ** 2))                           # Euclidean distance between film ends

                area_tria = 0.5 * abs(np.cross(coor_0 - coor_bubble, coor_1 - coor_bubble))

                if delta_p != 0:
                    r = 2.0 * gamma / abs(delta_p)  # Curvature radius
                    theta = np.arcsin(d / r)
                    area_temp = theta * r ** 2 - abs(0.5 * r ** 2 * np.sin(2 * theta))
                    area_segm = area_temp if delta_p > 0 else -area_temp
                else:
                    area_segm = 0

                area_slice = area_tria + area_segm

                area_bubble += area_slice

            mesh.set_property(bubble_area_handle, vh, area_bubble)


def calc_jacobian(mesh, h=0.01):
    """ Calculate jacobian for junction in mesh. Add as property to each cell vertex in mesh.
        TriMesh mesh: an existing OpenMesh mesh
        double h: describing step distance in approximation of derivative
    """

    # HANDLES
    global jacobian_handle          # Define and make property handle global, so that it is accessible from other functions
    jacobian_handle = VPropHandle()
    mesh.add_property(jacobian_handle, "vprop_float")

    for fh in mesh.faces():
        if not mesh.is_boundary(fh):
            # VARIABLES
            jacobian = np.empty(5, 5)
            x_j  = np.array([mesh.property(mid_point_handle, f_j)[0], mesh.property(mid_point_handle, f_j)[1]])
            x_k  = np.array([mesh.property(mid_point_handle, f_k)[0], mesh.property(mid_point_handle, f_k)[1]])
            x_ij = x_j - x_i
            x_ik = x_k - x_i

            theta_i_x = 1



# ----------------------------------------------------------------------------------------------------------------------
# Manipulate dual mesh
# ----------------------------------------------------------------------------------------------------------------------

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

        coor_iter          = 0                                                          # Initialise coor_iter as being different from coor_original
        heh_iter           = heh_original                                               # Initialise halfedge iterator
        to_be_flipped_list = []                                                         # Initialise list of edges to be flipped


        # Do edge flip, if iterator comes across more world facing faces (characterised by having a vertex with only 3 connecting edges)
        while coor_iter != coor_original:
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
    add_mid_point(mesh)                                                             # Update mid point face property


# ----------------------------------------------------------------------------------------------------------------------
# Helper functions
# ----------------------------------------------------------------------------------------------------------------------

def calc_angle(x_ij, x_ik, delta_pj, delta_pk, gamma=73.0):
    """ Calculate an angle given three points and two pressure differences
        np.array x_ij, x_ik: 2-dimensional vectors from x_i to x_j and x_k, respectively
        double delta_pj, delta_pk: pressure differences between cell i with j and k, respectively
        return double theta: angle
    """
    d_j = np.sqrt(np.sum((x_ij) ** 2))
    d_k = np.sqrt(np.sum((x_ik) ** 2))
    R_j = 2 * gamma / abs(delta_pj) if delta_pj != 0 else np.inf
    R_k = 2 * gamma / abs(delta_pk) if delta_pk != 0 else np.inf

    if np.linalg.det([x_ij, x_ik]) > 0:
        theta = np.arccos(np.dot(x_ij, x_ik) / (d_j * d_k)) + 0.5 * (np.sign(delta_pj) * np.arcsin(d_j / (2 * R_j))
                + np.sign(delta_pk) * np.arcsin(d_k / (2 * R_k)))
    else:
        theta = (2 * np.pi - np.arccos(np.dot(x_ij, x_ik) / (d_j * d_k))) + 0.5 * (np.sign(delta_pj) * np.arcsin(d_j / (2 * R_j))
                + np.sign(delta_pk) * np.arcsin(d_k / (2 * R_k)))

    return theta


# ----------------------------------------------------------------------------------------------------------------------
# Equilibriate foam mesh
# ----------------------------------------------------------------------------------------------------------------------

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
    for iter in range(max_iterations):
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

    for iter in range(max_iterations):
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

    add_pressure(mesh)


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


# ----------------------------------------------------------------------------------------------------------------------
# Topological processes
# ----------------------------------------------------------------------------------------------------------------------

def T1_process(mesh,edge_handle):
    """ Do a T1-process, i.e. perform an edge-flip on a specified edge
        TriMesh mesh: an existing OpenMesh mesh
        openmesh.EdgeHandle edge_handle: a handle to an edge in the trimesh
    """
    mesh.flip(edge_handle)


def T2_process(mesh,vertex_handle):
    """ Do a T2-process, i.e. perform a face collapse of a face in the primary mesh
        TriMesh mesh: an existing OpenMesh mesh
        openmesh.VertexHandle face_handle: a handle to an vertex in the dual trimesh
    """
    mesh.request_face_status()
    mesh.request_vertex_status()
    mesh.request_halfedge_status()
    mesh.request_edge_status()

    for heh in mesh.voh(vertex_handle):
        if mesh.is_collapse_ok(heh):
            mesh.collapse(heh)
            break

    mesh.garbage_collection()


# ----------------------------------------------------------------------------------------------------------------------
# Visualise mesh
# ----------------------------------------------------------------------------------------------------------------------

def draw_dual_mesh(mesh,title='dual_mesh',style='edges',scale=50, offset=10):
    """ Draw dual mesh and write svg
        TriMesh mesh: an existing dual OpenMesh mesh
        string style: 'edges' or 'edges+vertices' for mesh drawing style
        string title: title of svg-file to be saved to
        int scale: scale defining scale factor of visualization
        int offset: added to move drawing away from screen edges
    """
    dwg = svgwrite.Drawing('images/'+title+'.svg')

    # Draw edges
    if style == 'edges':
        for eh in mesh.edges():
            h0_h = mesh.halfedge_handle(eh,0)                                           # Fetch handle for halfedge in edge
            h1_h = mesh.halfedge_handle(eh,1)
            v0_h = mesh.to_vertex_handle(h0_h)                                          # Converts end point of halfedge handle to vertex handle
            v1_h = mesh.to_vertex_handle(h1_h)

            v0   = (mesh.point(v0_h)[0]*scale+offset,mesh.point(v0_h)[1]*scale+offset)  # Convert to tuple for coordinate
            v1   = (mesh.point(v1_h)[0]*scale+offset,mesh.point(v1_h)[1]*scale+offset)  # Convert to tuple for other coordinate

            dwg.add(dwg.line(v0, v1, stroke=svgwrite.rgb(80, 80, 80, '%')))             # Add to svg-file

    # Add vertices, if specified
    if style != 'edges':
        for vh in mesh.vertices():
            coor = (mesh.point(vh)[0]*scale+offset,mesh.point(vh)[1]*scale+offset)      # Fetch coordinate for point
            dwg.add(dwg.circle(center=coor, r=2))


    dwg.save()


def draw_primary_mesh(mesh,title='primary_mesh',draw_boundary=True,scale=50, offset=10, gamma=73.0):
    """ Draw dual mesh and write svg
        TriMesh mesh: an existing dual OpenMesh mesh
        string title: title of svg-file to be saved to
        bool draw_boundary: specifies, if boundary has to be drawn or not
        int scale: scale defining scale factor of visualization
        int offset: added to move drawing away from screen edges
        double gamma: surface tension
    """

    ## Draw dual mesh
    dwg  = svgwrite.Drawing('images/'+title+'.svg')

    # DRAW DUAL MESH
    for eh in mesh.edges():

        h0_h = mesh.halfedge_handle(eh,0)                                           # Fetch handle for halfedge in edge
        h1_h = mesh.halfedge_handle(eh,1)
        v0_h = mesh.to_vertex_handle(h0_h)                                          # Converts end point of halfedge handle to vertex handle
        v1_h = mesh.to_vertex_handle(h1_h)

        v0   = (mesh.point(v0_h)[0]*scale+offset,mesh.point(v0_h)[1]*scale+offset)  # Convert to tuple for coordinate
        v1   = (mesh.point(v1_h)[0]*scale+offset,mesh.point(v1_h)[1]*scale+offset)  # Convert to tuple for other coordinate

        dwg.add(dwg.line(v0, v1, stroke=svgwrite.rgb(80, 80, 80, '%')))             # Add to svg-file


    # DRAW PRIMARY MESH ON TOP
    for eh in mesh.edges():
        is_redundant_edge       = False                                             # Set boolean value, so it is always False in every iteration, unless edge belongs to *middle part* of boundary film in foam
        is_double_boundary_film = False                                             # Set boolean value, so it is always False in every iteration, unless edge belongs to *end* of boundary film in foam

        if mesh.is_boundary(eh) == False:                                           # Draw all films inside foam and boundary films not belonging to boundary face
            f0_h = mesh.face_handle(mesh.halfedge_handle(eh, 0))                    # Face handle for the one neighboring face
            f1_h = mesh.face_handle(mesh.halfedge_handle(eh, 1))                    # Face handle for the other neighboring face
            v0_h = mesh.to_vertex_handle(mesh.halfedge_handle(eh, 0))               # Vertex handle for the vertex at the end
            v1_h = mesh.to_vertex_handle(mesh.halfedge_handle(eh, 1))               # Vertex handle for other vertex at the end

            if mesh.is_boundary(v0_h) != mesh.is_boundary(v1_h):
                # Check, if edge belongs to boundary film in foam
                is_double_boundary_film = mesh.is_boundary(f0_h) != mesh.is_boundary(f1_h)

                # Check, if edge is redundant edge in the middle of a boundary film spanning more than 2 edges
                is_redundant_edge       = mesh.is_boundary(f0_h) and mesh.is_boundary(f1_h)

            if is_redundant_edge == False or draw_boundary==False:


                if is_double_boundary_film == False or draw_boundary==False:        # If film belongs to edge pointing to boundary, but not belonging to boundary face
                    delta_p = mesh.property(pressure_handle, v1_h) - mesh.property(pressure_handle, v0_h)

                    coor_0 = (mesh.property(mid_point_handle, f0_h)[0]*scale+offset, mesh.property(mid_point_handle, f0_h)[1]*scale+offset)
                    coor_1 = (mesh.property(mid_point_handle, f1_h)[0]*scale+offset, mesh.property(mid_point_handle, f1_h)[1]*scale+offset)
                else:                                                               # Treats special case, if boundary film spans several edges in dual film
                    if mesh.is_boundary(v1_h) == False:
                        temp = v1_h
                        v1_h = v0_h
                        v0_h = temp

                    delta_p = mesh.property(pressure_handle, v1_h) - mesh.property(pressure_handle, v0_h)

                    # Fetch coordinate for the one end of the multi-edge-crossing bubble film
                    fh_next  = f0_h if mesh.is_boundary(f0_h) else f1_h
                    fh_start = f0_h if not mesh.is_boundary(f0_h) else f1_h
                    coor_0   = (mesh.property(mid_point_handle,fh_start)[0]*scale+offset, mesh.property(mid_point_handle,fh_start)[1]*scale+offset)

                    # Calculate vector in boundry facing edge, describing the other non boundary facing edge in that face
                    current_edge_idx = eh.idx()
                    eh_last     = [eh0 for eh0 in mesh.fe(fh_next) if mesh.is_boundary(eh0)==False and current_edge_idx!=eh0.idx()][0]

                    eh_last_v0  = mesh.to_vertex_handle(mesh.halfedge_handle(eh_last, 0))
                    eh_last_v1  = mesh.to_vertex_handle(mesh.halfedge_handle(eh_last, 1))
                    (b1, a1)    = (mesh.point(eh_last_v1), mesh.point(eh_last_v0)) if mesh.is_boundary(eh_last_v1) else (mesh.point(eh_last_v0), mesh.point(eh_last_v1))
                    vec_1       = np.array([b1[0] - a1[0], -(b1[1] - a1[1])])

                    # Visit neighbouring faces to find the other coordinate
                    while mesh.is_boundary(fh_next):
                        for fh in mesh.ff(fh_next):
                            if fh != fh_start:
                                fh_start = fh_next                                  # Move one step forward
                                fh_next  = fh                                       # Move one step forward
                                break
                    coor_1 = (mesh.property(mid_point_handle,fh_next)[0]*scale+offset, mesh.property(mid_point_handle,fh_next)[1]*scale+offset)

                    ## Calculate edge vector in order to calculate cross product of the two edges. This way, the order of the edges can be calculated
                    # Calculate vector to current edge
                    (b0, a0)    = (mesh.point(v1_h),mesh.point(v0_h)) if mesh.is_boundary(v1_h) else (mesh.point(v0_h),mesh.point(v1_h))
                    vec_0       = np.array([b0[0] - a0[0], -(b0[1] - a0[1])])

                    order       = np.cross(vec_0,vec_1)

                    # Swap coor_0 and coor_1, if delta_p and the order of the coordinates don't match
                    if order>0:
                        temp   = coor_0
                        coor_0 = coor_1
                        coor_1 = temp

                # Draw between coordinates
                arc_path = dwg.path(d='M %f , %f' % coor_0, fill='none', stroke='black')

                if abs(delta_p) > 0:
                    rx = 2.0 * gamma / abs(delta_p) / scale

                    if delta_p > 0:
                        arc_path.push_arc(coor_1, 0, rx, False, '-', True)
                    else:
                        arc_path.push_arc(coor_1, 0, rx, False, '+', True)

                    dwg.add(arc_path)
                else:
                    dwg.add(dwg.line(coor_0, coor_1, stroke=svgwrite.rgb(0, 0, 0, '%')))
        else:
            pass


    for fh in mesh.faces():
        if fh.idx()==1019:
            coor = (mesh.property(mid_point_handle,fh)[0]*scale+offset,mesh.property(mid_point_handle,fh)[1]*scale+offset)            # Fetch coordinate for point
            dwg.add(dwg.circle(center=coor, r=2))

    dwg.save()


def draw_primary_mesh1(mesh,title='primary_mesh',style='edges',draw_boundary=True,scale=50, gamma=73.0):
    """ Draw dual mesh and write svg
        TriMesh mesh: an existing dual OpenMesh mesh
        string title: title of svg-file to be saved to
        string style: 'edges' or 'edges+vertices' for mesh drawing style
        bool draw_boundary: specifies, if boundary has to be drawn or not
        int scale: scale defining scale factor of visualization
        double gamma: surface tension
    """

    ## Draw dual mesh
    dwg  = svgwrite.Drawing('images/'+title+'.svg')

    # DRAW DUAL MESH
    for eh in mesh.edges():

        h0_h = mesh.halfedge_handle(eh,0)                                       # Fetch handle for halfedge in edge
        h1_h = mesh.halfedge_handle(eh,1)
        v0_h = mesh.to_vertex_handle(h0_h)                                      # Converts end point of halfedge handle to vertex handle
        v1_h = mesh.to_vertex_handle(h1_h)

        v0   = (mesh.point(v0_h)[0]*scale,mesh.point(v0_h)[1]*scale)            # Convert to tuple for coordinate
        v1   = (mesh.point(v1_h)[0]*scale,mesh.point(v1_h)[1]*scale)            # Convert to tuple for other coordinate

        dwg.add(dwg.line(v0, v1, stroke=svgwrite.rgb(80, 80, 80, '%')))         # Add to svg-file


    # DRAW PRIMARY MESH ON TOP
    for eh in mesh.edges():
        is_double_boundary_film = False                                         # Set boolean value, so it is always False in every iteration, unless edge belongs to *part* of boundary film in foam

        if mesh.is_boundary(eh) == False:                                       # Draw all films inside foam and boundary films not belonging to boundary face
            f0_h = mesh.face_handle(mesh.halfedge_handle(eh, 0))                # Face handle for the one neighboring face
            f1_h = mesh.face_handle(mesh.halfedge_handle(eh, 1))                # Face handle for the other neighboring face
            v0_h = mesh.to_vertex_handle(mesh.halfedge_handle(eh, 0))           # Vertex handle for the vertex at the end
            v1_h = mesh.to_vertex_handle(mesh.halfedge_handle(eh, 1))           # Vertex handle for other vertex at the end

            # Check, if edge belongs to boundary film in foam
            if mesh.is_boundary(v0_h) or mesh.is_boundary(v1_h):
                is_double_boundary_film = mesh.is_boundary(f0_h) or mesh.is_boundary(f1_h)

            delta_p = mesh.property(pressure_handle, v1_h) - mesh.property(pressure_handle, v0_h)

            if is_double_boundary_film == False:                                # If film belongs to edge pointing to boundary, but not belonging to boundary face
                coor_0 = (mesh.property(mid_point_handle, f0_h)[0] * scale, mesh.property(mid_point_handle, f0_h)[1] * scale)
                coor_1 = (mesh.property(mid_point_handle, f1_h)[0] * scale, mesh.property(mid_point_handle, f1_h)[1] * scale)

        else:

            # Continue, until correct face is found. OpenMesh is not consistent in this case, therefore this step is necessary
            correct_face = False
            i            = 0
            while correct_face == False:
                if mesh.is_boundary(mesh.face_handle(mesh.halfedge_handle(eh, i))):
                    fb_h                   = mesh.face_handle(mesh.halfedge_handle(eh, i))
                    neighbour_face_counter = 0
                    for _ in mesh.ff(fb_h):                                     # Count edge neighbouring faces
                        neighbour_face_counter += 1
                        if neighbour_face_counter > 2:
                            break                                               # No need to count more than 2 (plus, sometimes it counts to inf)
                    if neighbour_face_counter == 2:
                        correct_face = True
                i += 1

            fn_h    = [fh for fh in mesh.ff(fb_h)]                                  # Edge neighboring faces to boundary face

            coor_0  = (mesh.property(mid_point_handle, fn_h[0])[0] * scale, mesh.property(mid_point_handle, fn_h[0])[1] * scale)
            coor_1  = (mesh.property(mid_point_handle, fn_h[1])[0] * scale, mesh.property(mid_point_handle, fn_h[1])[1] * scale)

            for vh in mesh.fv(fb_h):
                if mesh.is_boundary(vh) == False:
                    vin_h = vh                                                      # Handle for vertex inside foam
                else:
                    vbn_h = vh                                                      # Handle for boundary vertex

            delta_p = mesh.property(pressure_handle, vin_h) - mesh.property(pressure_handle, vbn_h)

        # Draw between coordinates
        if is_double_boundary_film == False:
            if abs(delta_p) > 0:
                rx = 2.0 * gamma / abs(delta_p) / scale
                arc_path = dwg.path(d='M %f , %f' % coor_0, fill='none', stroke='black')

                if delta_p > 0:
                    arc_path.push_arc(coor_1, 0, rx, False, '-', True)
                else:
                    arc_path.push_arc(coor_1, 0, rx, False, '+', True)

                dwg.add(arc_path)
            else:
                dwg.add(dwg.line(coor_0, coor_1, stroke=svgwrite.rgb(0, 0, 0, '%')))


    # Add vertices, if specified
    if style!='edges':
        for vh in mesh.vertices():
            coor = (mesh.point(vh)[0]*scale,mesh.point(vh)[1]*scale)            # Fetch coordinate for point
            dwg.add(dwg.circle(center=coor, r=2))

    dwg.save()

    # # Produce PDF
    # drawing = svg2rlg('images/'+title+'.svg')
    # renderPDF.drawToFile(drawing, 'images/'+title+'.pdf')


def draw_primary_mesh2(mesh,title='primary_mesh',style='edges',draw_boundary=True,scale=50, gamma=73.0):
    """ Draw dual mesh and write svg
        TriMesh mesh: an existing dual OpenMesh mesh
        string title: title of svg-file to be saved to
        string style: 'edges' or 'edges+vertices' for mesh drawing style
        bool draw_boundary: specifies, if boundary has to be drawn or not
        int scale: scale defining scale factor of visualization
        double gamma: surface tension
    """

    ## Draw dual mesh
    dwg  = svgwrite.Drawing('images/'+title+'.svg')

    # Draw edges
    for eh in mesh.edges():
        h0_h = mesh.halfedge_handle(eh,0)                                       # Fetch handle for halfedge in edge
        h1_h = mesh.halfedge_handle(eh,1)
        v0_h = mesh.to_vertex_handle(h0_h)                                      # Converts end point of halfedge handle to vertex handle
        v1_h = mesh.to_vertex_handle(h1_h)

        v0   = (mesh.point(v0_h)[0]*scale,mesh.point(v0_h)[1]*scale)            # Convert to tuple for coordinate
        v1   = (mesh.point(v1_h)[0]*scale,mesh.point(v1_h)[1]*scale)            # Convert to tuple for other coordinate

        dwg.add(dwg.line(v0, v1, stroke=svgwrite.rgb(80, 80, 80, '%')))         # Add to svg-file

    # Add vertices, if specified
    if style!='edges':
        for vh in mesh.vertices():
            coor = (mesh.point(vh)[0]*scale,mesh.point(vh)[1]*scale)            # Fetch coordinate for point
            dwg.add(dwg.circle(center=coor, r=2))


    ## Add primary mesh on top (each bubble film crosses one edge in dual mesh)
    for eh in mesh.edges():

        if mesh.is_boundary(eh)==False:                                                 # Restriction before world bubble is defined

            f0_h    = mesh.face_handle(mesh.halfedge_handle(eh, 0))                     # Face handle for the one neighboring face
            f1_h    = mesh.face_handle(mesh.halfedge_handle(eh, 1))                     # Face handle for the other neighboring face
            v0_h    = mesh.to_vertex_handle(mesh.halfedge_handle(eh, 0))                # Vertex handle for the vertex at the end
            v1_h    = mesh.to_vertex_handle(mesh.halfedge_handle(eh, 1))                # Vertex handle for other vertex at the end

            delta_p = mesh.property(pressure_handle,v1_h)-mesh.property(pressure_handle,v0_h)

            bf0_h   = mesh.property(boundary_neighbour_faces_handle,f0_h)               # Neighbour points to current face (if neighbour is boundary)
            bf1_h   = mesh.property(boundary_neighbour_faces_handle,f1_h)               # Neighbour points to current face (if neighbour is boundary)face

            # Identify coordinates, between which to draw a line
            if len(bf0_h) == 0 and len(bf1_h) == 0 or draw_boundary==False:             # If film is inside foam
                coor_0 = (mesh.property(mid_point_handle, f0_h)[0] * scale, mesh.property(mid_point_handle, f0_h)[1] * scale)
                coor_1 = (mesh.property(mid_point_handle, f1_h)[0] * scale, mesh.property(mid_point_handle, f1_h)[1] * scale)
            else:                                                                       # If film is on boundary of foam
                nf_h   = bf0_h if len(bf0_h) > 0 else bf1_h

                if delta_p >0:
                    coor_0 = (nf_h[0][0] * scale, nf_h[0][1] * scale)
                    coor_1 = (nf_h[1][0] * scale, nf_h[1][1] * scale)
                else:
                    coor_1 = (nf_h[0][0] * scale, nf_h[0][1] * scale)
                    coor_0 = (nf_h[1][0] * scale, nf_h[1][1] * scale)


            # Draw between coordinates
            if abs(delta_p) > 0:
                rx       = 2.0*gamma/abs(delta_p)/scale
                arc_path = dwg.path(d='M %f , %f' % coor_0, fill='none', stroke='black')

                if delta_p>0:
                    arc_path.push_arc(coor_1, 0, rx, False, '-', True)
                else:
                    arc_path.push_arc(coor_1, 0, rx, False, '+', True)

                dwg.add(arc_path)
            else:
                dwg.add(dwg.line(coor_0, coor_1,stroke=svgwrite.rgb(0, 0, 0, '%')))

    dwg.save()

    # # Produce PDF
    # drawing = svg2rlg('images/'+title+'.svg')
    # renderPDF.drawToFile(drawing, 'images/'+title+'.pdf')


# ----------------------------------------------------------------------------------------------------------------------
# Main function
# ----------------------------------------------------------------------------------------------------------------------

if __name__ == "__main__":
    # Read variables from config-file (Initially for random grid)
    min = config.length['min']
    max = config.length['max']

    #---------------------------------------------------------------------------
    # Initialise mesh
    #---------------------------------------------------------------------------
    foam_mesh    = create_random_mesh(min,max)
    add_pressure(foam_mesh, equal=False)

    #---------------------------------------------------------------------------
    # Do parameter test (to be moved to separate file)
    #---------------------------------------------------------------------------
    # comprehensive_parameter_test(foam_mesh)


    #---------------------------------------------------------------------------
    # Equilibriate
    #---------------------------------------------------------------------------

    combine_equilibriate_vertices(foam_mesh)
    add_pressure(foam_mesh,equal=False)
    # calc_bubble_areas(foam_mesh)
    # calc_film_lengths(foam_mesh)
    add_face_connectivity(foam_mesh)
    calc_plateau_angles_new(foam_mesh)

    #draw_primary_mesh(foam_mesh, title='test_2')








# ex = []
# for eh in foam_mesh.edges():
#     if foam_mesh.is_boundary(eh) == False:
#         f0_h = foam_mesh.face_handle(foam_mesh.halfedge_handle(eh, 0))                        # Face handle for the one neighboring face
#         f1_h = foam_mesh.face_handle(foam_mesh.halfedge_handle(eh, 1))                        # Face handle for the other neighboring face
#         v0_h = foam_mesh.to_vertex_handle(foam_mesh.halfedge_handle(eh, 0))                   # Vertex handle for the vertex at the end
#         v1_h = foam_mesh.to_vertex_handle(foam_mesh.halfedge_handle(eh, 1))                   # Vertex handle for other vertex at the end
#
#         is_double_boundary_film = foam_mesh.is_boundary(f0_h) != foam_mesh.is_boundary(f1_h)
#
#         if is_double_boundary_film:
#             ex.append(foam_mesh.property(film_length_handle,eh))
#
# ex.sort()
#
#
# ex = []
# for vh in foam_mesh.vertices():
#     if foam_mesh.is_boundary(vh) == False:
#         ex.append(foam_mesh.property(bubble_area_handle,vh))
#
# ex.sort()

    for fh in foam_mesh.faces():
        if not foam_mesh.is_boundary(fh):
            ex = np.array(foam_mesh.property(plateau_angle_handle,fh))
            print sum(ex[:,0])



