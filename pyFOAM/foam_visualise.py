from foam_import_openmesh import *
import foam_properties as fp
import config

# ----------------------------------------------------------------------------------------------------------------------
# Import constants
# ----------------------------------------------------------------------------------------------------------------------
scale_imp  = config.draw['scale']
offset_imp = config.draw['offset']
gamma_imp  = config.physical['gamma']
min_imp    = config.length['min']
max_imp    = config.length['max']


def draw_dual_mesh(mesh,title='dual_mesh',style='edges',scale=scale_imp, offset=offset_imp):
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
    if style == 'edges' or style != 'edges':
        for vh in mesh.vertices():
            coor = (mesh.point(vh)[0]*scale+offset,mesh.point(vh)[1]*scale+offset)      # Fetch coordinate for point
            dwg.add(dwg.circle(center=coor, r=2))


    dwg.save()


def draw_primary_mesh(mesh,title='primary_mesh',draw_boundary=True,scale=scale_imp, offset=offset_imp, gamma=gamma_imp, draw_dual = False, pressure_colour=False, min=min_imp, max=max_imp):
    """ Draw dual mesh and write svg
        TriMesh mesh: an existing dual OpenMesh mesh
        string title: title of svg-file to be saved to
        bool draw_boundary: specifies, if boundary has to be drawn or not
        int scale: scale defining scale factor of visualization
        int offset: added to move drawing away from screen edges
        double gamma: surface tension
        bool draw_dual: draw dual mesh or not
        bool pressure_colour: draw colour map indicating pressures
    """

    # Initialise drawing
    dwg  = svgwrite.Drawing('images/'+title+'.svg')

    # Draw white corners of image in white to keep scaling of image
    coor_0 = (min * scale + offset-10, min * scale + offset-10)
    coor_1 = (min * scale + offset-10, max * scale + offset+10)
    coor_2 = (max * scale + offset+10, min * scale + offset-10)
    coor_3 = (max * scale + offset+10, max * scale + offset+10)
    dwg.add(dwg.circle(center=coor_0, r=1, fill='white'))
    dwg.add(dwg.circle(center=coor_1, r=1, fill='white'))
    dwg.add(dwg.circle(center=coor_2, r=1, fill='white'))
    dwg.add(dwg.circle(center=coor_3, r=1, fill='white'))

    # Draw dual mesh
    if draw_dual:
        for eh in mesh.edges():

            h0_h = mesh.halfedge_handle(eh,0)                                           # Fetch handle for halfedge in edge
            h1_h = mesh.halfedge_handle(eh,1)
            v0_h = mesh.to_vertex_handle(h0_h)                                          # Converts end point of halfedge handle to vertex handle
            v1_h = mesh.to_vertex_handle(h1_h)

            v0   = (mesh.point(v0_h)[0]*scale+offset,mesh.point(v0_h)[1]*scale+offset)  # Convert to tuple for coordinate
            v1   = (mesh.point(v1_h)[0]*scale+offset,mesh.point(v1_h)[1]*scale+offset)  # Convert to tuple for other coordinate

            dwg.add(dwg.line(v0, v1, stroke=svgwrite.rgb(80, 80, 80, '%')))             # Add to svg-file


    # Colour according to pressure, and draw mesh lines on top afterwards
    if pressure_colour:
        for vh in mesh.vertices():
            if not mesh.is_boundary(vh):
                i = 0
                pressure = mesh.property(fp.pressure_handle, vh)
                if pressure > 1.8: pressure = 1.8
                red      = (pressure - 0.90) / 0.9 * 255
                blue     = 255 - red
                for fh in mesh.vf(vh):
                    if not mesh.is_boundary(fh):
                        coor = (mesh.property(fp.mid_point_handle, fh)[0] * scale + offset, mesh.property(fp.mid_point_handle, fh)[1] * scale + offset)
                        if i == 0:
                            p = dwg.path(d='M %f , %f' % coor, fill=svgwrite.rgb(red, 0, blue), stroke='none')
                            i += 1
                        else:
                            p.push(coor)

                        dwg.add(p)


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

            # Initialise check, for when dual vertices inside foam have moved beyond boundary vertices
            order2 = -1

            if is_redundant_edge == False or draw_boundary==False:

                if is_double_boundary_film == False or draw_boundary==False:        # If film belongs to edge pointing to boundary, but not belonging to boundary face
                    delta_p = mesh.property(fp.pressure_handle, v1_h) - mesh.property(fp.pressure_handle, v0_h)
                    if abs(delta_p) < 1e-8: delta_p = 0                             # Prevent overflow in division, which results in line not being drawn

                    coor_0 = (mesh.property(fp.mid_point_handle, f0_h)[0]*scale+offset, mesh.property(fp.mid_point_handle, f0_h)[1]*scale+offset)
                    coor_1 = (mesh.property(fp.mid_point_handle, f1_h)[0]*scale+offset, mesh.property(fp.mid_point_handle, f1_h)[1]*scale+offset)
                else:                                                               # Treats special case, if boundary film spans several edges in dual film
                    if mesh.is_boundary(v1_h) == False:
                        temp = v1_h
                        v1_h = v0_h
                        v0_h = temp

                    delta_p = mesh.property(fp.pressure_handle, v1_h) - mesh.property(fp.pressure_handle, v0_h)
                    if abs(delta_p) < 1e-8: delta_p = 0                             # Prevent overflow in division, which results in line not being drawn


                    # Fetch coordinate for the one end of the multi-edge-crossing bubble film
                    fh_next  = f0_h if mesh.is_boundary(f0_h) else f1_h
                    fh_boundary = fh_next
                    fh_start = f0_h if not mesh.is_boundary(f0_h) else f1_h
                    coor_0   = (mesh.property(fp.mid_point_handle,fh_start)[0]*scale+offset, mesh.property(fp.mid_point_handle,fh_start)[1]*scale+offset)

                    # Calculate vector in boundry facing edge, describing the other non boundary facing edge in that face
                    current_edge_idx = eh.idx()
                    eh_last     = [eh0 for eh0 in mesh.fe(fh_next) if mesh.is_boundary(eh0)==False and current_edge_idx!=eh0.idx()][0]

                    eh_last_v0  = mesh.to_vertex_handle(mesh.halfedge_handle(eh_last, 0))
                    eh_last_v1  = mesh.to_vertex_handle(mesh.halfedge_handle(eh_last, 1))
                    (b1, a1)    = (mesh.point(eh_last_v1), mesh.point(eh_last_v0)) if mesh.is_boundary(eh_last_v1) else (mesh.point(eh_last_v0), mesh.point(eh_last_v1))
                    vec_1       = np.array([b1[0] - a1[0], b1[1] - a1[1]])

                    # Visit neighbouring faces to find the other coordinate
                    while mesh.is_boundary(fh_next):
                        for fh in mesh.ff(fh_next):
                            if fh != fh_start:
                                fh_start = fh_next                                  # Move one step forward
                                fh_next  = fh                                       # Move one step forward
                                break
                    coor_1 = (mesh.property(fp.mid_point_handle,fh_next)[0]*scale+offset, mesh.property(fp.mid_point_handle,fh_next)[1]*scale+offset)

                    ## Calculate edge vector in order to calculate cross product of the two edges. This way, the order of the edges can be calculated
                    # Calculate vector to current edge
                    (b0, a0)    = (mesh.point(v1_h),mesh.point(v0_h)) if mesh.is_boundary(v1_h) else (mesh.point(v0_h),mesh.point(v1_h))
                    vec_0       = np.array([b0[0] - a0[0], b0[1] - a0[1]])

                    order       = np.cross(vec_0,vec_1)

                    # Create second check, if dual vertices inside foam have moved beyond boundary vertices in update_vertex_positions
                    for eh_boundary in mesh.fe(fh_boundary):
                        eh_boundary_face = eh_boundary
                        break
                    if not mesh.is_boundary(mesh.halfedge_handle(eh_boundary_face, 0)):
                        heh0        = mesh.halfedge_handle(eh_boundary_face, 0)
                        heh1        = mesh.halfedge_handle(eh_boundary_face, 1)
                        heh2        = mesh.next_halfedge_handle(heh0)
                    else:
                        heh0        = mesh.halfedge_handle(eh_boundary_face, 1)
                        heh1        = mesh.halfedge_handle(eh_boundary_face, 0)
                        heh2        = mesh.next_halfedge_handle(heh0)
                    a           = mesh.point(mesh.to_vertex_handle(heh1))
                    b           = mesh.point(mesh.to_vertex_handle(heh0))
                    #b           = np.array([3.75, 3.75])
                    c           = mesh.point(mesh.to_vertex_handle(heh2))
                    order2      = np.cross(np.array([a[0] - b[0], a[1] - b[1]]), np.array([c[0] - b[0], c[1] - b[1]]))

                    # Swap coor_0 and coor_1, if delta_p and the order of the coordinates don't match
                    if order<0:
                        temp   = coor_0
                        coor_0 = coor_1
                        coor_1 = temp

                # Draw between coordinates
                if abs(delta_p) > 0:
                    if pressure_colour:
                        if delta_p > 0:
                            pressure = mesh.property(fp.pressure_handle, v1_h)
                        else:
                            pressure = mesh.property(fp.pressure_handle, v0_h)
                        if pressure > 1.8: pressure = 1.8
                        red = (pressure - 0.9) / 0.9 * 255
                        blue = 255 - red
                        arc_path = dwg.path(d='M %f , %f' % coor_0, fill=svgwrite.rgb(red, 0, blue), stroke='black')
                    else:
                        arc_path = dwg.path(d='M %f , %f' % coor_0, fill='none', stroke='black')

                    rx = 2.0 * gamma / abs(delta_p) * scale

                    # if v0_h.idx()==87 or v1_h.idx()==87:
                    #     print rx, 0.5*np.sqrt((coor_0[0]-coor_1[0])**2+(coor_0[1]-coor_1[1])**2), delta_p, mesh.property(fp.pressure_handle, v1_h), mesh.property(fp.pressure_handle, v0_h)

                    #rx = 2.0 * gamma / delta_p * scale

                    if delta_p > 0:
                        if order2 < 0:
                            arc_path.push_arc(coor_1, 0, rx, False, '-', True)
                        else:                                                                                           # If dual vertices inside foam have moved beyond boundary vertices
                            arc_path.push_arc(coor_1, 0, rx, False, '+', True)
                    else:
                        if order2 < 0:
                            arc_path.push_arc(coor_1, 0, rx, False, '+', True)
                        else:                                                                                           # If dual vertices inside foam have moved beyond boundary vertices
                            arc_path.push_arc(coor_1, 0, rx, False, '-', True)

                    dwg.add(arc_path)
                else:
                    dwg.add(dwg.line(coor_0, coor_1, stroke=svgwrite.rgb(0, 0, 0, '%')))
        else:
            pass



    # for vh in mesh.vertices():
    #     if vh.idx() == 105:
    #         coor = (mesh.point(vh)[0]*scale+offset,mesh.point(vh)[1]*scale+offset)  # Fetch coordinate for point
    #         dwg.add(dwg.circle(center=coor, r=3))
    #     elif vh.idx() == 63:
    #         coor = (mesh.point(vh)[0]*scale+offset,mesh.point(vh)[1]*scale+offset)  # Fetch coordinate for point
    #         dwg.add(dwg.circle(center=coor, r=2))
    #     elif vh.idx() == 37:
    #         coor = (mesh.point(vh)[0] * scale + offset, mesh.point(vh)[1] * scale + offset)  # Fetch coordinate for point
    #         dwg.add(dwg.circle(center=coor, r=1 ))
    #
    # coor = (6.15063445 * scale + offset, 6.46762186 * scale + offset)  # Fetch coordinate for point
    # dwg.add(dwg.circle(center=coor, r=3))
    # coor = (6.63375577 * scale + offset, 6.02992085 * scale + offset)  # Fetch coordinate for point
    # dwg.add(dwg.circle(center=coor, r=2))

    dwg.save()

