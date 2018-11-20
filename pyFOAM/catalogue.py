from openmesh_example import *

#---------------------------------------------------------------------------
# Create catalogue (to be moved to separate file)
#---------------------------------------------------------------------------


def two_bubble():
    mesh = TriMesh()

    vh0 = mesh.add_vertex(TriMesh.Point(5.5, 5, 0))
    vh1 = mesh.add_vertex(TriMesh.Point(6.5, 5, 0))
    vh2 = mesh.add_vertex(TriMesh.Point(6, 4, 0))
    vh3 = mesh.add_vertex(TriMesh.Point(6, 6, 0))
    vh4 = mesh.add_vertex(TriMesh.Point(5, 4, 0))
    vh5 = mesh.add_vertex(TriMesh.Point(5, 6, 0))
    vh6 = mesh.add_vertex(TriMesh.Point(7, 4, 0))
    vh7 = mesh.add_vertex(TriMesh.Point(7, 6, 0))

    fh0 = mesh.add_face(vh0, vh2, vh1)
    fh1 = mesh.add_face(vh0, vh1, vh3)
    fh2 = mesh.add_face(vh4, vh2, vh0)
    fh3 = mesh.add_face(vh5, vh4, vh0)
    fh4 = mesh.add_face(vh3, vh5, vh0)
    fh5 = mesh.add_face(vh2, vh6, vh1)
    fh6 = mesh.add_face(vh1, vh6, vh7)
    fh7 = mesh.add_face(vh3, vh1, vh7)

    add_pressure(mesh, equal=False, world_pressure_fraction=0.1)
    add_mid_point(mesh)

    draw_primary_mesh(mesh, title='2-bubble_equal')


def three_bubble():
    mesh = TriMesh()

    vh0  = mesh.add_vertex(TriMesh.Point(5, 5, 0))
    vh1  = mesh.add_vertex(TriMesh.Point(6, 5, 0))
    vh2  = mesh.add_vertex(TriMesh.Point(7, 5, 0))
    vh3  = mesh.add_vertex(TriMesh.Point(5.5, 6, 0))
    vh4  = mesh.add_vertex(TriMesh.Point(6.5, 6, 0))
    vh5  = mesh.add_vertex(TriMesh.Point(6, 7, 0))
    vh6  = mesh.add_vertex(TriMesh.Point(5.5, 4, 0))
    vh7  = mesh.add_vertex(TriMesh.Point(6.5, 4, 0))
    vh8  = mesh.add_vertex(TriMesh.Point(4.5, 6, 0))
    vh9  = mesh.add_vertex(TriMesh.Point(7.5, 6, 0))
    vh10 = mesh.add_vertex(TriMesh.Point(5.5, 7, 0))
    vh11 = mesh.add_vertex(TriMesh.Point(6.5, 7, 0))

    fh0  = mesh.add_face(vh0, vh1, vh3)
    fh1  = mesh.add_face(vh1, vh4, vh3)
    fh2  = mesh.add_face(vh1, vh2, vh4)
    fh3  = mesh.add_face(vh3, vh4, vh5)
    fh4  = mesh.add_face(vh6, vh7, vh1)
    fh5  = mesh.add_face(vh1, vh7, vh2)
    fh6  = mesh.add_face(vh2, vh9, vh4)
    fh7  = mesh.add_face(vh4, vh9, vh11)
    fh8  = mesh.add_face(vh5, vh4, vh11)
    fh9  = mesh.add_face(vh3, vh5, vh10)
    fh10 = mesh.add_face(vh3, vh10, vh8)
    fh11 = mesh.add_face(vh0, vh3, vh8)
    fh11 = mesh.add_face(vh0, vh6, vh1)

    add_pressure(mesh, equal=False, world_pressure_fraction=0.9)
    add_mid_point(mesh)

    draw_primary_mesh(mesh, title='3-bubble')


def three_bubble_chain():
    mesh = TriMesh()

    vh0  = mesh.add_vertex(TriMesh.Point(5.5, 5, 0))
    vh1  = mesh.add_vertex(TriMesh.Point(6.5, 5, 0))
    vh2  = mesh.add_vertex(TriMesh.Point(6, 3.5, 0))
    vh3  = mesh.add_vertex(TriMesh.Point(6, 6.5, 0))
    vh4  = mesh.add_vertex(TriMesh.Point(5, 4, 0))
    vh5  = mesh.add_vertex(TriMesh.Point(5, 6, 0))
    vh6  = mesh.add_vertex(TriMesh.Point(7, 3.5, 0))
    vh7  = mesh.add_vertex(TriMesh.Point(7, 6.5, 0))
    vh8  = mesh.add_vertex(TriMesh.Point(7.5, 5, 0))
    vh9  = mesh.add_vertex(TriMesh.Point(8, 4, 0))
    vh10 = mesh.add_vertex(TriMesh.Point(8, 6, 0))


    fh0  = mesh.add_face(vh0, vh2, vh1)
    fh1  = mesh.add_face(vh0, vh1, vh3)
    fh2  = mesh.add_face(vh4, vh2, vh0)
    fh3  = mesh.add_face(vh5, vh4, vh0)
    fh4  = mesh.add_face(vh3, vh5, vh0)
    fh5  = mesh.add_face(vh1, vh2, vh6)
    fh6  = mesh.add_face(vh1, vh7, vh3)
    fh7  = mesh.add_face(vh1, vh6, vh8)
    fh8  = mesh.add_face(vh1, vh8, vh7)
    fh9  = mesh.add_face(vh8, vh9, vh10)
    fh10 = mesh.add_face(vh8, vh10, vh7)
    fh11 = mesh.add_face(vh8, vh6, vh9)


    add_pressure(mesh, equal=False, world_pressure_fraction=0.9)
    add_mid_point(mesh)

    draw_primary_mesh(mesh, title='3-bubble_chain')


def H_bubble():
    mesh = TriMesh()

    vh0  = mesh.add_vertex(TriMesh.Point(5, 5, 0))
    vh1  = mesh.add_vertex(TriMesh.Point(6, 4, 0))
    vh2  = mesh.add_vertex(TriMesh.Point(6, 6, 0))
    vh3  = mesh.add_vertex(TriMesh.Point(7, 5, 0))
    vh4  = mesh.add_vertex(TriMesh.Point(4, 3, 0))
    vh5  = mesh.add_vertex(TriMesh.Point(8, 3, 0))
    vh6  = mesh.add_vertex(TriMesh.Point(4, 7, 0))
    vh7  = mesh.add_vertex(TriMesh.Point(8, 7, 0))

    fh0  = mesh.add_face(vh0, vh1, vh2)
    fh1  = mesh.add_face(vh1, vh3, vh2)
    fh2  = mesh.add_face(vh0, vh4, vh1)
    fh3  = mesh.add_face(vh1, vh4, vh5)
    fh4  = mesh.add_face(vh1, vh5, vh3)
    fh5  = mesh.add_face(vh3, vh5, vh7)
    fh6  = mesh.add_face(vh3, vh7, vh2)
    fh7  = mesh.add_face(vh2, vh7, vh6)
    fh8  = mesh.add_face(vh0, vh2, vh6)
    fh9  = mesh.add_face(vh4, vh0, vh6)

    add_pressure(mesh, equal=False, world_pressure_fraction=0.9)
    add_mid_point(mesh)

    draw_primary_mesh(mesh, title='H-bubble')


def H_bubble_v2():
    mesh = TriMesh()

    vh0  = mesh.add_vertex(TriMesh.Point(5, 5, 0))
    vh1  = mesh.add_vertex(TriMesh.Point(6, 4, 0))
    vh2  = mesh.add_vertex(TriMesh.Point(6, 6, 0))
    vh3  = mesh.add_vertex(TriMesh.Point(7, 5, 0))
    vh4  = mesh.add_vertex(TriMesh.Point(4, 3, 0))
    vh5  = mesh.add_vertex(TriMesh.Point(8, 3, 0))
    vh6  = mesh.add_vertex(TriMesh.Point(4, 7, 0))
    vh7  = mesh.add_vertex(TriMesh.Point(8, 7, 0))
    vh8  = mesh.add_vertex(TriMesh.Point(4, 4.33, 0))
    vh9  = mesh.add_vertex(TriMesh.Point(4, 5.67, 0))
    vh10 = mesh.add_vertex(TriMesh.Point(8, 3.8, 0))
    vh11 = mesh.add_vertex(TriMesh.Point(8, 4.6, 0))
    vh12 = mesh.add_vertex(TriMesh.Point(8, 5.4, 0))
    vh13 = mesh.add_vertex(TriMesh.Point(8, 6.2, 0))

    fh0  = mesh.add_face(vh0, vh1, vh2)
    fh1  = mesh.add_face(vh1, vh3, vh2)
    fh2  = mesh.add_face(vh0, vh4, vh1)
    fh3  = mesh.add_face(vh1, vh4, vh5)
    fh4  = mesh.add_face(vh1, vh5, vh3)
    fh5  = mesh.add_face(vh3, vh7, vh2)
    fh6  = mesh.add_face(vh2, vh7, vh6)
    fh7  = mesh.add_face(vh0, vh2, vh6)
    fh8  = mesh.add_face(vh0, vh8, vh4)
    fh9  = mesh.add_face(vh0, vh9, vh8)
    fh10 = mesh.add_face(vh0, vh6, vh9)
    fh11 = mesh.add_face(vh3, vh5, vh10)
    fh12 = mesh.add_face(vh3, vh10, vh11)
    fh13 = mesh.add_face(vh3, vh11, vh12)
    fh14 = mesh.add_face(vh3, vh12, vh13)
    fh15 = mesh.add_face(vh3, vh13, vh7)

    add_pressure(mesh, equal=False, world_pressure_fraction=0.9)
    add_mid_point(mesh)

    draw_primary_mesh(mesh, title='H-bubble_v2')


if __name__ == "__main__":
    # two_bubble()
    # three_bubble()
    # three_bubble_chain()
    # H_bubble()
    # H_bubble_v2()

    # mesh = TriMesh()
    # vh_0 = mesh.add_vertex(TriMesh.Point(0, 0, 0))
    # vh_1 = mesh.add_vertex(TriMesh.Point(1, 0, 0))
    # vh_2 = mesh.add_vertex(TriMesh.Point(0, 1, 0))
    # vh_3 = mesh.add_vertex(TriMesh.Point(1, 1, 0))
    # fh_0 = mesh.add_face(vh_0, vh_1, vh_2)
    # fh_1 = mesh.add_face(vh_1, vh_3, vh_2)
    #
    # draw_dual_mesh(mesh, style='edges', title='dual_mesh_0')
    #
    # mesh.request_face_status()
    # mesh.request_vertex_status()
    # mesh.request_edge_status()
    # mesh.delete_face(fh_0)
    # mesh.garbage_collection()
    #
    # draw_dual_mesh(mesh, style='edges', title='dual_mesh_1')