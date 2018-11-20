import sys

sys.path.append("/Users/jtbondorf/bin/OpenMesh/build/Build/python")
sys.path.append("/Users/jtbondorf/Documents/Universitet/Speciale/2016_JOHAN_TEICHERT_BONDORF/src")

from openmesh import *
from openmesh_example import *

mesh = TriMesh()

vh0 = mesh.add_vertex(TriMesh.Point(1, 1, 0))
vh1 = mesh.add_vertex(TriMesh.Point(1, 2, 0))
vh2 = mesh.add_vertex(TriMesh.Point(2, 1.5, 0))
vh3 = mesh.add_vertex(TriMesh.Point(1.5,1.5, 0))
vh4 = mesh.add_vertex(TriMesh.Point(2,2, 0))

fh0 = mesh.add_face(vh0, vh1, vh3)
fh1 = mesh.add_face(vh1, vh2, vh3)
fh2 = mesh.add_face(vh2, vh0, vh3)
fh3 = mesh.add_face(vh4, vh2, vh1)


for eh in mesh.edges():
    if mesh.is_boundary(eh):
        print mesh.halfedge_handle(eh)
        print mesh.halfedge_handle(eh, 0), mesh.face_handle(mesh.halfedge_handle(eh, 0))
        print mesh.halfedge_handle(eh, 0), mesh.face_handle(mesh.halfedge_handle(eh, 0))
        print mesh.halfedge_handle(eh, 0), mesh.face_handle(mesh.halfedge_handle(eh, 0))
        print mesh.halfedge_handle(eh, 0), mesh.face_handle(mesh.halfedge_handle(eh, 0))
        break
        # print mesh.face_handle(mesh.next_halfedge_handle(mesh.halfedge_handle(eh, 0)))
        # print mesh.face_handle(mesh.next_halfedge_handle(mesh.next_halfedge_handle(mesh.halfedge_handle(eh, 0))))



# mesh.request_face_status()
# mesh.request_vertex_status()
# mesh.request_halfedge_status()
# mesh.request_edge_status()
#
draw_dual_mesh(mesh,title='test')
#
# for heh in mesh.voh(vh3):
#     if mesh.is_collapse_ok(heh):
#         mesh.collapse(heh)
#         break
#
# mesh.garbage_collection()
#
# draw_dual_mesh(mesh,title='test_after_T2_process')










