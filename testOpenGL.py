import numpy as np
import pyqtgraph.opengl as gl
from PyQt5.QtWidgets import QApplication
from PyQt5.QtGui import QSurfaceFormat
import sys

# Force OpenGL Compatibility Profile
fmt = QSurfaceFormat()
fmt.setRenderableType(QSurfaceFormat.OpenGL)
fmt.setProfile(QSurfaceFormat.CompatibilityProfile)  # Allows legacy OpenGL
fmt.setVersion(2, 1)  # Use OpenGL 2.1 (widely compatible)
QSurfaceFormat.setDefaultFormat(fmt)

app = QApplication(sys.argv)

# Create GLViewWidget window
w = gl.GLViewWidget()
w.show()

# Add a 3D grid for reference
grid = gl.GLGridItem()
w.addItem(grid)

# Define simple tetrahedron mesh
verts = np.array([
    [0, 0, 0],
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 1]
], dtype=np.float32)

faces = np.array([
    [0, 1, 2],
    [0, 1, 3],
    [0, 2, 3],
    [1, 2, 3]
], dtype=np.int32)

colors = np.ones((faces.shape[0], 4), dtype=np.float32)  # RGBA colors

print("Vertices:\n", verts)
print("Faces:\n", faces)

# Create mesh item
mesh = gl.GLMeshItem(vertexes=verts, faces=faces, faceColors=colors, drawEdges=True)
print("Adding mesh...")
w.addItem(mesh)

# Start Qt application event loop
sys.exit(app.exec_())


# import pyqtgraph.opengl as gl
# from PyQt5.QtWidgets import QApplication
# import sys

# from PyQt5.QtGui import QSurfaceFormat

# fmt = QSurfaceFormat()
# fmt.setRenderableType(QSurfaceFormat.OpenGL)
# fmt.setProfile(QSurfaceFormat.CompatibilityProfile)  # Allows legacy OpenGL
# fmt.setVersion(2, 1)  # Use OpenGL 2.1 (widely compatible)
# QSurfaceFormat.setDefaultFormat(fmt)

# app = QApplication(sys.argv)
# w = gl.GLViewWidget()
# w.show()

# # Add a 3D grid
# grid = gl.GLGridItem()
# w.addItem(grid)
# sys.exit(app.exec_())

# import numpy as np
# import pyqtgraph.opengl as gl
# from PyQt5.QtWidgets import QApplication
# from PyQt5.QtGui import QSurfaceFormat

# # Force OpenGL Compatibility Profile
# fmt = QSurfaceFormat()
# fmt.setRenderableType(QSurfaceFormat.OpenGL)
# fmt.setProfile(QSurfaceFormat.CompatibilityProfile)
# fmt.setVersion(2, 1)  
# QSurfaceFormat.setDefaultFormat(fmt)

# app = QApplication([])

# w = gl.GLViewWidget()
# w.show()

# # Add grid for reference
# grid = gl.GLGridItem()
# w.addItem(grid)

# # Define a simple tetrahedron
# verts = np.array([
#     [0, 0, 0],
#     [1, 0, 0],
#     [0, 1, 0],
#     [0, 0, 1]
# ], dtype=np.float32)

# faces = np.array([
#     [0, 1, 2],
#     [0, 1, 3],
#     [0, 2, 3],
#     [1, 2, 3]
# ], dtype=np.int32)

# colors = np.ones((faces.shape[0], 4), dtype=np.float32)  # RGBA colors

# print("Vertices:\n", verts)
# print("Faces:\n", faces)

# # Create mesh item
# mesh = gl.GLMeshItem(vertexes=verts, faces=faces, faceColors=colors, drawEdges=True)
# print("Adding mesh...")
# w.addItem(mesh)

# app.exec_()

# # import numpy as np
# # import pyqtgraph.opengl as gl

# # verts = np.array([
# #     [0, 0, 0],
# #     [1, 0, 0],
# #     [0, 1, 0],
# #     [0, 0, 1]
# # ])

# # faces = np.array([
# #     [0, 1, 2],
# #     [0, 1, 3],
# #     [0, 2, 3],
# #     [1, 2, 3]
# # ])

# # colors = np.ones((faces.shape[0], 4))  # RGBA colors

# # print("Vertices:\n", verts)
# # print("Faces:\n", faces)


# # mesh = gl.GLMeshItem(vertexes=verts, faces=faces, faceColors=colors, drawEdges=True)

# # app = QApplication([])
# # w = gl.GLViewWidget()
# # w.show()
# # grid = gl.GLGridItem()
# # w.addItem(grid)
# # w.addItem(mesh)
# # app.exec_()



