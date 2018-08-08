include("${CMAKE_CURRENT_LIST_DIR}/blender_lite.cmake")

set(WITH_PYTHON_INSTALL_NUMPY           OFF CACHE BOOL "" FORCE)
set(WITH_SYSTEM_GLEW                    OFF CACHE BOOL "" FORCE)
# Brecht says: KEEP THIS ENABLED. *SOMETHING* breaks otherwise.
# Something == DPI. Go figure.
set(WITH_GHOST_XDND                     ON  CACHE BOOL "" FORCE)
