######################################################################
# Automatically generated by qmake (3.0) mer. nov. 6 14:03:42 2019
######################################################################

TEMPLATE = app
TARGET = project
INCLUDEPATH += . \
               src/external/perlin \
               src/lib/3d \
               src/lib/common \
               src/lib/interface \
               src/lib/intersection \
               src/lib/mesh \
               src/lib/mesh/format \
               src/lib/opengl \
               src/lib/perlin \
               src/local/interface \
               src/local/scene

# Input
HEADERS += src/external/perlin/simplexnoise1234.hpp \
           src/lib/3d/mat1x4.hpp \
           src/lib/3d/mat2.hpp \
           src/lib/3d/mat3.hpp \
           src/lib/3d/mat4.hpp \
           src/lib/3d/mat4x1.hpp \
           src/lib/3d/quaternion.hpp \
           src/lib/3d/vec2.hpp \
           src/lib/3d/vec3.hpp \
           src/lib/3d/vec4.hpp \
           src/lib/common/backtrace.hpp \
           src/lib/common/error_handling.hpp \
           src/lib/common/exception_cpe.hpp \
           src/lib/interface/application_qt.hpp \
           src/lib/interface/camera_matrices.hpp \
           src/lib/interface/navigator_tool.hpp \
           src/lib/interface/picking_data.hpp \
           src/lib/interface/selected_index.hpp \
           src/lib/interface/trackball.hpp \
           src/lib/intersection/intersection.hpp \
           src/lib/mesh/mesh.hpp \
           src/lib/mesh/mesh_basic.hpp \
           src/lib/mesh/mesh_io.hpp \
           src/lib/mesh/mesh_parametric.hpp \
           src/lib/mesh/triangle_index.hpp \
           src/lib/opengl/axes_helper.hpp \
           src/lib/opengl/glutils.hpp \
           src/lib/opengl/line_opengl.hpp \
           src/lib/opengl/mesh_opengl.hpp \
           src/lib/perlin/perlin.hpp \
           src/local/interface/myWidgetGL.hpp \
           src/local/interface/myWindow.hpp \
           src/local/scene/scene.hpp \
           src/lib/mesh/format/mesh_io_obj.hpp \
           src/lib/mesh/format/mesh_io_off.hpp
FORMS += src/local/interface/mainwindow.ui
SOURCES += src/external/perlin/simplexnoise1234.cpp \
           src/lib/3d/mat1x4.cpp \
           src/lib/3d/mat2.cpp \
           src/lib/3d/mat3.cpp \
           src/lib/3d/mat4.cpp \
           src/lib/3d/mat4x1.cpp \
           src/lib/3d/quaternion.cpp \
           src/lib/3d/vec2.cpp \
           src/lib/3d/vec3.cpp \
           src/lib/3d/vec4.cpp \
           src/lib/common/backtrace.cpp \
           src/lib/common/exception_cpe.cpp \
           src/lib/interface/application_qt.cpp \
           src/lib/interface/navigator_tool.cpp \
           src/lib/interface/picking_data.cpp \
           src/lib/interface/selected_index.cpp \
           src/lib/interface/trackball.cpp \
           src/lib/intersection/intersection.cpp \
           src/lib/mesh/mesh.cpp \
           src/lib/mesh/mesh_basic.cpp \
           src/lib/mesh/mesh_io.cpp \
           src/lib/mesh/mesh_parametric.cpp \
           src/lib/mesh/triangle_index.cpp \
           src/lib/opengl/axes_helper.cpp \
           src/lib/opengl/glutils.cpp \
           src/lib/opengl/line_opengl.cpp \
           src/lib/opengl/mesh_opengl.cpp \
           src/lib/perlin/perlin.cpp \
           src/local/interface/main.cpp \
           src/local/interface/myWidgetGL.cpp \
           src/local/interface/myWindow.cpp \
           src/local/scene/scene.cpp \
           src/lib/mesh/format/mesh_io_obj.cpp \
           src/lib/mesh/format/mesh_io_off.cpp
