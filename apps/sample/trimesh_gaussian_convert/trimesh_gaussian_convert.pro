include(../common.pri)
TARGET = trimesh_gaussian_convert
SOURCES += trimesh_gaussian_convert.cpp ../../../wrap/ply/plylib.cpp

HEADERS += \
    ../trimesh_gaussian/export_ply_GS.h \
    ../trimesh_gaussian/export_splat.h \
    ../trimesh_gaussian/gaussian_splat.h \
    ../trimesh_gaussian/import_ply_GS.h \
    ../trimesh_gaussian/import_splat.h
