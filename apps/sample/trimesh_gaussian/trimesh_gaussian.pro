include(../common.pri)
TARGET = trimesh_gaussian
SOURCES += trimesh_gaussian.cpp ../../../wrap/ply/plylib.cpp

HEADERS += \
    export_ply_GS.h \
    gaussian_splat.h \
    import_ply_GS.h
