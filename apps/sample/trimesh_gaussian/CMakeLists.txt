cmake_minimum_required(VERSION 3.13)
project(trimesh_gaussian)

if (VCG_HEADER_ONLY)
	set(SOURCES
		trimesh_gaussian.cpp
		${VCG_INCLUDE_DIRS}/wrap/ply/plylib.cpp)
endif()

add_executable(trimesh_gaussian
	${SOURCES})

target_link_libraries(
	trimesh_gaussian
	PUBLIC
		vcglib
	)
