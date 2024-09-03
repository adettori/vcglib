cmake_minimum_required(VERSION 3.13)
project(gs_box_filter)

if (VCG_HEADER_ONLY)
	set(SOURCES
		gs_box_filter.cpp
		${VCG_INCLUDE_DIRS}/wrap/ply/plylib.cpp)
endif()

add_executable(gs_box_filter
	${SOURCES})

target_link_libraries(
	gs_box_filter
	PUBLIC
		vcglib
	)
