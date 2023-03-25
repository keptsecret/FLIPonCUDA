cmake_minimum_required(VERSION 3.17)

###########################################################################
# openvdb

set(OPENVDB_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR}/openvdb/ PARENT_SCOPE)

add_subdirectory(openvdb)

set_property(TARGET openvdb_static PROPERTY FOLDER "thirdparty")
