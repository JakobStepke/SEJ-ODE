cmake_minimum_required(VERSION 3.20)
project(ASC_ode)

set (CMAKE_CXX_STANDARD 17)


# use linalg from NGSolve, otherwise comment it out
#find_package(NGSolve CONFIG REQUIRED
# HINTS /usr $ENV{NETGENDIR}/.. /tmp/inst/usr /opt/netgen/ /Applications/Netgen.app/Contents/Resources/CMake C:/netgen
#)


include_directories(src)
include_directories(ESJ-bla/src)
add_subdirectory (src)

add_executable (test_ode demos/test_ode.cc)
find_package(LAPACK REQUIRED)
target_link_libraries (test_ode PUBLIC ${LAPACK_LIBRARIES})

add_executable (test_newmark demos/test_newmark.cc)
target_link_libraries (test_newmark PUBLIC ${LAPACK_LIBRARIES})

add_executable (test_alpha demos/test_alpha.cc)
target_link_libraries (test_alpha PUBLIC ${LAPACK_LIBRARIES})


#add_subdirectory (mass_spring)
