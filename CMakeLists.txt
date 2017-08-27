cmake_minimum_required(VERSION 2.6)
project(directPsi CXX Fortran)


MESSAGE("Hello, world")
IF(DEFINED APPLE)
    MESSAGE("APPLE")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
ENDIF(DEFINED APPLE)

set(CMAKE_OSX_DEPLOYMENT_TARGET "")

include(FindROOT.cmake)
MESSAGE(${ROOT_LIBRARIES})
include(FindEvtGen.cmake)

include_directories(src)
include_directories(${ROOT_INCLUDE_DIR})
include_directories(${EVTGEN_INCLUDE_DIR})

## set fortran flags
INCLUDE(SetCompileFlag.cmake)
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                 Fortran "-fno-automatic")



include_directories(/Users/luchinsky/Work/clhep/install/include/)

include_directories(/Users/luchinsky/Work/clhep/install/include/)
include_directories(/Users/luchinsky/Work/EvtGen/external/HepMC/include/)

# adding sources for EvtGen
include_directories(src/)
file(GLOB_RECURSE evtgen_sources src/*.cc src/*.for src/src/*)
add_executable(chiC_all.exe ${evtgen_sources})
target_link_libraries(chiC_all.exe ${ROOT_LIBRARIES} ${EVTGEN_LIBRARY_DIR}/archive/libEvtGen.a)
target_link_libraries(chiC_all.exe /Users/luchinsky/Work/EvtGen/external/HepMC/lib/libHepMC.a)
