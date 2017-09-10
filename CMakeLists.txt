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

include(FindHepMC.cmake)
MESSAGE(${HepMC_LIBRARY_DIR})


include_directories(src)
include_directories(${ROOT_INCLUDE_DIR})
include_directories(${EVTGEN_INCLUDE_DIR})
include_directories(${HepMC_INCLUDE_DIR})

## set fortran flags
INCLUDE(SetCompileFlag.cmake)
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                 Fortran "-fno-automatic")




# adding sources for EvtGen
include_directories(src/)
file(GLOB_RECURSE evtgen_sources src/chic_all.cc src/*.for src/src/*)
add_executable(chiC_all.exe ${evtgen_sources})
target_link_libraries(chiC_all.exe ${ROOT_LIBRARIES} ${EVTGEN_LIBRARY_DIR}/archive/libEvtGen.a)
target_link_libraries(chiC_all.exe ${HepMC_LIBRARY_DIR}/libHepMC.a)

# adding sources for read_data
add_executable(read_data.exe src/read_data.cc)
target_link_libraries(read_data.exe ${ROOT_LIBRARIES} ${EVTGEN_LIBRARY_DIR}/archive/libEvtGen.a)
target_link_libraries(read_data.exe ${HepMC_LIBRARY_DIR}/libHepMC.a)
