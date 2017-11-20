cmake_minimum_required(VERSION 2.6)
project(directPsi CXX)
enable_language(Fortran)


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



include_directories(src/)

# adding sources for EvtGen
file(GLOB_RECURSE evtgen_sources src/chic_all.cc src/a.for src/src/EvtGenModels/*.cpp)
add_executable(chiC_all.exe ${evtgen_sources})
target_link_libraries(chiC_all.exe ${ROOT_LIBRARIES} ${EVTGEN_LIBRARY_DIR}/archive/libEvtGen.a)
target_link_libraries(chiC_all.exe ${HepMC_LIBRARY_DIR}/libHepMC.a)
SET_TARGET_PROPERTIES(chiC_all.exe  PROPERTIES LINK_FLAGS -lgfortran)

# adding sources for test.exe
file(GLOB_RECURSE test_sources src/test.cc src/a.for src/src/EvtGenModels/*.cpp)
add_executable(test.exe ${test_sources})
target_link_libraries(test.exe ${ROOT_LIBRARIES} ${EVTGEN_LIBRARY_DIR}/archive/libEvtGen.a)
target_link_libraries(test.exe ${HepMC_LIBRARY_DIR}/libHepMC.a)
SET_TARGET_PROPERTIES(test.exe  PROPERTIES LINK_FLAGS -lgfortran)


# adding sources for read_data
add_executable(read_data.exe src/read_data.cc)
target_link_libraries(read_data.exe ${ROOT_LIBRARIES} ${EVTGEN_LIBRARY_DIR}/archive/libEvtGen.a)
target_link_libraries(read_data.exe ${HepMC_LIBRARY_DIR}/libHepMC.a)


# adding sources for MMSh
add_executable(MMSh.exe src/MMSh/MMSh.c++ src/MMSh/rambo_.for src/MMSh/kinematics.c++ src/MMSh/matr2.c++)
target_link_libraries(MMSh.exe ${ROOT_LIBRARIES})
#target_link_libraries(MMSh.exe ${HepMC_LIBRARY_DIR}/libHepMC.a)
SET_TARGET_PROPERTIES(MMSh.exe  PROPERTIES LINK_FLAGS -lgfortran)

# adding sources to gen_pol
add_executable(gen_pol.exe src/matr2_pol/gen_pol.cc)
target_link_libraries(gen_pol.exe ${ROOT_LIBRARIES} ${EVTGEN_LIBRARY_DIR}/archive/libEvtGen.a)
