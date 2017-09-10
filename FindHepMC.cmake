################################################################################
 #    Copyright (C) 2014 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH    #
 #                                                                              #
 #              This software is distributed under the terms of the             #
 #         GNU Lesser General Public Licence version 3 (LGPL) version 3,        #
 #                  copied verbatim in the file "LICENSE"                       #
 ################################################################################
# - Try to find HepMC instalation
# Once done this will define
#

MESSAGE(STATUS "Looking for HepMC ...")

FIND_PATH(HepMC_INCLUDE_DIR NAMES HepMC/Version.h PATHS
  /Users/luchinsky/Work/EvtGen/external/HepMC/include/
  /afs/ihep.su/user/a/aluchins/local/EvtGen/external/HepMC/include/
  NO_DEFAULT_PATH
)


FIND_PATH(HepMC_LIBRARY_DIR NAMES libHepMC.a PATHS
  /Users/luchinsky/Work/EvtGen/external/HepMC/lib/
 /afs/ihep.su/user/a/aluchins/local/EvtGen/external/HepMC/lib/
  NO_DEFAULT_PATH
)


if (HepMC_INCLUDE_DIR AND HepMC_LIBRARY_DIR)
   set(HepMC_FOUND TRUE)
endif (HepMC_INCLUDE_DIR AND HepMC_LIBRARY_DIR)

if (HepMC_FOUND)
  if (NOT HepMC_FOUND_QUIETLY)
    MESSAGE(STATUS "Looking for HepMC... - found ${HepMC_LIBRARY_DIR}")
    SET(LD_LIBRARY_PATH ${LD_LIBRARY_PATH} ${HepMC_LIBRARY_DIR})
  endif (NOT HepMC_FOUND_QUIETLY)
else (HepMC_FOUND)
  if (HepMC_FOUND_REQUIRED)
    message(FATAL_ERROR "Looking for HepMC... - Not found")
  endif (HepMC_FOUND_REQUIRED)
endif (HepMC_FOUND)