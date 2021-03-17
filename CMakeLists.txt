# CMakeLists.txt for event package. It creates a library with dictionary and a main program
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(helloworld)

# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED)
find_package(Ceres REQUIRED)

include(${ROOT_USE_FILE})
include_directories(${CERES_INCLUDE_DIRS})
include_directories(${ROOT_INCLUDE_DIRS})

#and here
#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)

#set(root_lib libCore.so libPhysics.so libTree.so libRIO.so )
#ROOT_GENERATE_DICTIONARY(G__helloworld helloworld.h LINKDEF helloworldLinkDef.h)


#---Create a shared library with geneated dictionary
#add_library(helloworld SHARED helloworld.cxx G__helloworld.cxx)
add_executable(helloworld helloworld.cxx ${OBJS})
target_link_libraries(helloworld ${CERES_LIBRARIES} ${ROOT_LIBRARIES})

#---Create  a main program using the library


