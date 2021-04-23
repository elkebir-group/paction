cmake_minimum_required( VERSION 2.8 )

project( SPRUCE-PACTION )

set( enumerate_src
  src/enum.cpp
  src/stategraph.cpp
  src/character.cpp
  src/charactermatrix.cpp
  src/utils.cpp
  src/matrix.cpp
  src/realmatrix.cpp
  src/tensor.cpp
  src/realtensor.cpp
  src/rootedcladisticancestrygraph.cpp
  src/rootedcladisticnoisyancestrygraph.cpp
  src/perfectphylomatrix.cpp
  src/perfectphylotree.cpp
  src/perfectphylograph.cpp
  src/statetree.cpp
  src/rootedcladisticenumeration.cpp
  src/rootedcladisticnoisyenumeration.cpp
  src/solution.cpp
  src/solutionset.cpp
)

set( enumerate_hdr
  src/stategraph.h
  src/character.h
  src/charactermatrix.h
  src/utils.h
  src/matrix.h
  src/realmatrix.h
  src/tensor.h
  src/realtensor.h
  src/rootedcladisticancestrygraph.h
  src/rootedcladisticnoisyancestrygraph.h
  src/perfectphylomatrix.h
  src/perfectphylotree.h
  src/perfectphylograph.h
  src/statetree.h
  src/rootedcladisticenumeration.h
  src/rootedcladisticnoisyenumeration.h
  src/solution.h
  src/solutionset.h
)

set( LIBLEMON_ROOT "$ENV{HOME}/lemon" CACHE PATH "Additional search directory for lemon library" )

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

if( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-local-typedefs" )
endif()

add_definitions( -DIL_STD )

set( Boost_USE_STATIC_LIBS ON )
find_package( Boost REQUIRED COMPONENTS thread system filesystem )

include_directories( "${LIBLEMON_ROOT}/include" "src" ${Boost_INCLUDE_DIRS} )
link_directories( "${LIBLEMON_ROOT}/lib" )

if(NOT CMAKE_BUILD_TYPE)
set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: Debug Release Profile" FORCE)
endif(NOT CMAKE_BUILD_TYPE)

set( CommonLibs
  m
  emon
  pthread
  ${Boost_LIBRARIES}
)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")

add_executable( enumerate ${enumerate_src} ${enumerate_hdr} )
target_link_libraries( enumerate ${CommonLibs} )

