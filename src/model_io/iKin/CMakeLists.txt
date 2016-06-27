# Copyright: (C) 2016 Fondazione Istituto Italiano di Tecnologia
# Authors: Silvio Traversaro
# CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT

project(iDynTree_ModelIO_iKin CXX)

set(IDYNTREE_MODELIO_IKIN_HEADERS include/iDynTree/iKinModelImport.h)
set(IDYNTREE_MODELIO_IKIN_SOURCES src/iKinModelImport.cpp)

source_group("Source Files" FILES ${IDYNTREE_MODELIO_IKIN_SOURCES})
source_group("Header Files" FILES ${IDYNTREE_MODELIO_IKIN_HEADERS})

# share headers with all iDynTree targets
set(libraryname idyntree-modelio-ikin)

add_library(${libraryname} ${IDYNTREE_MODELIO_IKIN_SOURCES} ${IDYNTREE_MODELIO_IKIN_HEADERS})

if(NOT CMAKE_VERSION VERSION_LESS 2.8.12)
    target_include_directories(${libraryname} SYSTEM PUBLIC ${orocos_kdl_INCLUDE_DIRS} ${iKin_INCLUDE_DIRS} ${YARP_INCLUDE_DIRS})
    target_include_directories(${libraryname} SYSTEM PRIVATE ${EIGEN3_INCLUDE_DIR})
else()
    target_include_directories(${libraryname} PUBLIC ${orocos_kdl_INCLUDE_DIRS} ${iKin_INCLUDE_DIRS} ${YARP_INCLUDE_DIRS})
    target_include_directories(${libraryname} PRIVATE ${EIGEN3_INCLUDE_DIR})
endif()
target_include_directories(${libraryname} PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
                                                 "$<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}>"
                                                  ${IDYNTREE_TREE_INCLUDE_DIRS}
                                          PRIVATE ${EIGEN3_INCLUDE_DIR})

target_link_libraries(${libraryname} LINK_PUBLIC idyntree-model idyntree-sensors idyntree-yarp iKin)

target_include_directories(${libraryname} PRIVATE ${EIGEN3_INCLUDE_DIR})

target_include_directories(${libraryname} INTERFACE ${iKin_INCLUDE_DIRS})

include_directories(SYSTEM ${iKin_INCLUDE_DIRS})

if(NOT CMAKE_VERSION VERSION_LESS 2.8.12)
    target_compile_options(${libraryname} PRIVATE ${IDYNTREE_WARNING_FLAGS})
endif()

set_property(TARGET ${libraryname} PROPERTY PUBLIC_HEADER ${IDYNTREE_MODELIO_IKIN_HEADERS})

install(TARGETS ${libraryname}
        EXPORT iDynTree
        COMPONENT runtime
        RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
        LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
        ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT lib
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/iDynTree)

set_property(GLOBAL APPEND PROPERTY ${VARS_PREFIX}_TARGETS ${libraryname})

set_property(GLOBAL APPEND PROPERTY IDYNTREE_TREE_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include)
