#============================================================
#======================== Initialize ========================
#============================================================
cmake_minimum_required(VERSION 3.17.0)

INCLUDE(${CMAKE_SOURCE_DIR}/userConfig.cmake)
project(${PROJECT_NAME} VERSION 5.0)

#============================================================
#===================== Options and flags ====================
#============================================================
# Options
option(BUILD_SHARED_LIBS "Build hiperlife with shared libraries." ON)

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${HL_BASE_PATH}/apps" CACHE PATH "Install directory" FORCE)
  message(STATUS "Setting default install path at: ${CMAKE_INSTALL_PREFIX}")
  message(STATUS "To change the install path, add the option -DCMAKE_INSTALL_PREFIX=/path/")
else()
  message(STATUS "Install path set at: ${CMAKE_INSTALL_PREFIX}")
endif()

if(CMAKE_BUILD_TYPE STREQUAL "")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
  message(STATUS "Hiperlife build type not set, default build type: ${CMAKE_BUILD_TYPE}")
else()
  message(STATUS "Hiperlife build type set to: ${CMAKE_BUILD_TYPE}")
endif()

#============================================================
#============ HIPERLIFE AND DEPENDENCIES ====================
#============================================================
message(STATUS "Third party libraries:")

find_package(hiperlife QUIET HINTS ${HL_BASE_PATH})
if(hiperlife_FOUND)
  message(STATUS "  hiperlife:\t Found version ${hiperlife_VERSION} at ${hiperlife_DIR}")
  message(STATUS "  hiperlife include dirs: ${hiperlife_INCLUDE_DIR}")
  message(STATUS "  hiperlife libs: ${hiperlife_LIBRARIES}")
else()
  message(FATAL_ERROR "hiperlife not found at '${HL_BASE_PATH}'. Please add a valid hiperlife installation with the cmake option -DHL_BASE_PATH=/path/to/hiperlife")
endif(hiperlife_FOUND)

if(hiperlife_VERSION VERSION_GREATER ${PROJECT_NAME}_VERSION)
  message(FATAL_ERROR "The version of hiperlife (${hiperlife_VERSION}) is greater than this project version (${hl_BaseProject_VERSION}). Please update")
elseif(NOT hiperlife_VERSION VERSION_EQUAL ${PROJECT_NAME}_VERSION)
  message(WARNING "hiperlife version (${hiperlife_VERSION}) and this project version (${${PROJECT_NAME}_VERSION}) do not match. This may cause issues. Please consider updating")
endif()

#============================================================
#=================== PROJECT SETTINGS =======================
#============================================================
set(PROJECT_INSTALL_PATH ${CMAKE_INSTALL_PREFIX})
include_directories(${hiperlife_INCLUDE_DIR})

#Library directories
set(APPS_LIBS_INSTALL ${CMAKE_INSTALL_PREFIX}/lib)
set(APPS_INCLUDE_INSTALL ${CMAKE_INSTALL_PREFIX}/include)
set(APPS_LIBRARIES CACHE INTERNAL "APPS_LIBRARIES")

#============================================================
#======================== LIBRARIES =========================
#============================================================

set(APPS_INSTALL_DIR ${CMAKE_CURRENT_LIST_DIR}/bin/)
#Add (or modify) applications that do not depend on vtk and gmsh here:
#add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/shell_model)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/active_gel_bilayer_shell)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/shell_model)
#Add applications that depend on vtk and/or gmsh here: viscous_inext
#if(VTK_FOUND AND GMSH_LIB)

#elseif(GMSH_LIB)
  #add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/MyApp_with_GMSH)
#endif()
