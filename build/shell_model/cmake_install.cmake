# Install script for directory: /home/ubuntu/mechanics_of_epithelial_domes/shell_model

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/ubuntu/mechanics_of_epithelial_domes/source_compiled")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/home/ubuntu/mechanics_of_epithelial_domes/source_compiled/bin/hlshell_model" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/ubuntu/mechanics_of_epithelial_domes/source_compiled/bin/hlshell_model")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/home/ubuntu/mechanics_of_epithelial_domes/source_compiled/bin/hlshell_model"
         RPATH "/home/ubuntu/local/hiperlife/lib:/home/ubuntu/local/hiperlife/third-party/Mumps/lib:/usr/lib/x86_64-linux-gnu/lapack:/usr/lib/x86_64-linux-gnu/blas:/usr/lib/x86_64-linux-gnu/openmpi/lib:/home/ubuntu/local/hiperlife/third-party/Gmsh/lib:/home/ubuntu/local/hiperlife/third-party/Trilinos/lib:/home/ubuntu/local/hiperlife/third-party/vtk/lib")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/ubuntu/mechanics_of_epithelial_domes/source_compiled/bin/hlshell_model")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/ubuntu/mechanics_of_epithelial_domes/source_compiled/bin" TYPE EXECUTABLE FILES "/home/ubuntu/mechanics_of_epithelial_domes/build/shell_model/hlshell_model")
  if(EXISTS "$ENV{DESTDIR}/home/ubuntu/mechanics_of_epithelial_domes/source_compiled/bin/hlshell_model" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/ubuntu/mechanics_of_epithelial_domes/source_compiled/bin/hlshell_model")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/home/ubuntu/mechanics_of_epithelial_domes/source_compiled/bin/hlshell_model"
         OLD_RPATH "/home/ubuntu/local/hiperlife/lib:/home/ubuntu/local/hiperlife/third-party/Mumps/lib:/usr/lib/x86_64-linux-gnu/lapack:/usr/lib/x86_64-linux-gnu/blas:/usr/lib/x86_64-linux-gnu/openmpi/lib:/home/ubuntu/local/hiperlife/third-party/Gmsh/lib:/home/ubuntu/local/hiperlife/third-party/Trilinos/lib:/home/ubuntu/local/hiperlife/third-party/vtk/lib:"
         NEW_RPATH "/home/ubuntu/local/hiperlife/lib:/home/ubuntu/local/hiperlife/third-party/Mumps/lib:/usr/lib/x86_64-linux-gnu/lapack:/usr/lib/x86_64-linux-gnu/blas:/usr/lib/x86_64-linux-gnu/openmpi/lib:/home/ubuntu/local/hiperlife/third-party/Gmsh/lib:/home/ubuntu/local/hiperlife/third-party/Trilinos/lib:/home/ubuntu/local/hiperlife/third-party/vtk/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/home/ubuntu/mechanics_of_epithelial_domes/source_compiled/bin/hlshell_model")
    endif()
  endif()
endif()

