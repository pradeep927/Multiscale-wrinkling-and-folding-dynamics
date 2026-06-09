# Install script for directory: /home/ubuntu/mechanics_of_epithelial_domes/vertex_model

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/ubuntu/shell_models_continuum/source_compiled")
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
  if(EXISTS "$ENV{DESTDIR}/home/ubuntu/shell_models_continuum/source_compiled/bin/hlvertex_model" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/ubuntu/shell_models_continuum/source_compiled/bin/hlvertex_model")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/home/ubuntu/shell_models_continuum/source_compiled/bin/hlvertex_model"
         RPATH "hlCore:hlUtils:hlAdjcyLists:hlFields:hlBasisFunctions:hlMeshCreator:hlMeshMap:hlCubature:hlCubatureSet:hlBasisFunctionsSet:hlDistMesh:hlIO:hlDOFsHandler:hlIntegration:hlHiPerProblem:hlLinearSolvers:hlNonlinearSolvers:hlModelUtils:hlPostProc:hlRemesher:hlTimeIterators:hlPredefinedElementFillings:Amesos2::all_libs:ShyLU_Node::all_libs:ShyLU_NodeTacho::all_libs:Belos::all_libs:ML::all_libs:Ifpack::all_libs:Amesos::all_libs:Galeri::all_libs:AztecOO::all_libs:Isorropia::all_libs:Xpetra::all_libs:Thyra::all_libs:ThyraTpetraAdapters::all_libs:ThyraEpetraExtAdapters::all_libs:ThyraEpetraAdapters::all_libs:ThyraCore::all_libs:TrilinosSS::all_libs:Tpetra::all_libs:TpetraCore::all_libs:TpetraTSQR::all_libs:EpetraExt::all_libs:Triutils::all_libs:Zoltan::all_libs:Epetra::all_libs:RTOp::all_libs:KokkosKernels::all_libs:Teuchos::all_libs:TeuchosKokkosComm::all_libs:TeuchosKokkosCompat::all_libs:TeuchosRemainder::all_libs:TeuchosNumerics::all_libs:TeuchosComm::all_libs:TeuchosParameterList::all_libs:TeuchosParser::all_libs:TeuchosCore::all_libs:Kokkos::all_libs:/home/ubuntu/local/hiperlife/lib:/home/ubuntu/local/hiperlife/third-party/Trilinos/lib:/home/ubuntu/local/hiperlife/third-party/Gmsh/lib:/home/ubuntu/local/hiperlife/third-party/Mumps/lib:/usr/lib/x86_64-linux-gnu/lapack:/usr/lib/x86_64-linux-gnu/blas:/usr/lib/x86_64-linux-gnu/openmpi/lib:/home/ubuntu/local/hiperlife/third-party/vtk/lib")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/ubuntu/shell_models_continuum/source_compiled/bin/hlvertex_model")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/ubuntu/shell_models_continuum/source_compiled/bin" TYPE EXECUTABLE FILES "/home/ubuntu/mechanics_of_epithelial_domes/build/vertex_model/hlvertex_model")
  if(EXISTS "$ENV{DESTDIR}/home/ubuntu/shell_models_continuum/source_compiled/bin/hlvertex_model" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/ubuntu/shell_models_continuum/source_compiled/bin/hlvertex_model")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/home/ubuntu/shell_models_continuum/source_compiled/bin/hlvertex_model"
         OLD_RPATH "/home/ubuntu/local/hiperlife/lib:/home/ubuntu/local/hiperlife/third-party/Trilinos/lib:/home/ubuntu/local/hiperlife/third-party/Gmsh/lib:/home/ubuntu/local/hiperlife/third-party/Mumps/lib:/usr/lib/x86_64-linux-gnu/lapack:/usr/lib/x86_64-linux-gnu/blas:/usr/lib/x86_64-linux-gnu/openmpi/lib:/home/ubuntu/local/hiperlife/third-party/vtk/lib::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
         NEW_RPATH "hlCore:hlUtils:hlAdjcyLists:hlFields:hlBasisFunctions:hlMeshCreator:hlMeshMap:hlCubature:hlCubatureSet:hlBasisFunctionsSet:hlDistMesh:hlIO:hlDOFsHandler:hlIntegration:hlHiPerProblem:hlLinearSolvers:hlNonlinearSolvers:hlModelUtils:hlPostProc:hlRemesher:hlTimeIterators:hlPredefinedElementFillings:Amesos2::all_libs:ShyLU_Node::all_libs:ShyLU_NodeTacho::all_libs:Belos::all_libs:ML::all_libs:Ifpack::all_libs:Amesos::all_libs:Galeri::all_libs:AztecOO::all_libs:Isorropia::all_libs:Xpetra::all_libs:Thyra::all_libs:ThyraTpetraAdapters::all_libs:ThyraEpetraExtAdapters::all_libs:ThyraEpetraAdapters::all_libs:ThyraCore::all_libs:TrilinosSS::all_libs:Tpetra::all_libs:TpetraCore::all_libs:TpetraTSQR::all_libs:EpetraExt::all_libs:Triutils::all_libs:Zoltan::all_libs:Epetra::all_libs:RTOp::all_libs:KokkosKernels::all_libs:Teuchos::all_libs:TeuchosKokkosComm::all_libs:TeuchosKokkosCompat::all_libs:TeuchosRemainder::all_libs:TeuchosNumerics::all_libs:TeuchosComm::all_libs:TeuchosParameterList::all_libs:TeuchosParser::all_libs:TeuchosCore::all_libs:Kokkos::all_libs:/home/ubuntu/local/hiperlife/lib:/home/ubuntu/local/hiperlife/third-party/Trilinos/lib:/home/ubuntu/local/hiperlife/third-party/Gmsh/lib:/home/ubuntu/local/hiperlife/third-party/Mumps/lib:/usr/lib/x86_64-linux-gnu/lapack:/usr/lib/x86_64-linux-gnu/blas:/usr/lib/x86_64-linux-gnu/openmpi/lib:/home/ubuntu/local/hiperlife/third-party/vtk/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/home/ubuntu/shell_models_continuum/source_compiled/bin/hlvertex_model")
    endif()
  endif()
endif()

