#------------------------------------------------------------
# CMakeLists for one hiperlife App
#------------------------------------------------------------

#Name
get_filename_component(programname ${CMAKE_CURRENT_LIST_DIR} NAME)
string(REGEX REPLACE "^[0-9]" "" programname ${programname})
set(programname "hl${programname}")

#Info
message(STATUS "programname:  ${programname} ")

#Sources
file(GLOB PROGRAM_CPP *.cpp)
file(GLOB PROGRAM_H *.h)

#Add executable
add_executable(${programname} ${PROGRAM_CPP} ${PROGRAM_H})

#Link libs
target_link_libraries(${programname} ${Trilinos_LIBRARIES})
target_link_libraries(${programname} ${hiperlife_LIBRARIES})
target_link_libraries(${programname} ${APPS_LIBRARIES})
if(VTK_FOUND)
  target_link_libraries(${programname} ${VTK_LIBRARIES})
endif()
if(GMSH_LIB)
  target_link_libraries(${programname} ${GMSH_LIBRARIES})
endif()
if(GINAC_LIBRARIES)
  target_link_libraries(${programname} ${GINAC_LIBRARIES})
endif()

#Install
install(TARGETS ${programname} DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
