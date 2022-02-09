target_sources(mrchem PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/QMFunction.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Orbital.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/OrbitalIterator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Density.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/qmfunction_utils.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/orbital_utils.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/density_utils.cpp
    )

get_filename_component(_dirname ${CMAKE_CURRENT_LIST_DIR} NAME)

list(APPEND ${_dirname}_h
    ${CMAKE_CURRENT_SOURCE_DIR}/ComplexFunction.h
    ${CMAKE_CURRENT_SOURCE_DIR}/Density.h
    ${CMAKE_CURRENT_SOURCE_DIR}/density_utils.h
    ${CMAKE_CURRENT_SOURCE_DIR}/Orbital.h
    ${CMAKE_CURRENT_SOURCE_DIR}/OrbitalIterator.h
    ${CMAKE_CURRENT_SOURCE_DIR}/orbital_utils.h
    ${CMAKE_CURRENT_SOURCE_DIR}/QMFunction.h
    ${CMAKE_CURRENT_SOURCE_DIR}/qmfunction_utils.h
    )

# Sets install directory for all the headers in the list
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/${_dirname})
foreach(_header IN LISTS ${_dirname}_h)
  configure_file(
    ${_header}
    ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/${_dirname}
    COPYONLY
    )
  install(
    FILES ${_header}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/qmfunctions
    )
endforeach()

