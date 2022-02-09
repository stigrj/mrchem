# <<< Build MRChem library >>>
add_library(mrchem
    driver.cpp
    parallel.cpp
    mrenv.cpp
  )

add_dependencies(mrchem mrchem-info)

target_compile_definitions(mrchem
  PUBLIC
    $<$<AND:$<TARGET_EXISTS:OpenMP::OpenMP_CXX>,$<BOOL:${ENABLE_OPENMP}>>:MRCHEM_HAS_OMP>
    $<$<AND:$<TARGET_EXISTS:MPI::MPI_CXX>,$<BOOL:${ENABLE_MPI}>>:MRCHEM_HAS_MPI>
  )

set_target_properties(mrchem PROPERTIES EXPORT_PROPERTIES "MRCHEM_HAS_OMP;MRCHEM_HAS_MPI")

target_include_directories(mrchem
  PUBLIC
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  )

list(APPEND _public_headers
  ${CMAKE_CURRENT_SOURCE_DIR}/driver.h
  ${CMAKE_CURRENT_SOURCE_DIR}/mrchem.h
  ${CMAKE_CURRENT_SOURCE_DIR}/mrenv.h
  ${CMAKE_CURRENT_SOURCE_DIR}/parallel.h
  ${PROJECT_BINARY_DIR}/config.h
  ${PROJECT_BINARY_DIR}/version.h
  )

# RPATH fixing
file(RELATIVE_PATH _rel ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR} ${CMAKE_INSTALL_PREFIX})
if(APPLE)
  set(_rpath "@loader_path/${_rel}")
else()
  set(_rpath "\$ORIGIN/${_rel}")
endif()
file(TO_NATIVE_PATH "${_rpath}/${CMAKE_INSTALL_LIBDIR}" MRChem_RPATH)

add_subdirectory(analyticfunctions)
add_subdirectory(chemistry)
add_subdirectory(environment)
add_subdirectory(mrdft)
add_subdirectory(qmfunctions)
add_subdirectory(qmoperators)
add_subdirectory(scf_solver)
add_subdirectory(tensor)
add_subdirectory(utils)

target_link_libraries(mrchem
  PRIVATE
    Eigen3::Eigen
  PUBLIC
    XCFun::xcfun
    MRCPP::mrcpp
    nlohmann_json::nlohmann_json
    $<$<AND:$<TARGET_EXISTS:OpenMP::OpenMP_CXX>,$<BOOL:${ENABLE_OPENMP}>>:OpenMP::OpenMP_CXX>
    $<$<AND:$<TARGET_EXISTS:MPI::MPI_CXX>,$<BOOL:${ENABLE_MPI}>>:MPI::MPI_CXX>
  )

set_target_properties(mrchem
  PROPERTIES
    VERSION ${MRCHEM_VERSION_MAJOR}
    SOVERSION ${MRCHEM_VERSION_MAJOR}
    MACOSX_RPATH ON
    OUTPUT_NAME "mrchem"
    EXPORT_NAME "mrchem"
    PUBLIC_HEADER "${_public_headers}"
  )

if(APPLE)
  set_target_properties(mrchem
    PROPERTIES
      LINK_FLAGS "-undefined dynamic_lookup"
    )
endif()

install(
  TARGETS
    mrchem
  EXPORT
    "MRChemTargets"
  RUNTIME DESTINATION
    ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION
    ${CMAKE_INSTALL_LIBDIR}
  PUBLIC_HEADER DESTINATION
    ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
  )

# <<<  Export Config  >>>
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  ${PROJECT_BINARY_DIR}/${CMAKECONFIG_INSTALL_DIR}/MRChemConfigVersion.cmake
  VERSION ${MRCHEM_VERSION}
  COMPATIBILITY SameMajorVersion
  )

configure_package_config_file(
  ${PROJECT_SOURCE_DIR}/cmake/MRChemConfig.cmake.in
  ${PROJECT_BINARY_DIR}/${CMAKECONFIG_INSTALL_DIR}/MRChemConfig.cmake
  INSTALL_DESTINATION
    ${CMAKECONFIG_INSTALL_DIR}
  )

install(
  FILES
    ${PROJECT_BINARY_DIR}/${CMAKECONFIG_INSTALL_DIR}/MRChemConfig.cmake
    ${PROJECT_BINARY_DIR}/${CMAKECONFIG_INSTALL_DIR}/MRChemConfigVersion.cmake
  DESTINATION
    ${CMAKECONFIG_INSTALL_DIR}
  )

install(
  EXPORT
    "MRChemTargets"
  NAMESPACE
    "MRChem::"
  DESTINATION
    ${CMAKECONFIG_INSTALL_DIR}
    )

# Alias useful for use with FetchContent
add_library(MRChem::mrchem ALIAS mrchem)

