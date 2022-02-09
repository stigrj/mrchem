list(APPEND api_h
  ${CMAKE_CURRENT_LIST_DIR}/QMFunctions
  )

# Sets install directory for all the headers in the list
foreach(_header IN LISTS api_h)
  configure_file(
    ${_header}
    ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
    COPYONLY
    )
endforeach()

# Append to mrchem_headers list.
# These will be installed as the public headers of the target
list(APPEND mrchem_headers "${api_h}")
set(mrchem_headers "${mrchem_headers}" PARENT_SCOPE)

