target_sources(mrchem PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/DeltaOperator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/IdentityOperator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/NuclearOperator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/NuclearGradientOperator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/PositionOperator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/DistanceOperator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/QMDerivative.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/QMPotential.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/QMSpin.cpp
    )
