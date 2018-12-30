cmake_minimum_required(VERSION 2.6)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(/usr/local/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../kiss/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../HepPID/include)

add_library(crpropa STATIC
	src/base64.cpp
    src/Candidate.cpp
    src/Common.cpp
    src/GridTools.cpp
    src/Module.cpp
    src/ParticleID.cpp
    src/ParticleMass.cpp
    src/ParticleState.cpp
    src/Random.cpp
    src/Variant.cpp
    src/magneticField/MagneticField.cpp
    src/magneticField/MagneticFieldGrid.cpp
    src/module/Boundary.cpp
    src/module/PropagationCK.cpp
)

SET_TARGET_PROPERTIES(crpropa PROPERTIES COMPILE_FLAGS -O3)
