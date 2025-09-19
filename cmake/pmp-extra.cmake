if(NOT TARGET pmp::pmp)
  FetchContent_GetProperties(pmp)

  file(GLOB pmp-sources ${pmp_SOURCE_DIR}/src/pmp/*.cpp ${pmp_SOURCE_DIR}/src/pmp/algorithms/*.cpp ${pmp_SOURCE_DIR}/src/pmp/io/*.cpp)
  file(GLOB pmp-headers ${pmp_SOURCE_DIR}/src/pmp/*.h ${pmp_SOURCE_DIR}/src/pmp/algorithms/*.h ${pmp_SOURCE_DIR}/src/pmp/io/*.h)
  file(GLOB_RECURSE eigen-headers ${pmp_SOURCE_DIR}/external/eigen-3.4.0/Eigen/*)

  if(HEXMESHER_SHARED_LIB)
    add_library(hexmesher-extern-pmp SHARED)
  else()
    add_library(hexmesher-extern-pmp STATIC)
  endif()

  target_sources(hexmesher-extern-pmp PRIVATE ${pmp-sources})

  # NOTE(mmuegge): pmp includes the Eigen 3.4.0 headers, but does not link
  # against Eigen. This is intended.  pmp uses the type definitions of eigen to
  # support converting from Eigen vectors/matrices to pmp vectors/matrices, but
  # does not use the implementation of eigen in any way.
  target_sources(hexmesher-extern-pmp
    PUBLIC FILE_SET pmp_headers
    TYPE HEADERS
    BASE_DIRS ${pmp_SOURCE_DIR}/src/ ${pmp_SOURCE_DIR}/external/eigen-3.4.0/
    FILES ${pmp-headers} ${eigen-headers})


  target_compile_features(hexmesher-extern-pmp PRIVATE cxx_std_20)
  target_compile_options(hexmesher-extern-pmp PRIVATE -w)

  if(HEXMESHER_HAVE_OMP)
    find_package(OpenMP REQUIRED)
    if(OpenMP_CXX_FOUND)
      target_link_libraries(hexmesher-extern-pmp PUBLIC OpenMP::OpenMP_CXX)
    endif()
  endif()

  add_library(pmp::pmp ALIAS hexmesher-extern-pmp)

  install(
    TARGETS hexmesher-extern-pmp
    EXPORT hexmesherv2Targets
    FILE_SET pmp_headers
  )
endif()
