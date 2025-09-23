set(MESHHEXER_SOURCES
  meshhexer.cpp
  io.cpp
  meshing.cpp
  properties.cpp
)

set(MESHHEXER_PUBLIC_HEADERS
  meshhexer.hpp
  types.hpp
)

set(MESHHEXER_PRIVATE_HEADERS
  types.hpp
  io.hpp
  meshing.hpp
  properties.hpp
  warnings.hpp
)

set(MESHHEXER_TESTS
  properties-test
  io-test
)

if(MESHHEXER_SHARED_LIB)
  add_library(MeshHexer SHARED)
else()
  add_library(MeshHexer STATIC)
endif()

target_link_libraries(MeshHexer PRIVATE meshhexer-settings)
target_sources(MeshHexer PRIVATE ${MESHHEXER_SOURCES})
target_sources(MeshHexer
  PUBLIC FILE_SET public_headers
  TYPE HEADERS
  FILES ${MESHHEXER_PUBLIC_HEADERS})
target_sources(MeshHexer
  PRIVATE FILE_SET private_headers
  TYPE HEADERS
  FILES ${MESHHEXER_PRIVATE_HEADERS})

find_package(CGAL 6.0.1)
find_package(pmp 3.0.0)

if(MESHHEXER_SHARED_LIB)
  target_link_libraries(MeshHexer PRIVATE CGAL::CGAL pmp::pmp)
else()
  # Per default CMake exports the private dependencies of static library, as they might contains symbols
  # that are needed by the main library.
  # We know that CGAL is a header only library and thus supress that behaviour here
  # See: https://gitlab.kitware.com/cmake/cmake/-/issues/15415
  target_link_libraries(MeshHexer PRIVATE $<BUILD_LOCAL_INTERFACE:CGAL::CGAL> pmp::pmp)
endif()


if(MESHHEXER_TESTING)
  find_package(Catch2)

  include(Catch)
  foreach (test ${MESHHEXER_TESTS})
    add_executable(${test} ${test}.cpp)

    target_link_libraries(${test} PRIVATE MeshHexer Catch2::Catch2WithMain CGAL::CGAL pmp::pmp)
    target_include_directories(${test} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

    add_dependencies(tests ${test})

    catch_discover_tests(${test})
  endforeach()
endif()

if(MESHHEXER_DOCUMENTATION)
  find_package(Doxygen)
  set(DOXYGEN_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/doc)
  set(DOXYGEN_EXTRACT_ALL YES)
  if(MESHHEXER_DOCUMENTATION_INTERNAL)
    doxygen_add_docs(doc ${MESHHEXER_PUBLIC_HEADERS} ${MESHHEXER_PRIVATE_HEADERS})
  else()
    doxygen_add_docs(doc ${MESHHEXER_PUBLIC_HEADERS})
  endif()
endif()
