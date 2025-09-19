set(HEXMESHER_SOURCES
  hexmesher.cpp
  io.cpp
  meshing.cpp
  properties.cpp
)

set(HEXMESHER_PUBLIC_HEADERS
  hexmesher.hpp
  types.hpp
)

set(HEXMESHER_PRIVATE_HEADERS
  types.hpp
  io.hpp
  meshing.hpp
  properties.hpp
  warnings.hpp
)

set(HEXMESHER_TESTS
  properties-test
  io-test
)

if(HEXMESHER_SHARED_LIB)
  add_library(HexMesher SHARED)
else()
  add_library(HexMesher STATIC)
endif()

target_link_libraries(HexMesher PRIVATE hexmesher-settings)
target_sources(HexMesher PRIVATE ${HEXMESHER_SOURCES})
target_sources(HexMesher
  PUBLIC FILE_SET public_headers
  TYPE HEADERS
  FILES ${HEXMESHER_PUBLIC_HEADERS})
target_sources(HexMesher
  PRIVATE FILE_SET private_headers
  TYPE HEADERS
  FILES ${HEXMESHER_PRIVATE_HEADERS})

find_package(CGAL 6.0.1)
find_package(pmp 3.0.0)

if(HEXMESHER_SHARED_LIB)
  target_link_libraries(HexMesher PRIVATE CGAL::CGAL pmp::pmp)
else()
  # Per default CMake exports the private dependencies of static library, as they might contains symbols
  # that are needed by the main library.
  # We know that CGAL is a header only library and thus supress that behaviour here
  # See: https://gitlab.kitware.com/cmake/cmake/-/issues/15415
  target_link_libraries(HexMesher PRIVATE $<BUILD_LOCAL_INTERFACE:CGAL::CGAL> pmp::pmp)
endif()


if(HEXMESHER_TESTING)
  find_package(Catch2)

  include(Catch)
  foreach (test ${HEXMESHER_TESTS})
    add_executable(${test} ${test}.cpp)

    target_link_libraries(${test} PRIVATE HexMesher Catch2::Catch2WithMain CGAL::CGAL pmp::pmp)
    target_include_directories(${test} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

    add_dependencies(tests ${test})

    catch_discover_tests(${test})
  endforeach()
endif()
