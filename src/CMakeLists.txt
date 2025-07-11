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

add_library(HexMesher)

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

find_package(CGAL)
find_package(pmp)
target_link_libraries(HexMesher PRIVATE CGAL::CGAL pmp::pmp)

if(HEXMESHER_TESTING)
  find_package(Catch2)

  include(Catch)
  foreach (test ${HEXMESHER_TESTS})
    add_executable(${test} ${test}.cpp)

    target_link_libraries(${test} PRIVATE HexMesher Catch2::Catch2WithMain)
    target_include_directories(${test} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

    add_dependencies(tests ${test})

    catch_discover_tests(${test})
  endforeach()
endif()
