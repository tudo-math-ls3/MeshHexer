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

target_compile_options(HexMesher PRIVATE -Wall)

find_package(CGAL)
find_package(pmp)
target_link_libraries(HexMesher PRIVATE CGAL::CGAL pmp::pmp)
