set(HEXMESHER_SOURCES
  hexmesher.cpp
  io.cpp
)

set(HEXMESHER_HEADERS
  hexmesher.hpp
  types.hpp
  io.hpp
)

add_library(HexMesher)

target_link_libraries(HexMesher PRIVATE hexmesher-settings)
target_sources(HexMesher PRIVATE ${HEXMESHER_SOURCES})
target_sources(HexMesher PUBLIC FILE_SET HEADERS FILES ${HEXMESHER_HEADERS})

find_package(cgal)
target_link_libraries(HexMesher PRIVATE CGAL::CGAL Boost::boost)
