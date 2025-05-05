set(HEXMESHER_SOURCES
  hexmesher.cpp
)

set(HEXMESHER_HEADERS
  hexmesher.hpp
)

add_library(HexMesher)
target_link_libraries(HexMesher PRIVATE hexmesher-settings)
target_sources(HexMesher PRIVATE ${HEXMESHER_SOURCES})
target_sources(HexMesher PUBLIC FILE_SET HEADERS FILES ${HEXMESHER_HEADERS})

