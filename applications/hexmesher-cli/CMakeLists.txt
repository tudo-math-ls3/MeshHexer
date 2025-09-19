set(HEXMESHER_CLI_SOURCES
  hexmesher-cli.cpp
)

set(HEXMESHER_CLI_HEADERS
  timing.hpp
)

add_executable(hexmesher-cli ${HEXMESHER_CLI_SOURCES})
target_sources(hexmesher-cli
  PUBLIC FILE_SET headers
  TYPE HEADERS
  FILES ${HEXMESHER_CLI_HEADERS})
target_link_libraries(hexmesher-cli PRIVATE HexMesher)
