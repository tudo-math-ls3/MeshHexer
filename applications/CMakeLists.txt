add_subdirectory(hexmesher-cli)

find_package(cgal CONFIG)

add_executable(dbg-min-gap dbg_min_gap.cpp)
target_link_libraries(dbg-min-gap PRIVATE HexMesher CGAL::CGAL Boost::boost)
