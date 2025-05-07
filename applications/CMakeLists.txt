find_package(cgal)

add_executable(hexmesher hexmesher.cpp)
target_link_libraries(hexmesher PRIVATE HexMesher CGAL::CGAL Boost::boost)
