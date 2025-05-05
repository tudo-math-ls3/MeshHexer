if(NOT TARGET CGAL::CGAL)
  FetchContent_GetProperties(cgal)
  find_package(Boost 1.81 REQUIRED)

  add_library(hexmesher-cgal-extern INTERFACE)
  target_include_directories(hexmesher-cgal-extern INTERFACE "${cgal_SOURCE_DIR}/include")
  target_link_libraries(hexmesher-cgal-extern INTERFACE Boost::boost)

  add_library(CGAL::CGAL ALIAS hexmesher-cgal-extern)
endif()