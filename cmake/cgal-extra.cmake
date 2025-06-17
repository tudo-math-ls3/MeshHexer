if(NOT TARGET CGAL::CGAL)
  FetchContent_GetProperties(cgal)
  find_package(Boost 1.81 REQUIRED)

  add_library(hexmesher-cgal-extern INTERFACE IMPORTED)
  target_include_directories(hexmesher-cgal-extern INTERFACE "${cgal_SOURCE_DIR}/include")
  target_link_libraries(hexmesher-cgal-extern INTERFACE Boost::boost)
  target_compile_definitions(hexmesher-cgal-extern INTERFACE CGAL_DISABLE_GMP)

  add_library(CGAL::CGAL ALIAS hexmesher-cgal-extern)
else()
  message(STATUS "Target CGAL::CGAL already exists!")
endif()
