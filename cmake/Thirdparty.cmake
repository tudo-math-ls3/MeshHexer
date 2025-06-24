# Setup thirdparty dependencies
include(FetchContent)
include(FindPackageOverrideHelper)
include(PrintPackageConfig)

set(MAKE_AVAIL_LIST)
mark_as_advanced(MAKE_AVAIL_LIST)

if(HEXMESHER_PREFER_EXTERNAL_TPL OR Eigen3_DIR)
  if(NOT HEXMESHER_ALLOW_EXTERNAL_DOWNLOAD)
    find_package(Eigen3 3.4.0 CONFIG REQUIRED)
  else()
    find_package(Eigen3 3.4.0 CONFIG)
  endif()
endif()

if(NOT Eigen3_FOUND)
  FetchContent_declare(
    Eigen3
    URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
    URL_HASH MD5=4c527a9171d71a72a9d4186e65bea559
    UPDATE_DISCONNECTED ON
    OVERRIDE_FIND_PACKAGE
    EXCLUDE_FROM_ALL
  )
  list(APPEND MAKE_AVAIL_LIST Eigen3)
  find_package_override_helper(Eigen3 3.4.0 AnyNewerVersion)
else()
  print_package_info(Eigen3)
endif()

if(HEXMESHER_PREFER_EXTERNAL_TPL OR Boost_DIR)
  if(NOT HEXMESHER_ALLOW_EXTERNAL_DOWNLOAD)
    find_package(Boost 1.88 REQUIRED)
  else()
    find_package(Boost 1.88)
  endif()
endif()

if(NOT Boost_FOUND)
  #set(BOOST_INCLUDE_LIBRARIES graph heap logic)

  # Make boost create Boost::boost target
  set(BOOST_ENABLE_COMPATIBILITY_TARGETS ON)
  # Fetch Boost from the github repository. For some reason the github releases
  # contain a root CMakeLists.txt and the website releases do not.
  FetchContent_declare(
    Boost
    URL https://github.com/boostorg/boost/releases/download/boost-1.88.0/boost-1.88.0-cmake.tar.xz
    URL_HASH MD5=3edffaacd2cfe63c240ef1b99497c74f
    UPDATE_DISCONNECTED ON
    # FIND_PACKAGE_ARGS CONFIG GLOBAL
    OVERRIDE_FIND_PACKAGE
    # EXCLUDE_FROM_ALL #TODO: probably a good idea
  )
  list(APPEND MAKE_AVAIL_LIST Boost)
  find_package_override_helper(boost 1.88 AnyNewerVersion)
else()
  print_package_info(Boost)
endif()

if(HEXMESHER_PREFER_EXTERNAL_TPL OR CGAL_DIR)
  if(NOT HEXMESHER_ALLOW_EXTERNAL_DOWNLOAD)
    find_package(CGAL 6.0.1 CONFIG REQUIRED)
  else()
    find_package(CGAL 6.0.1 CONFIG)
  endif()
endif()

message(STATUS "EIGEN3_INCLUDE_DIR: ${EIGEN3_INCLUDE_DIR}")
if(NOT CGAL_FOUND)
  FetchContent_declare(
    CGAL
    URL https://github.com/CGAL/cgal/releases/download/v6.0.1/CGAL-6.0.1.tar.xz
    URL_HASH MD5=944c789615bff14a56d78b398ec2cc49
    SOURCE_SUBDIR Non-Existing # Use non-existing source dir to disable add_subdirectory call of MakeAvailable
    UPDATE_DISCONNECTED ON
    # FIND_PACKAGE_ARGS NAMES CGAL GLOBAL
    OVERRIDE_FIND_PACKAGE
  )

  set(CGAL_SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/_deps/cgal-src)
  set(CGAL_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/_deps/cgal-build)

  list(APPEND MAKE_AVAIL_LIST CGAL)
  find_package_override_helper(cgal 6.0.1 AnyNewerVersion)
else()
  print_package_info(CGAL)
endif()

if(HEXMESHER_PREFER_EXTERNAL_TPL OR pmp_DIR)
  if(NOT HEXMESHER_ALLOW_EXTERNAL_DOWNLOAD)
    find_package(pmp 3.0.0 CONFIG REQUIRED)
  else()
    find_package(pmp 3.0.0 CONFIG)
  endif()
endif()

if(NOT pmp_FOUND)
  FetchContent_declare(
    pmp
    URL https://github.com/pmp-library/pmp-library/archive/refs/tags/3.0.0.zip
    URL_HASH MD5=7b7f9ce07a7a687c9d78a6583cf64a2c
    SOURCE_SUBDIR Non-Existing # Use non-existing source dir to disable add_subdirectory call of MakeAvailable
    UPDATE_DISCONNECTED ON
    # FIND_PACKAGE_ARGS NAMES CGAL GLOBAL
    OVERRIDE_FIND_PACKAGE
  )

  list(APPEND MAKE_AVAIL_LIST pmp)
  find_package_override_helper(pmp 3.0.0 AnyNewerVersion)
else()
  print_package_info(pmp)
endif()

if(MAKE_AVAIL_LIST)
  message(STATUS "Fetching the following third-party dependencies:  ${MAKE_AVAIL_LIST}")
  foreach(tpl IN LISTS MAKE_AVAIL_LIST)
    message(STATUS "-------------------------------------------------------------------")
    message(STATUS "- Configuring ${tpl}...")
    message(STATUS "-------------------------------------------------------------------")
    FetchContent_MakeAvailable(${tpl})

    message(STATUS "-------------------------------------------------------------------")
    message(STATUS "- Configuring ${tpl}...done")
    message(STATUS "-------------------------------------------------------------------")
  endforeach()
endif()
