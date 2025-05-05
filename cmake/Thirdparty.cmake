# Setup thirdparty dependencies
include(FetchContent)
include(FindPackageOverrideHelper)
include(PrintPackageConfig)

set(MAKE_AVAIL_LIST)
mark_as_advanced(MAKE_AVAIL_LIST)

if(HEXMESHER_PREFER_EXTERNAL_TPL OR Boost_DIR)
  if(NOT HEXMESHER_ALLOW_EXTERNAL_DOWNLOAD)
    find_package(Boost 1.88 REQUIRED)
  else()
    find_package(Boost 1.88)
  endif()
endif()

if(NOT Boost_FOUND)
  set(BOOST_INCLUDE_LIBRARIES graph heap logic)

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
  list(APPEND MAKE_AVAIL_LIST CGAL)
  find_package_override_helper(cgal 6.0.1 AnyNewerVersion)
else()
  print_package_info(CGAL)
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