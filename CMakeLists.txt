cmake_minimum_required(VERSION 3.28.0)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "HexMesherV2 build type")
endif()

################
# Setup project
################
option(HEXMESHER_BUILD_APPLICATIONS "should targets for applications be added?" ON)
option(HEXMESHER_INSTALL "should installation of HexMesher be set up?" ON)
option(HEXMESHER_TESTING "should testing be enabled?" ON)
option(HEXMESHER_PREFER_EXTERNAL_TPL "try to find external libraries before fetching" ON)
option(HEXMESHER_ALLOW_EXTERNAL_DOWNLOAD "enable download of TPLs" ON)
option(HEXMESHER_TPL_CACHE_DIRECTORY "directory for caching thirdparty dependencies" "")

project(
  HexMesherV2
  VERSION 0.1
  DESCRIPTION "Meshing tool for hexahedral meshes"
  LANGUAGES CXX
)

if(HEXMESHER_TESTING)
  enable_testing()
endif()

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)


# Setup thirdparty dependencies.
# Ensures targets Boost::boost and CGAL::CGAL are available to link against
include(Thirdparty)

###################
# Compiler options
###################

# Add library for common compiler options
add_library(hexmesher-settings INTERFACE EXCLUDE_FROM_ALL)

set_target_properties(hexmesher-settings PROPERTIES
  CXX_STANDARD 17
)

target_compile_options(hexmesher-settings INTERFACE
  -Wall
  -Wextra
  -Wundef
  -Wshadow
  -Woverloaded-virtual
  -Wuninitialized
  -Wvla
  -Wdouble-promotion
  -Wformat=2
  -Wnonnull
  -Wswitch-bool
  -Wsizeof-array-argument
  -Wbool-compare
  -Wsuggest-override
  -Wnon-virtual-dtor
  -Wdelete-non-virtual-dtor
  -Wshift-negative-value
  -Wduplicated-cond
  -Wduplicated-branches
  -Wrestrict
  -Wdangling-else
  -Wnonnull
  -Wrestrict
  -Walloc-zero
  -Wparentheses
  -Wpedantic
)

#########
# OpenMP
#########

option(HEXMESHER_HAVE_OMP "Build HexMesherV2 with OpenMP" ON)
if(HEXMESHER_HAVE_OMP)
  find_package(OpenMP REQUIRED)
  target_link_libraries(hexmesher-settings INTERFACE OpenMP::OpenMP_CXX)
endif()

##########
# Testing
##########

if(HEXMESHER_TESTING)
  # Add Catch2 scripts to module path if integrated via FetchContent
  if(NOT Catch2_FOUND)
    message(STATUS "Adding ${catch2_SOURCE_DIR}/extras to module path")
    list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
  endif()

  add_custom_target(tests)
endif()

#####################
# Add subdirectories
#####################
add_subdirectory(src)

if(HEXMESHER_BUILD_APPLICATIONS)
  add_subdirectory(applications)
endif()

##############
# Config file
##############

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/hexmesher_config.hpp.in ${CMAKE_CURRENT_BINARY_DIR}/hexmesher_config.hpp)
target_sources(HexMesher
  PRIVATE FILE_SET private_headers
  TYPE HEADERS
  BASE_DIRS ${CMAKE_CURRENT_BINARY_DIR}
  FILES ${CMAKE_CURRENT_BINARY_DIR}/hexmesher_config.hpp)

############
# Packaging
############

# Setup alias for inclusion via FetchContent
add_library(HexMesher::HexMesher ALIAS HexMesher)

# Setup install
if(HEXMESHER_INSTALL)
  include(GNUInstallDirs)
  include(CMakePackageConfigHelpers)

  # Create package version file
  write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/hexmesherv2-config-version.cmake
    COMPATIBILITY AnyNewerVersion
  )

  # Install package version file and package config file
  install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/hexmesherv2-config-version.cmake
    ${CMAKE_CURRENT_SOURCE_DIR}/hexmesherv2-config.cmake
    DESTINATION ${CMAKE_INSTALL_DATADIR}/hexmesherv2
  )

  # Add target to export
  install(
    TARGETS HexMesher hexmesher-settings
    EXPORT hexmesherv2Targets
    FILE_SET public_headers
  )

  # Generate hexmesherTargets.cmake
  install(
    EXPORT hexmesherv2Targets
    NAMESPACE HexMesher::
    DESTINATION ${CMAKE_INSTALL_DATADIR}/hexmesherv2
  )
endif()
