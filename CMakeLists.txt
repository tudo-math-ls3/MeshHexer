cmake_minimum_required(VERSION 3.28.0)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "HexMesherV2 build type")
endif()

################
# Setup project
################
project(
  HexMesherV2
  VERSION 0.1
  DESCRIPTION "Meshing tool for hexahedral meshes"
  LANGUAGES CXX
)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

option(HEXMESHER_BUILD_APPLICATIONS "should targets for application be added?" ON)
option(HEXMESHER_INSTALL "should installation of HexMesher be set up?" OFF)
option(HEXMESHER_PREFER_EXTERNAL_TPL "try to find external libraries before fetching" ON)
option(HEXMESHER_ALLOW_EXTERNAL_DOWNLOAD "enable download of TPLs" ON)

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

#####################
# Add subdirectories
#####################
add_subdirectory(src)

if(HEXMESHER_BUILD_APPLICATIONS)
  add_subdirectory(applications)
endif()


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
