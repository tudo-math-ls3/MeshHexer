# HexMesherV2

A geometry grab-bag library.

## Dependencies
To build HexMesherV2 you need:
- A C++17 compatible C++ compiler
- CMake version 3.28.3 or newer
- CGAL version 6.0.1 or newer
- Boost version 1.81 or newer
 
CGAL and Boost will be downloaded and built automatically if they are not found by CMake.

## Building

Assuming you cloned this repo to ``~/hexmesherv2`` you can build the library via:
```
mkdir build
cd build
cmake ~/hexmesherv2
cmake --build . --target HexMesher
```

Additionally run
```
cmake --build . --target hexmesher-cli
```
to build the corresponding command-line application. Or just run ``cmake --build .``.

## Installing

To install the library (and its dependencies, if they were built during the configure step) run
```
cmake --install . --prefix=path/to/install/dir
```

HexMesher will be installed as a CMake package and can be found via ``find_package(HexMesherV2)`` if the install directory is seen by CMake.
After finding HexMesher a HexMesher::HexMesher target is available to link against.
See src/hexmesher.hpp for the public API.

## hexmesher-cli

``hexmesher-cli`` is a command-line application that gives access to HexMeshers's features without having to integrate the library into your own application.

Call as
```
hexmesher-cli [-h|--help] <mesh> <command> [<args>]
```
Running ``hexmesher-cli --help`` gives more detailed information.

The possible commands are

- ``min-gap``: Heuristically try to find the smallest "real" gap between opposite faces of a surface mesh. Tries to discard gaps caused by local geometry, such as in that sharp edges of the mesh, and gaps caused by badly meshes surfaces.

- ``report``: Generate a report about the given mesh. Includes topological information like whether the mesh is closed and the default orientation of surface normals, as well as defects like degenerate triangles.

- ``warnings``: Generates a list of warnings about bad triangles of the mesh. Currently checks for self-intersecting triangles, degenerate triangles, and highly anisotropic triangles. Pass ``--summarize`` to get a summary of warnings, instead of printing each warning individually.


