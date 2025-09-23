#include <catch2/catch_test_macros.hpp>

#include <meshhexer.hpp>
#include <io.hpp>
#include <types.hpp>

#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

using namespace MeshHexer;

TEST_CASE("Write FEAT mesh format", "[IO]")
{
  VolumeMesh vmesh(2, 2, 2);

  vmesh.xs_begin()[0] = Slice(0.0, 0);
  vmesh.xs_begin()[1] = Slice(1.0, 0);

  vmesh.ys_begin()[0] = Slice(0.0, 0);
  vmesh.ys_begin()[1] = Slice(1.0, 0);

  vmesh.zs_begin()[0] = Slice(0.0, 0);
  vmesh.zs_begin()[1] = Slice(1.0, 0);

  std::ostringstream output;

  vmesh.write_feat_xml(output);

  constexpr std::string_view expected(
    "<FeatMeshFile version=\"1\" mesh=\"conformal:hypercube:3:3\">\n"
    "<Mesh type=\"conformal:hypercube:3:3\" size=\"8 12 6 1\">\n"
    "<Vertices>\n"
    "0 0 0\n"
    "1 0 0\n"
    "0 1 0\n"
    "1 1 0\n"
    "0 0 1\n"
    "1 0 1\n"
    "0 1 1\n"
    "1 1 1\n"
    "</Vertices>\n"
    "<Topology dim=\"1\">\n"
    "0 1\n"
    "2 3\n"
    "4 5\n"
    "6 7\n"
    "0 2\n"
    "1 3\n"
    "4 6\n"
    "5 7\n"
    "0 4\n"
    "1 5\n"
    "2 6\n"
    "3 7\n"
    "</Topology>\n"
    "<Topology dim=\"2\">\n"
    "0 1 2 3\n"
    "4 5 6 7\n"
    "0 1 4 5\n"
    "2 3 6 7\n"
    "0 2 4 6\n"
    "1 3 5 7\n"
    "</Topology>\n"
    "<Topology dim=\"3\">\n"
    "0 1 2 3 4 5 6 7\n"
    "</Topology>\n"
    "</Mesh>\n"
    "</FeatMeshFile>\n"
  );

  REQUIRE(output.str() == expected);
}
