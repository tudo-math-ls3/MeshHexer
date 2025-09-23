#include <algorithm>
#include <array>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <string>

#include <meshhexer.hpp>

namespace MeshHexerCLI::Markdown
{
  static std::string h1(const std::string& heading)
  {
    return heading + "\n" + std::string(heading.size(), '=');
  }

  static std::string h2(const std::string& heading)
  {
    return heading + "\n" + std::string(heading.size(), '-');
  }

  static std::string li(const std::string& content)
  {
    return "- " + content;
  }
} // namespace MeshHexerCLI::Markdown

namespace MeshHexerCLI
{
  namespace
  {
    void print_min_gap(const MeshHexer::Gap& min_gap)
    {
      std::cout << "Min-gap of " << min_gap.diameter << " between faces " << min_gap.face << " and "
                << min_gap.opposite_face << "\n";
      std::cout << "Use `SelectIDs(IDs=[0, " << min_gap.face << ", 0, " << min_gap.opposite_face
                << "], FieldType='CELL')` to select the chosen triangles in ParaView\n";
    }

    // NOTE(mmuegge): This really belongs in the io header. Maybe introduce a public io header?
    template<typename Iter>
    void write_range_as_mtx(std::ostream& stream, Iter begin, Iter end)
    {
      // NOTE(mmuegge): For compatability with FEAT3 we write everything as
      // real. We could inspect the type produced by the iterator and set
      // integer as type if appropriate, but then we would also need to update
      // the parsing logic in FEAT3.
      stream << "%%MatrixMarket matrix array real general\n";

      const auto size = std::distance(begin, end);
      stream << size << " 1\n";

      for(Iter it = begin; it != end; it++)
      {
        stream << *it << "\n";
      }
    }
  } // namespace

  const char* usage = "meshhexer-cli: Commandline tool for the MeshHexer library\n"
                      "\n"
                      "Usage:\n"
                      "meshhexer-cli [-h|--help] <mesh> <command> [<args>]\n"
                      "\n"
                      "Options:\n"
                      "\t-h, --help\n"
                      "\t\t Produce this help text\n"
                      "\n"
                      "Commands:\n"
                      "\tfbm-mesh\n"
                      "\t\tGenerate a non-fitting volume mesh from the surface mesh.\n"
                      "\t\tThe mesh is output as fbm_mesh.xml in FEAT3's mesh file format.\n"
                      "\t\tA separate fbm_mesh.mtx file is created with recommended adaptive\n"
                      "\t\trefinement levels for all vertices.\n"
                      "\t\tThe output mesh is constructed such that all cells are about the same\n"
                      "\t\tsize as their local min-gaps.\n"
                      "\n"
                      "\t\tOptions:\n"
                      "\t\t--levels: Set size of multigrid-hierarchy. If passed, the mesh will be constructed\n"
                      "\t\tsuch that the finest level of the hierarchy matches the min-gaps.\n"
                      "\n"
                      "\tmin-gap\n"
                      "\t\tCalculate smallest inside gap between opposite faces of the mesh\n"
                      "\n"
                      "\treport\n"
                      "\t\tPrint information about the mesh\n"
                      "\n"
                      "\twarnings\n"
                      "\t\tPrint warnings about the mesh. Warns about self-intersections,\n"
                      "\t\tdegenerate triangles, and anisotropic triangles.\n"
                      "\n"
                      "\t\tOptions:\n"
                      "\n"
                      "\t\t--summarize: Summarize warnings\n";

  int main(int argc, char* argv[])
  {
    if(argc < 2)
    {
      std::cout << "Invalid number of arguments\n";
      std::cout << usage;
      return 1;
    }

    if(strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0)
    {
      std::cout << usage;
      return 0;
    }

    if(argc < 3)
    {
      std::cout << "Invalid number of arguments\n";
      std::cout << usage;
      return 1;
    }

    const std::string filename(argv[1]);
    const std::string mode(argv[2]);

    MeshHexer::Result<MeshHexer::SurfaceMesh, std::string> result = MeshHexer::load_from_file(filename, true);

    if(result.is_err())
    {
      std::cout << "Reading mesh failed with error: " << result.err_ref() << "\n";
      exit(1);
    }

    MeshHexer::SurfaceMesh mesh = std::move(result).take_ok();

    if(mode == "fbm-mesh")
    {
      std::uint64_t levels = 0;
      if(argc > 3)
      {
        if(strcmp(argv[3], "--levels") == 0)
        {
          levels = std::strtoull(argv[3], nullptr, 10);
        }
      }
      MeshHexer::VolumeMesh vmesh = mesh.fbm_mesh(levels);

      std::ofstream mesh_file("fbm_mesh.xml");

      if(mesh_file.fail())
      {
        std::cerr << "Error opening fbm_mesh.xml for writing\n";
        return 1;
      }

      vmesh.write_feat_xml(mesh_file);

      std::ofstream mtx_file("fbm_mesh.mtx");

      if(mtx_file.fail())
      {
        std::cerr << "Error opening fbm_mesh.mtx for writing\n";
      }

      std::vector<std::uint64_t> sdls;
      sdls.reserve(vmesh.num_vertices());
      for(std::size_t i(0); i < vmesh.num_vertices(); i++)
      {
        sdls.push_back(vmesh.subdivision_level(i));
      }

      write_range_as_mtx(mtx_file, sdls.begin(), sdls.end());
    }

    if(mode == "min-gap")
    {
      MeshHexer::Gap min_gap = mesh.min_gap();
      print_min_gap(min_gap);

      MeshHexer::Result<void, std::string> result = mesh.write_to_file("thickness.ply");

      if(result.is_err())
      {
        std::cout << "Writing mesh failed with error: " << result.err_ref() << "\n";
      }
    }

    if(mode == "report")
    {
      std::filesystem::path absolute_path = std::filesystem::canonical(filename);

      MeshHexer::MeshWarnings warnings = mesh.warnings();

      std::cout << Markdown::h1("Mesh-Report for " + absolute_path.filename().string()) << "\n";

      std::cout << "\n";

      std::cout << Markdown::h2("Metadata") << "\n";
      std::cout << Markdown::li("Path: " + absolute_path.string()) << "\n";

      std::cout << "\n";

      std::cout << Markdown::h2("Topology") << "\n";
      std::cout << Markdown::li("Number of vertices: " + std::to_string(mesh.num_vertices())) << "\n";
      std::cout << Markdown::li("Number of edges: " + std::to_string(mesh.num_edges())) << "\n";
      std::cout << Markdown::li("Number of faces: " + std::to_string(mesh.num_faces())) << "\n";

      MeshHexer::BoundingBox bb = mesh.bounding_box();
      // clang-format off
      std::cout << Markdown::li("Extent (x, y, z): [" +
        std::to_string(bb.min.x) + ", " +
        std::to_string(bb.max.x) + "] x [" +
        std::to_string(bb.min.y) + ", " +
        std::to_string(bb.max.y) + "] x [" +
        std::to_string(bb.min.z) + ", " +
        std::to_string(bb.max.z) + "]") << "\n";
      // clang-format on

      std::cout << Markdown::li("Is closed: " + std::string(mesh.is_closed() ? "True" : "False")) << "\n";
      std::cout << Markdown::li(
                     "Is wound consistently: " + std::string(mesh.is_wound_consistently() ? "True" : "False"))
                << "\n";
      std::cout << Markdown::li("Is oriented outward: " + std::string(mesh.is_outward_oriented() ? "True" : "False"))
                << "\n";
      std::cout << Markdown::li("Minimal triangle aspect ratio: " + std::to_string(mesh.minimal_aspect_ratio()))
                << "\n";
      std::cout << Markdown::li("Maximal triangle aspect ratio: " + std::to_string(mesh.maximal_aspect_ratio()))
                << "\n";

      std::cout << "\n";

      std::cout << Markdown::h2("Defects") << "\n";
      std::cout << Markdown::li("Self-intersections: " + std::to_string(warnings.self_intersections.size())) << "\n";
      std::cout << Markdown::li("Degenerate triangles: " + std::to_string(warnings.degenerate_triangles.size()))
                << "\n";
      std::cout << Markdown::li("Anisotropic triangles: " + std::to_string(warnings.anisotropic_triangles.size()))
                << "\n";

      std::cout << "\n";
    }

    if(mode == "warnings")
    {
      MeshHexer::MeshWarnings warnings = mesh.warnings();

      bool summarize = false;
      if(argc > 3)
      {
        if(strcmp(argv[3], "--summarize") == 0)
        {
          summarize = true;
        }
      }

      if(summarize)
      {
        std::cout << std::to_string(warnings.self_intersections.size()) << " x Self-intersection of mesh ["
                  << MeshHexer::SelfIntersectionWarning::name << "]\n";
      }
      else
      {
        for(MeshHexer::SelfIntersectionWarning& warning : warnings.self_intersections)
        {
          std::cout << "Self-intersection of mesh between triangle " << std::to_string(warning.tri_a)
                    << " and triangle " << std::to_string(warning.tri_b) << " ["
                    << MeshHexer::SelfIntersectionWarning::name << "]\n";
        }
      }

      if(summarize)
      {
        std::cout << std::to_string(warnings.degenerate_triangles.size()) << " x Triangle with colinear coordinates ["
                  << MeshHexer::DegenerateTriangleWarning::name << "]\n";
      }
      else
      {
        for(MeshHexer::DegenerateTriangleWarning& warning : warnings.degenerate_triangles)
        {
          std::cout << "Coordinates of triangle " << std::to_string(warning.idx) << " are colinear ["
                    << MeshHexer::DegenerateTriangleWarning::name << "]\n";
        }
      }

      if(summarize)
      {
        std::cout << std::to_string(warnings.anisotropic_triangles.size())
                  << " x Triangle with aspect ratio greater than 15 [" << MeshHexer::DegenerateTriangleWarning::name
                  << "]\n";
      }
      else
      {
        for(MeshHexer::AnisotropicTriangleWarning& warning : warnings.anisotropic_triangles)
        {
          std::cout << "Triangle " << std::to_string(warning.idx) << " has aspect ratio greater than 15 ["
                    << MeshHexer::DegenerateTriangleWarning::name << "]\n";
        }
      }
    }

    return 0;
  }
} // namespace MeshHexerCLI

int main(int argc, char* argv[])
{
  MeshHexerCLI::main(argc, argv);
}
