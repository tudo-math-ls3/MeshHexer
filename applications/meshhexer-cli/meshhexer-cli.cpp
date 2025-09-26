#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

#include <meshhexer/meshhexer.hpp>
#include <string_view>

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
      std::cout << min_gap.diameter << "\n";
    }

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

  const static char* const usage =
    "meshhexer-cli: Commandline tool for the MeshHexer library\n"
    "\n"
    "Usage:\n"
    "meshhexer-cli [<global args>] <command> <mesh> [<args>]\n"
    "\n"
    "Global options:\n"
    "\t-h, --help\n"
    "\t\tProduce this help text\n"
    "\t--checkpoint-path\n"
    "\t\tWrite checkpoint file with intermediary values. File extension must be.ply\n"
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

  struct GlobalParameters
  {
    /// If true, show help text and end program
    bool show_help = false;

    /// Checkpoint file location. Checkpoint file is written if path is not empty.
    std::filesystem::path checkpoint_path;
  };


  static bool cmp_argument(const char* parameter, char* arg)
  {
    std::size_t n = std::min(std::strlen(parameter), std::strlen(arg));
    return std::strncmp(parameter, arg, n) == 0;
  }

  static char* parse_argument(const char* parameter, int& argc, char*** argv)
  {
    char* arg = (*argv)[0];

    while((*arg != 0) && *arg != '=')
    {
      arg++;
    }

    if(*arg == 0)
    {
      // --param arg case

      if(argc < 2)
      {
        std::cerr << "Error parsing parameter " << parameter << ". Expected additional argument\n";
        std::cerr << usage;
        std::exit(1);
      }

      argc--;
      (*argv)++;

      char* result = (*argv)[0];

      argc--;
      (*argv)++;

      return result;
    }
    else
    {
      // --param=arg case

      // We used up one argument
      argc--;
      (*argv)++;

      return (arg + 1);
    }
  }

  static void handle_flag(bool& b, int& argc, char*** argv)
  {
    b = true;
    argc--;
    (*argv)++;
  }

  /**
   * \brief Parse global parameters
   */
  static GlobalParameters parse_args(int& argc, char*** argv)
  {
    GlobalParameters result;

    // Skip first parameter
    argc--;
    *argv = &(*argv)[1];

    while(argc > 0)
    {
      char* cmd = (*argv)[0];

      if(cmp_argument("--help", cmd) || cmp_argument("-h", cmd))
      {
        handle_flag(result.show_help, argc, argv);
      }
      else if(cmp_argument("--checkpoint-path", cmd))
      {
        char* path = parse_argument("--checkpoint-path", argc, argv);
        result.checkpoint_path = path;
      }
      else
      {
        // Unknown argument
        break;
      }
    }

    return result;
  }

  int main(int argc, char* argv[])
  {
    int orig_argc = argc;
    char** orig_argv = argv;

    GlobalParameters params = parse_args(argc, &argv);

    if(params.show_help)
    {
      std::cout << usage;
      return 0;
    }

    // Handle positional arguments
    std::string command;
    if(argc > 0)
    {
      command = argv[0];
      if(command != "fbm-mesh" && command != "min-gap" && command != "report" && command != "warnings")
      {
        std::cerr << "Unknown command " << command << "\n\n";
        std::cerr << usage;
        exit(1);
      }
      argc--;
      argv++;
    }
    else
    {
      std::cerr << "Expected <command> argument\n\n";
      std::cerr << usage;
      exit(1);
    }

    std::string filename;
    if(argc > 0)
    {
      filename = argv[0];
      argc--;
      argv++;
    }
    else
    {
      std::cerr << "Expected <mesh> argument\n\n";
      std::cerr << usage;
      exit(1);
    }

    MeshHexer::Result<MeshHexer::SurfaceMesh, std::string> result = MeshHexer::load_from_file(filename, true);
    if(result.is_err())
    {
      std::cout << "Reading mesh failed with error: " << result.err_ref() << "\n";
      exit(1);
    }
    MeshHexer::SurfaceMesh mesh = std::move(result).take_ok();

    if(command == "fbm-mesh")
    {
      std::uint64_t levels = 0;
      if(argc > 1)
      {
        if(strcmp(argv[0], "--levels") == 0)
        {
          levels = std::strtoull(argv[1], nullptr, 10);
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

    if(command == "min-gap")
    {
      MeshHexer::Gap min_gap = mesh.min_gap();
      print_min_gap(min_gap);

      if(!params.checkpoint_path.empty())
      {
        MeshHexer::Result<void, std::string> result = mesh.write_to_file(params.checkpoint_path);

        if(result.is_err())
        {
          std::cout << "Writing checkpoint failed with error: " << result.err_ref() << "\n";
        }
      }
    }

    if(command == "report")
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

    if(command == "warnings")
    {
      MeshHexer::MeshWarnings warnings = mesh.warnings();

      bool summarize = false;
      if(argc > 0)
      {
        if(strcmp(argv[0], "--summarize") == 0)
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
