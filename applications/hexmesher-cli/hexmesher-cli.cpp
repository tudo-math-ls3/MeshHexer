#include <array>
#include <cstring>
#include <filesystem>
#include <iostream>
#include <string>
#include <algorithm>
#include <functional>
#include <fstream>
#include <iostream>

#include <hexmesher.hpp>

namespace HexMesherCLI::Markdown
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
} // namespace HexMesherCLI::Markdown

namespace HexMesherCLI
{
  namespace
  {
    void print_min_gap(const HexMesher::Gap& min_gap)
    {
      std::cout << "Min-gap of " << min_gap.diameter << " between faces " << min_gap.face << " and " << min_gap.opposite_face
                << "\n";
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
  }

  const char* usage = "hexmesher-cli: Commandline tool for the hexmesher library\n"
                      "\n"
                      "Usage:\n"
                      "hexmesher-cli [-h|--help] <mesh> <command> [<args>]\n"
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

    HexMesher::Result<HexMesher::SurfaceMesh, std::string> result = HexMesher::load_from_file(filename, true);

    if(result.is_err())
    {
      std::cout << "Reading mesh failed with error: " << result.err_ref() << "\n";
      exit(1);
    }

    HexMesher::SurfaceMesh mesh = std::move(result).take_ok();

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
      HexMesher::VolumeMesh vmesh = mesh.fbm_mesh(levels);

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
      HexMesher::Gap min_gap = mesh.min_gap();
      print_min_gap(min_gap);

      HexMesher::Result<void, std::string> result = mesh.write_to_file("thickness.ply");

      if(result.is_err())
      {
        std::cout << "Writing mesh failed with error: " << result.err_ref() << "\n";
      }
    }

    if(mode == "report")
    {
      std::filesystem::path absolute_path = std::filesystem::canonical(filename);

      HexMesher::MeshWarnings warnings = mesh.warnings();

      std::cout << Markdown::h1("Mesh-Report for " + absolute_path.filename().string()) << "\n";

      std::cout << "\n";

      std::cout << Markdown::h2("Metadata") << "\n";
      std::cout << Markdown::li("Path: " + absolute_path.string()) << "\n";


      std::cout << "\n";

      std::cout << Markdown::h2("Topology") << "\n";
      std::cout << Markdown::li("Number of vertices: " + std::to_string(mesh.num_vertices())) << "\n";
      std::cout << Markdown::li("Number of edges: " + std::to_string(mesh.num_edges())) << "\n";
      std::cout << Markdown::li("Number of faces: " + std::to_string(mesh.num_faces())) << "\n";

      HexMesher::BoundingBox bb = mesh.bounding_box();
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
      HexMesher::MeshWarnings warnings = mesh.warnings();

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
                  << HexMesher::SelfIntersectionWarning::name << "]\n";
      }
      else
      {
        for(HexMesher::SelfIntersectionWarning& warning : warnings.self_intersections)
        {
          std::cout << "Self-intersection of mesh between triangle " << std::to_string(warning.tri_a)
                    << " and triangle " << std::to_string(warning.tri_b) << " ["
                    << HexMesher::SelfIntersectionWarning::name << "]\n";
        }
      }

      if(summarize)
      {
        std::cout << std::to_string(warnings.degenerate_triangles.size()) << " x Triangle with colinear coordinates ["
                  << HexMesher::DegenerateTriangleWarning::name << "]\n";
      }
      else
      {
        for(HexMesher::DegenerateTriangleWarning& warning : warnings.degenerate_triangles)
        {
          std::cout << "Coordinates of triangle " << std::to_string(warning.idx) << " are colinear ["
                    << HexMesher::DegenerateTriangleWarning::name << "]\n";
        }
      }

      if(summarize)
      {
        std::cout << std::to_string(warnings.anisotropic_triangles.size())
                  << " x Triangle with aspect ratio greater than 15 [" << HexMesher::DegenerateTriangleWarning::name
                  << "]\n";
      }
      else
      {
        for(HexMesher::AnisotropicTriangleWarning& warning : warnings.anisotropic_triangles)
        {
          std::cout << "Triangle " << std::to_string(warning.idx) << " has aspect ratio greater than 15 ["
                    << HexMesher::DegenerateTriangleWarning::name << "]\n";
        }
      }
    }

    return 0;
  }
} // namespace HexMesherCLI

int main(int argc, char* argv[])
{
  HexMesherCLI::main(argc, argv);
}

/*
int main(int argc, char* argv[])
{
  std::cout << "Running hexmesher\n";

  if(argc > 1)
  {
    const std::string mode(argv[1]);
    const std::string filename(argv[2]);

    HexMesher::Mesh mesh;
    std::cout << "Reading mesh " << filename << "\n";
    if(ends_with(filename, ".ply"))
    {
      std::cout << "Reading .ply\n";
      std::ifstream mesh_file(filename);
      std::string comment("");
      if(!CGAL::IO::read_PLY(mesh_file, mesh, comment, true))
      {
        std::cerr << "Could not read mesh!\n";
        return 1;
      }

      std::cout << "Read mesh with double properties:\n";
      for(auto prop : mesh.properties<HexMesher::FaceIndex>())
      {
        std::cout << prop << "\n";
      }
    }
    else if(!CGAL::Polygon_mesh_processing::IO::read_polygon_mesh(filename,
mesh))
    {
      std::cerr << "Could not read mesh!\n";
      return 1;
    }

    if(CGAL::is_empty(mesh))
    {
      std::cerr << "Input mesh is empty!\n";
      return 1;
    }

    if(!CGAL::is_triangle_mesh(mesh))
    {
      CGAL::Polygon_mesh_processing::triangulate_faces(mesh);
    }

    std::vector<HexMesher::PolygonWithHoles2D> union_components;

    if(mode == "radial")
    {
      HexMesher::Point3D p(std::atof(argv[3]), std::atof(argv[4]),
std::atof(argv[5])); HexMesher::Vector3D up(std::atof(argv[6]),
std::atof(argv[7]), std::atof(argv[8])); HexMesher::Vector3D
u(std::atof(argv[9]), std::atof(argv[10]), std::atof(argv[11]));

      int steps = std::atoi(argv[12]);

      HexMesher::RadialCrossSectionSampler sampler(steps, p, u, up);

      union_components = HexMesher::union_of_cross_sections(mesh, sampler);
    }

    if(mode == "line")
    {
      HexMesher::Point3D start(std::atof(argv[3]), std::atof(argv[4]),
std::atof(argv[5])); HexMesher::Point3D end(std::atof(argv[6]),
std::atof(argv[7]), std::atof(argv[8])); HexMesher::Vector3D
normal(std::atof(argv[9]), std::atof(argv[10]), std::atof(argv[11]));
      HexMesher::Vector3D up(std::atof(argv[12]), std::atof(argv[13]),
std::atof(argv[14]));

      int steps = std::atoi(argv[15]);

      HexMesher::LineCrossSectionSampler sampler(steps, start, end, normal, up);

      union_components = HexMesher::union_of_cross_sections(mesh, sampler);
    }

    if(mode == "report")
    {
      std::cout << "Mesh quality report for " << filename << "\n";

      bool is_closed = CGAL::is_closed(mesh);
      std::cout << "Is closed: " << is_closed << "\n";
      if(is_closed)
      {
        std::cout << "Is outward-oriented: " <<
CGAL::Polygon_mesh_processing::is_outward_oriented(mesh) << "\n";
      }

      std::cout << "Is wound consistently: " <<
HexMesher::is_wound_consistently(mesh) << "\n";

      std::vector<std::pair<HexMesher::FaceIndex, HexMesher::FaceIndex>>
self_intersections; CGAL::Polygon_mesh_processing::self_intersections(mesh,
std::back_inserter(self_intersections)); std::cout << "Self-intersections: " <<
self_intersections.size() << "\n"; if(self_intersections.size() > 0)
      {
        std::cout << "List of self-intersections: ";
        for(int i(0); i < self_intersections.size(); i++)
        {
          //std::cout << "(" << self_intersections[i].first << ", " <<
self_intersections[i].second << ")";

          if(i < self_intersections.size() - 1)
          {
            //std::cout << ", ";
          }
        }
        std::cout << "\n";
      }

      std::vector<HexMesher::FaceIndex> degenerate_faces;
      CGAL::Polygon_mesh_processing::degenerate_faces(mesh,
std::back_inserter(degenerate_faces)); std::cout << "Number of degenerate
triangles: " << degenerate_faces.size() << "\n";

      if(degenerate_faces.size() > 0)
      {
        std::cout << "List of degenerate triangles: ";
        for(int i(0); i < degenerate_faces.size(); i++)
        {
          std::cout << degenerate_faces[i];

          if(i < degenerate_faces.size() - 1)
          {
            std::cout << ", ";
          }
        }
        std::cout << "\n";
      }

      // Minimal and maximal aspect ratio
      double min_aspect_ratio = std::numeric_limits<double>::max();
      double max_aspect_ratio = 0;

      for(HexMesher::FaceIndex f : mesh.faces())
      {
        double ratio = CGAL::Polygon_mesh_processing::face_aspect_ratio(f,
mesh); min_aspect_ratio = std::min(min_aspect_ratio, ratio); max_aspect_ratio =
std::max(max_aspect_ratio, ratio);
      }

      std::cout << "Minimum aspect ratio: " << min_aspect_ratio << "\n";
      std::cout << "Maximum aspect ratio: " << max_aspect_ratio << "\n";

      return 0;
    }

    if(mode == "thickness")
    {
      std::cout << "Computing vertex normals...\n";
      HexMesher::compute_vertex_normals(mesh);

      //std::cout << "Computing principal vertex curvatures...\n";
      //HexMesher::compute_curvature(mesh);

      std::cout << "Computing maximum inscribed spheres...\n";
      HexMesher::StopWatch thickness_stopwatch;
      thickness_stopwatch.start();
      HexMesher::compute_mesh_thickness(mesh);
      thickness_stopwatch.stop();

      std::cout << "Finished computing maximum inscribed spheres. Took "
        << thickness_stopwatch.elapsed_string() << " ("
        << (double)mesh.num_faces() / thickness_stopwatch.elapsed() << "
Elements per second)\n";

      std::ofstream output("thickness.ply");

      if(output)
      {
        CGAL::IO::write_PLY(output, mesh);
      }
      return 0;
    }

    if(mode == "min-gap")
    {
      double h = std::numeric_limits<double>::max();
      for(auto edge_iter = mesh.edges_begin(); edge_iter != mesh.edges_end();
edge_iter++)
      {
        h =
std::min(CGAL::to_double(CGAL::Polygon_mesh_processing::edge_length(*edge_iter,
mesh)), h);
      }

      // NOTE(mmuegge): Tried to store these as HexMesher::Points,
      // but CGAL produced a stack overflow on cleanup,
      // because it reuses values of the points in the exact kernel
      // and cleans them up recursively.
      // Try again once we support the inexact kernel.
      double min_x;
      double min_y;
      double min_z;

      double max_x;
      double max_y;
      double max_z;

      for(const HexMesher::Point3D& p : mesh.points())
      {
        min_x = std::min(CGAL::to_double(p.x()), min_x);
        min_y = std::min(CGAL::to_double(p.y()), min_y);
        min_z = std::min(CGAL::to_double(p.z()), min_z);

        max_x = std::max(CGAL::to_double(p.x()), max_x);
        max_y = std::max(CGAL::to_double(p.y()), max_y);
        max_z = std::max(CGAL::to_double(p.z()), max_z);
      }

      double max_mesh_delta = std::max({max_x - min_x, max_y - min_y, max_z -
min_z});

      HexMesher::compute_vertex_normals(mesh);

      if(!mesh.property_map<HexMesher::FaceIndex, double>("f:MIS_diameter") ||
         !mesh.property_map<HexMesher::FaceIndex, std::uint32_t>("f:MIS_id") ||
         !mesh.property_map<HexMesher::FaceIndex,
double>("f:similarity_of_normals"))
        {
        std::cout << "Missing mesh properties f:MIS_diameter and/or f:MIS_id
and/or f:similarity_of_normals. Computing maximum inscribed spheres...";
std::cout.flush();

        HexMesher::StopWatch thickness_stopwatch;
        thickness_stopwatch.start();
        HexMesher::compute_mesh_thickness(mesh);
        thickness_stopwatch.stop();

        std::cout << " done. Took "
          << thickness_stopwatch.elapsed_string() << " ("
          << (double)mesh.num_faces() / thickness_stopwatch.elapsed() << "
Elements per second)\n";
      }

      if(!mesh.property_map<HexMesher::FaceIndex,
double>("f:topological_distance"))
      {
        std::cout << "Missing mesh property f:topological_distance. Computing
topological distances..."; std::cout.flush();

        HexMesher::StopWatch sw;
        sw.start();
        HexMesher::topological_distances(mesh, "f:MIS_id", "f:MIS_diameter");
        sw.stop();

        std::cout << " done. Took "
          << sw.elapsed_string() << " ("
          << (double)mesh.num_faces() / sw.elapsed() << " Elements per
second)\n";
      }

      HexMesher::Mesh::Property_map<HexMesher::FaceIndex, double> topo_dist =
        mesh.property_map<HexMesher::FaceIndex,
double>("f:topological_distance").value();

      HexMesher::Mesh::Property_map<HexMesher::FaceIndex, double> diameter =
        mesh.property_map<HexMesher::FaceIndex,
double>("f:MIS_diameter").value();

      HexMesher::Mesh::Property_map<HexMesher::FaceIndex, double> similarity =
        mesh.property_map<HexMesher::FaceIndex,
double>("f:similarity_of_normals").value();

      // Find maximum topological distance for normalization
      double max_topo_dist = *std::max_element(topo_dist.begin(),
topo_dist.end()); double max_diameter = *std::max_element(diameter.begin(),
diameter.end());

      HexMesher::MinGap min_gap = HexMesher::determine_min_gap_direct(mesh,
[&](HexMesher::FaceIndex idx) { if(topo_dist[idx] == -1.0)
        {
          return diameter[idx];
        }

        if(topo_dist[idx] != 0)
        {
          return diameter[idx] * (1.0 + max_mesh_delta /
std::exp(topo_dist[idx]));
        }
        else
        {
          return max_diameter;
        }
      }, std::string("f:MIS_id"), std::string("f:gap_exp"));
      std::cout << "Min-gap (gap = diameter * (1 + max_topo_dist/e^topo_dist))
is " << min_gap.gap << " between faces "
<< min_gap.origin << " and " << min_gap.limiting << "\n"; std::cout << "Use
`SelectIDs(IDs=[0, " << min_gap.origin.idx()
<< ", 0, " << min_gap.limiting.idx() << "], FieldType='CELL')` to select the
chosen triangles in ParaView\n\n";

      min_gap = HexMesher::determine_min_gap_direct(mesh,
[&](HexMesher::FaceIndex idx) {
        // NOTE (mmuegge): This whole thing is really some kind of fuzzy logic,
        // where some decisions are hard and some are soft via penalties.
        // Im am sure there is formal literature about that.
        if(topo_dist[idx] == -1.0)
        {
          return diameter[idx];
        }

        if(diameter[idx] != 0)
        {
          // Relative size of gap and topological distance. Higher is better.
          double relative_distance = topo_dist[idx] / diameter[idx];
          return diameter[idx] * (1.0 + 10.0 / relative_distance);
        }
        else
        {
          return max_diameter;
        }
      }, std::string("f:MIS_id"), std::string("f:gap_exp"));
      std::cout << "Min-gap (gap = diameter * (1 + 10 / (topo_dist / diameter)))
is " << min_gap.gap << " between faces " << min_gap.origin << " and " <<
min_gap.limiting << "\n"; std::cout << "Use `SelectIDs(IDs=[0, " <<
min_gap.origin.idx() << ", 0, " << min_gap.limiting.idx() << "],
FieldType='CELL')` to select the chosen triangles in ParaView\n\n";

      min_gap = HexMesher::determine_min_gap_direct(mesh,
[&](HexMesher::FaceIndex idx) { if(topo_dist[idx] == 0)
        {
          return max_diameter;
        }

        double similarity_penalty = std::pow(1.0 - similarity[idx], 0.2) *
max_mesh_delta; if(topo_dist[idx] == -1.0)
        {
          return diameter[idx] + similarity_penalty;
        }

        double topo_penalty = max_mesh_delta / std::exp(topo_dist[idx]);
        return diameter[idx] * (1.0 + topo_penalty) + similarity_penalty;
      }, std::string("f:MIS_id"), std::string("f:gap_exp_similarity"));
      std::cout << "Min-gap (gap = diameter * (1 + max_topo_dist/e^topo_dist) +
similarity * mesh_size) is " << min_gap.gap << " between faces " <<
min_gap.origin << " and " << min_gap.limiting << "\n"; std::cout << "Use
`SelectIDs(IDs=[0, " << min_gap.origin.idx() << ", 0, " <<
min_gap.limiting.idx() << "], FieldType='CELL')` to select the chosen triangles
in ParaView\n\n";

      min_gap = HexMesher::determine_min_gap_weighted(mesh,
[&](HexMesher::FaceIndex idx) { double topo_factor(0.0);

        if(topo_dist[idx] > 0)
        {
          topo_factor = topo_dist[idx] / max_topo_dist;
        }

        return (1.0 / 2.0) * topo_factor + (1.0 / 2.0) * similarity[idx];
      }, std::string("f:MIS_diameter"), std::string("f:MIS_id"),
std::string("f:gap_score")); std::cout << "Min-gap (weighted score) is " <<
min_gap.gap << " between faces " << min_gap.origin << " and " <<
min_gap.limiting << "\n"; std::cout << "Use `SelectIDs(IDs=[0, " <<
min_gap.origin.idx() << ", 0, " << min_gap.limiting.idx() << "],
FieldType='CELL')` to select the chosen triangles in ParaView\n\n";

      HexMesher::compute_max_dihedral_angle(mesh);
      HexMesher::score_gaps(mesh);

      std::cout << "\n";

      HexMesher::MinGap first_percentile = HexMesher::select_min_gap(mesh,
0.01); std::cout << "1st percentile min gap by score: ";
      print_min_gap(first_percentile);

      HexMesher::MinGap fifth_percentile = HexMesher::select_min_gap(mesh,
0.05); std::cout << "5th percentile min gap by score: ";
      print_min_gap(fifth_percentile);

      HexMesher::MinGap tenth_percentile = HexMesher::select_min_gap(mesh, 0.1);
      std::cout << "10th percentile min gap by score: ";
      print_min_gap(tenth_percentile);

      HexMesher::MinGap twentieth_percentile = HexMesher::select_min_gap(mesh,
0.20); std::cout << "20th percentile min gap by score: ";
      print_min_gap(twentieth_percentile);

      HexMesher::MinGap fiftieth_percentile = HexMesher::select_min_gap(mesh,
0.50); std::cout << "50th percentile min gap by score: ";
      print_min_gap(fiftieth_percentile);

      HexMesher::MinGap fixed = HexMesher::select_min_gap2(mesh);
      std::cout << "0th percentile min gap by score: ";
      print_min_gap(fixed);

      auto edge_lengths = mesh.property_map<HexMesher::EdgeIndex,
double>("e:lengths"); if(edge_lengths)
      {
        auto map = edge_lengths.value();
        mesh.remove_property_map(map);
      }

      std::ofstream output("thickness.ply");

      if(output)
      {
        CGAL::IO::write_PLY(output, mesh);
      }
      return 0;
    }

    HexMesher::Real h = 100;
    for(auto edge_iter = mesh.edges_begin(); edge_iter != mesh.edges_end();
edge_iter++)
    {
      h = std::min(CGAL::Polygon_mesh_processing::edge_length(*edge_iter, mesh),
h);
    }

    for(int i(0); i < union_components.size(); i++)
    {
      HexMesher::PolygonWithHoles2D& component = union_components[i];

      int idx;
      int total_vertices = component.outer_boundary().size();
      for(const HexMesher::Polygon2D& hole : component.holes())
      {
        std::cout << "  Hole " << idx << " with " << hole.size() << "
vertices\n"; total_vertices += hole.size();
      }
      std::cout << "Component " << i << ":\n";
      std::cout << "  Boundary Vertices: " << component.outer_boundary().size()
<< "\n"; std::cout << "  Number of holes: " << component.holes().size() << "\n";
      std::cout << "  Total vertices: " << total_vertices << "\n";

      // Output the found shadow

      auto pred = [](const std::vector<HexMesher::Vector2D>& normals, const
HexMesher::Vector2D& next)
      {
        return std::abs(HexMesher::angle(normals.back(), next)) < 10.0 &&
std::abs(HexMesher::angle(normals.front(), next)) < 45.0;
      };

      HexMesher::Polygon2D simplified_boundary =
HexMesher::simplify_by_normal(component.outer_boundary(), pred).first;
      HexMesher::Polygon2D grid_sampled_boundary =
HexMesher::grid_sample(component.outer_boundary(), h);

      HexMesher::write_polygon("shadow_" + std::to_string(i) + ".vtp",
component); HexMesher::write_polygon("simplified_" + std::to_string(i) + ".vtp",
simplified_boundary); HexMesher::write_geo("union_" + std::to_string(i) +
".geo", simplified_boundary);
      //HexMesher::write_geo_compound_2d("union_2d_" + std::to_string(i) +
".geo", simplified_boundary);

      HexMesher::write_polygon("grid_sampled_" + std::to_string(i) + ".vtp",
grid_sampled_boundary); HexMesher::write_geo("grid_sampled_" + std::to_string(i)
+ ".geo", grid_sampled_boundary);
    }
  }
  else
  {
    const std::string
filename("/home/user/mmuegge/nobackup/repos/feat/data/models/scalexa_gendie_simple.off");
    //const std::string
filename("/home/user/mmuegge/nobackup/projects/hexmesher/meshes/surface_22630.off");
    //const std::string
filename("/home/user/mmuegge/nobackup/projects/hexmesher/meshes/otto_surface_10088.off");
    //const std::string
filename("/home/user/mmuegge/nobackup/projects/hexmesher/meshes/otto_surface_11598.off");

    HexMesher::Mesh mesh;
    if(!CGAL::Polygon_mesh_processing::IO::read_polygon_mesh(filename, mesh))
    {
      std::cerr << "Could not read mesh!\n";
      return 1;
    }

    if(CGAL::is_empty(mesh))
    {
      std::cerr << "Input mesh is empty!\n";
      return 1;
    }

    if(!CGAL::is_triangle_mesh(mesh))
    {
      CGAL::Polygon_mesh_processing::triangulate_faces(mesh);
    }

    HexMesher::compute_max_dihedral_angle(mesh);
    HexMesher::compute_vertex_normals(mesh);
    HexMesher::compute_mesh_thickness(mesh);
    HexMesher::topological_distances(mesh, "f:MIS_id", "f:MIS_diameter");

    HexMesher::score_gaps(mesh);

    HexMesher::MinGap fifth_percentile = HexMesher::select_min_gap(mesh, 0.05);
    std::cout << "5th percentile min gap by score:\n";
    print_min_gap(fifth_percentile);

    HexMesher::MinGap twentieth_percentile = HexMesher::select_min_gap(mesh,
0.20); std::cout << "20th percentile min gap by score:\n";
    print_min_gap(twentieth_percentile);

    HexMesher::MinGap fiftieth_percentile = HexMesher::select_min_gap(mesh,
0.50); std::cout << "50th percentile min gap by score:\n";
    print_min_gap(fiftieth_percentile);
    std::ofstream output("thickness.ply");

    if(output)
    {
      CGAL::IO::write_PLY(output, mesh);
    }

    HexMesher::Point start(0.0, 0.0, 80.0);
    HexMesher::Point end(0.0, 0.0, 60.0);
    HexMesher::Vector normal(0.0, 0.0, 1.0);
    HexMesher::Vector up(0.0, 1.0, 0.0);

    HexMesher::Point origin(0.0, 0.0, 0.0);
    HexMesher::Vector up(0.0, 0.0, 1.0);
    HexMesher::Vector normal(1.0, 0.0, 0.0);

    HexMesher::LineCrossSectionSampler sampler(50, start, end, normal, up);
    HexMesher::RadialCrossSectionSampler sampler(50, origin, normal, up);

    std::vector<HexMesher::PolygonWithHoles> union_components =
HexMesher::union_of_cross_sections(mesh, sampler);


    HexMesher::Real h = 100;
    for(auto edge_iter = mesh.edges_begin(); edge_iter != mesh.edges_end();
edge_iter++)
    {
      h = std::min(CGAL::Polygon_mesh_processing::edge_length(*edge_iter, mesh),
h);
    }

    for(int i(0); i < union_components.size(); i++)
    {
      HexMesher::PolygonWithHoles& component = union_components[i];

      int idx;
      int total_vertices = component.outer_boundary().size();
      for(const HexMesher::Polygon& hole : component.holes())
      {
        std::cout << "  Hole " << idx << " with " << hole.size() << "
vertices\n"; total_vertices += hole.size();
      }
      std::cout << "Component " << i << ":\n";
      std::cout << "  Boundary Vertices: " << component.outer_boundary().size()
<< "\n"; std::cout << "  Number of holes: " << component.holes().size() << "\n";
      std::cout << "  Total vertices: " << total_vertices << "\n";

      // Output the found shadow

      auto pred = [](const std::vector<HexMesher::Vector2D>& normals, const
HexMesher::Vector2D& next)
      {
        return HexMesher::angle(normals.back(), next) < 15.0 &&
angle(normals.front(), next) < 45.0;
      };

      //HexMesher::Polygon simplified_boundary =
HexMesher::simplify_by_normal(component.outer_boundary(), pred).first;
      HexMesher::Polygon grid_sampled_boundary =
HexMesher::grid_sample(component.outer_boundary(), h);

      HexMesher::write_polygon("shadow_" + std::to_string(i) + ".vtp",
component);
      //HexMesher::write_polygon("simplified_" + std::to_string(i) + ".vtp",
simplified_boundary);
      //HexMesher::write_geo("union_" + std::to_string(i) + ".geo",
simplified_boundary);
      //HexMesher::write_geo_compound_2d("union_2d_" + std::to_string(i) +
".geo", simplified_boundary);

      HexMesher::write_polygon("grid_sampled_" + std::to_string(i) + ".vtp",
grid_sampled_boundary); HexMesher::write_geo("grid_sampled_" + std::to_string(i)
+ ".geo", grid_sampled_boundary);
    }
  }

  return 0;
}
*/
