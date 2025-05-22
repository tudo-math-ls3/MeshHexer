#include <iostream>

#include <hexmesher.hpp>
#include <io.hpp>

#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/measure.h>

int main(int argc, char* argv[])
{
  std::cout << "Running hexmesher\n";

  if(argc > 1)
  {
    const std::string mode(argv[1]);
    const std::string filename(argv[2]);

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
      std::cerr << "Input mesh is not a triangle mesh\n";
      return 1;
    }

    std::vector<HexMesher::PolygonWithHoles> union_components;

    if(mode == "radial")
    {
      HexMesher::Point p(std::atof(argv[3]), std::atof(argv[4]), std::atof(argv[5]));
      HexMesher::Vector up(std::atof(argv[6]), std::atof(argv[7]), std::atof(argv[8]));
      HexMesher::Vector u(std::atof(argv[9]), std::atof(argv[10]), std::atof(argv[11]));

      int steps = std::atoi(argv[12]);

      HexMesher::RadialCrossSectionSampler sampler(steps, p, u, up);

      union_components = HexMesher::union_of_cross_sections(mesh, sampler);
    }

    if(mode == "line")
    {
      HexMesher::Point start(std::atof(argv[3]), std::atof(argv[4]), std::atof(argv[5]));
      HexMesher::Point end(std::atof(argv[6]), std::atof(argv[7]), std::atof(argv[8]));
      HexMesher::Vector normal(std::atof(argv[9]), std::atof(argv[10]), std::atof(argv[11]));
      HexMesher::Vector up(std::atof(argv[12]), std::atof(argv[13]), std::atof(argv[14]));

      int steps = std::atoi(argv[15]);

      HexMesher::LineCrossSectionSampler sampler(steps, start, end, normal, up);

      union_components = HexMesher::union_of_cross_sections(mesh, sampler);
    }

    if(mode == "thickness")
    {
      std::cout << "Computing thickness\n";
      HexMesher::compute_mesh_thickness(mesh);
      std::cout << "Computing topological distances\n";
      HexMesher::topological_distances(mesh, "f:MIS_id");
      std::ofstream output("thickness.ply");

      if(output)
      {
        CGAL::IO::write_PLY(output, mesh);
      }
      return 0;
    }

    HexMesher::Real h = 100;
    for(auto edge_iter = mesh.edges_begin(); edge_iter != mesh.edges_end(); edge_iter++)
    {
      h = std::min(CGAL::Polygon_mesh_processing::edge_length(*edge_iter, mesh), h);
    }

    for(int i(0); i < union_components.size(); i++)
    {
      HexMesher::PolygonWithHoles& component = union_components[i];

      int idx;
      int total_vertices = component.outer_boundary().size();
      for(const HexMesher::Polygon& hole : component.holes())
      {
        std::cout << "  Hole " << idx << " with " << hole.size() << " vertices\n";
        total_vertices += hole.size();
      }
      std::cout << "Component " << i << ":\n";
      std::cout << "  Boundary Vertices: " << component.outer_boundary().size() << "\n";
      std::cout << "  Number of holes: " << component.holes().size() << "\n";
      std::cout << "  Total vertices: " << total_vertices << "\n";

      // Output the found shadow

      auto pred = [](const std::vector<HexMesher::Vector2D>& normals, const HexMesher::Vector2D& next)
      {
        return std::abs(HexMesher::angle(normals.back(), next)) < 10.0 && std::abs(HexMesher::angle(normals.front(), next)) < 45.0;
      };

      HexMesher::Polygon simplified_boundary = HexMesher::simplify_by_normal(component.outer_boundary(), pred).first;
      HexMesher::Polygon grid_sampled_boundary = HexMesher::grid_sample(component.outer_boundary(), h);

      HexMesher::write_polygon("shadow_" + std::to_string(i) + ".vtp", component);
      HexMesher::write_polygon("simplified_" + std::to_string(i) + ".vtp", simplified_boundary);
      HexMesher::write_geo("union_" + std::to_string(i) + ".geo", simplified_boundary);
      //HexMesher::write_geo_compound_2d("union_2d_" + std::to_string(i) + ".geo", simplified_boundary);

      HexMesher::write_polygon("grid_sampled_" + std::to_string(i) + ".vtp", grid_sampled_boundary);
      HexMesher::write_geo("grid_sampled_" + std::to_string(i) + ".geo", grid_sampled_boundary);
    }
  }
  else
  {
    const std::string filename("/home/user/mmuegge/nobackup/repos/feat/data/models/scalexa_gendie_simple.off");

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
      std::cerr << "Input mesh is not a triangle mesh\n";
      return 1;
    }

    HexMesher::compute_mesh_thickness(mesh);
    std::ofstream output("thickness.ply");

    if(output)
    {
      CGAL::IO::write_PLY(output, mesh);
    }

    /*
    HexMesher::Point start(0.0, 0.0, 80.0);
    HexMesher::Point end(0.0, 0.0, 60.0);
    HexMesher::Vector normal(0.0, 0.0, 1.0);
    HexMesher::Vector up(0.0, 1.0, 0.0);

    HexMesher::Point origin(0.0, 0.0, 0.0);
    HexMesher::Vector up(0.0, 0.0, 1.0);
    HexMesher::Vector normal(1.0, 0.0, 0.0);

    HexMesher::LineCrossSectionSampler sampler(50, start, end, normal, up);
    HexMesher::RadialCrossSectionSampler sampler(50, origin, normal, up);

    std::vector<HexMesher::PolygonWithHoles> union_components = HexMesher::union_of_cross_sections(mesh, sampler);


    HexMesher::Real h = 100;
    for(auto edge_iter = mesh.edges_begin(); edge_iter != mesh.edges_end(); edge_iter++)
    {
      h = std::min(CGAL::Polygon_mesh_processing::edge_length(*edge_iter, mesh), h);
    }

    for(int i(0); i < union_components.size(); i++)
    {
      HexMesher::PolygonWithHoles& component = union_components[i];

      int idx;
      int total_vertices = component.outer_boundary().size();
      for(const HexMesher::Polygon& hole : component.holes())
      {
        std::cout << "  Hole " << idx << " with " << hole.size() << " vertices\n";
        total_vertices += hole.size();
      }
      std::cout << "Component " << i << ":\n";
      std::cout << "  Boundary Vertices: " << component.outer_boundary().size() << "\n";
      std::cout << "  Number of holes: " << component.holes().size() << "\n";
      std::cout << "  Total vertices: " << total_vertices << "\n";

      // Output the found shadow

      auto pred = [](const std::vector<HexMesher::Vector2D>& normals, const HexMesher::Vector2D& next)
      {
        return HexMesher::angle(normals.back(), next) < 15.0 && angle(normals.front(), next) < 45.0;
      };

      //HexMesher::Polygon simplified_boundary = HexMesher::simplify_by_normal(component.outer_boundary(), pred).first;
      HexMesher::Polygon grid_sampled_boundary = HexMesher::grid_sample(component.outer_boundary(), h);

      HexMesher::write_polygon("shadow_" + std::to_string(i) + ".vtp", component);
      //HexMesher::write_polygon("simplified_" + std::to_string(i) + ".vtp", simplified_boundary);
      //HexMesher::write_geo("union_" + std::to_string(i) + ".geo", simplified_boundary);
      //HexMesher::write_geo_compound_2d("union_2d_" + std::to_string(i) + ".geo", simplified_boundary);

      HexMesher::write_polygon("grid_sampled_" + std::to_string(i) + ".vtp", grid_sampled_boundary);
      HexMesher::write_geo("grid_sampled_" + std::to_string(i) + ".geo", grid_sampled_boundary);
    }
    */
  }

  return 0;
}
