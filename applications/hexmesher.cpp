#include <iostream>

#include <hexmesher.hpp>
#include <io.hpp>
#include <util.hpp>

#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Surface_mesh/IO.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

static bool ends_with(const std::string& string, const std::string& ending)
{
  if(ending.size() > string.size())
  {
    return false;
  }
  return std::equal(ending.rbegin(), ending.rend(), string.rbegin());
}

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
    else if(!CGAL::Polygon_mesh_processing::IO::read_polygon_mesh(filename, mesh))
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
      HexMesher::Point3D p(std::atof(argv[3]), std::atof(argv[4]), std::atof(argv[5]));
      HexMesher::Vector3D up(std::atof(argv[6]), std::atof(argv[7]), std::atof(argv[8]));
      HexMesher::Vector3D u(std::atof(argv[9]), std::atof(argv[10]), std::atof(argv[11]));

      int steps = std::atoi(argv[12]);

      HexMesher::RadialCrossSectionSampler sampler(steps, p, u, up);

      union_components = HexMesher::union_of_cross_sections(mesh, sampler);
    }

    if(mode == "line")
    {
      HexMesher::Point3D start(std::atof(argv[3]), std::atof(argv[4]), std::atof(argv[5]));
      HexMesher::Point3D end(std::atof(argv[6]), std::atof(argv[7]), std::atof(argv[8]));
      HexMesher::Vector3D normal(std::atof(argv[9]), std::atof(argv[10]), std::atof(argv[11]));
      HexMesher::Vector3D up(std::atof(argv[12]), std::atof(argv[13]), std::atof(argv[14]));

      int steps = std::atoi(argv[15]);

      HexMesher::LineCrossSectionSampler sampler(steps, start, end, normal, up);

      union_components = HexMesher::union_of_cross_sections(mesh, sampler);
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
        << (double)mesh.num_faces() / thickness_stopwatch.elapsed() << " Elements per second)\n";

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
      for(auto edge_iter = mesh.edges_begin(); edge_iter != mesh.edges_end(); edge_iter++)
      {
        h = std::min(CGAL::to_double(CGAL::Polygon_mesh_processing::edge_length(*edge_iter, mesh)), h);
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

      double max_mesh_delta = std::max({max_x - min_x, max_y - min_y, max_z - min_z});

      HexMesher::compute_vertex_normals(mesh);

      if(!mesh.property_map<HexMesher::FaceIndex, double>("f:MIS_diameter") ||
         !mesh.property_map<HexMesher::FaceIndex, std::uint32_t>("f:MIS_id") ||
         !mesh.property_map<HexMesher::FaceIndex, double>("f:similarity_of_normals"))
      {
        std::cout << "Missing mesh properties f:MIS_diameter and/or f:MIS_id and/or f:similarity_of_normals. Computing maximum inscribed spheres...\n";

        HexMesher::StopWatch thickness_stopwatch;
        thickness_stopwatch.start();
        HexMesher::compute_mesh_thickness(mesh);
        thickness_stopwatch.stop();

        std::cout << "Finished computing maximum inscribed spheres. Took "
          << thickness_stopwatch.elapsed_string() << " ("
          << (double)mesh.num_faces() / thickness_stopwatch.elapsed() << " Elements per second)\n";
      }

      if(!mesh.property_map<HexMesher::FaceIndex, double>("f:topological_distance"))
      {
        std::cout << "Missing mesh property f:topological_distance. Computing topological distances...\n";

        HexMesher::StopWatch sw;
        sw.start();
        HexMesher::topological_distances(mesh, "f:MIS_id");
        sw.stop();

        std::cout << "Finished computing topological distances. Took "
          << sw.elapsed_string() << " ("
          << (double)mesh.num_faces() / sw.elapsed() << " Elements per second)\n";
      }

      HexMesher::Mesh::Property_map<HexMesher::FaceIndex, double> topo_dist =
        mesh.property_map<HexMesher::FaceIndex, double>("f:topological_distance").value();

      HexMesher::Mesh::Property_map<HexMesher::FaceIndex, double> diameter =
        mesh.property_map<HexMesher::FaceIndex, double>("f:MIS_diameter").value();

      HexMesher::Mesh::Property_map<HexMesher::FaceIndex, double> similarity =
        mesh.property_map<HexMesher::FaceIndex, double>("f:similarity_of_normals").value();

      // Find maximum topological distance for normalization
      double max_topo_dist = *std::max_element(topo_dist.begin(), topo_dist.end());
      double max_diameter = *std::max_element(diameter.begin(), diameter.end());

      std::cout << "Computing possible min-gaps (width * 1/topo_dist)\n";
      double min_gap = HexMesher::determine_min_gap_weighted(mesh, [&](HexMesher::FaceIndex idx) {
        if(topo_dist[idx] != 0)
        {
          return diameter[idx] / topo_dist[idx];
        }
        else
        {
          return max_diameter;
        }},
        std::string("f:MIS_diameter"), std::string("f:gap_ratio"));
      std::cout << "Min-gap (width * 1 / topo_dist) is " << min_gap << "\n";

      std::cout << "Computing possible min-gaps (cutoff topo_dist > 5h)\n";
      min_gap = HexMesher::determine_min_gap_weighted(mesh, [&](HexMesher::FaceIndex idx) {
        if(topo_dist[idx] > 5.0 * h)
        {
          return diameter[idx];
        }
        else
        {
          return max_diameter;
        }
      }, std::string("f:MIS_diameter"), std::string("f:gap_cutoff"));
      std::cout << "Min-gap (cutoff) is " << min_gap << "\n";

      std::cout << "Computing possible min-gaps (gap = diameter * (1 + 1/topo_dist^3))\n";
      min_gap = HexMesher::determine_min_gap_direct(mesh, [&](HexMesher::FaceIndex idx) {
        if(topo_dist[idx] != 0)
        {
          return diameter[idx] * (1.0 + 1.0 / std::pow(topo_dist[idx], 3.0));
        }
        else
        {
          return max_diameter;
        }
      }, std::string("f:gap_cubed"));
      std::cout << "Min-gap (gap = diameter * (1 + 1/topo_dist^3)) is " << min_gap << "\n";

      std::cout << "Computing possible min-gaps (gap = diameter * (1 + max_topo_dist/topo_dist^5))\n";
      min_gap = HexMesher::determine_min_gap_direct(mesh, [&](HexMesher::FaceIndex idx) {
        if(topo_dist[idx] != 0)
        {
          return diameter[idx] * (1.0 + max_topo_dist / std::pow(topo_dist[idx], 5.0));
        }
        else
        {
          return max_diameter;
        }
      }, std::string("f:gap_^5"));
      std::cout << "Min-gap (gap = diameter * (1 + max_topo_dist/topo_dist^5)) is " << min_gap << "\n";

      std::cout << "Computing possible min-gaps (gap = diameter * (1 + max_topo_dist/e^topo_dist))\n";
      min_gap = HexMesher::determine_min_gap_direct(mesh, [&](HexMesher::FaceIndex idx) {
        if(topo_dist[idx] != 0)
        {
          return diameter[idx] * (1.0 + max_topo_dist / std::exp(topo_dist[idx]));
        }
        else
        {
          return max_diameter;
        }
      }, std::string("f:gap_exp"));
      std::cout << "Min-gap (gap = diameter * (1 + max_topo_dist/e^topo_dist)) is " << min_gap << "\n";

      std::cout << "Computing possible min-gaps (gap = diameter * (1 + max_topo_dist/e^topo_dist) + similarity * mesh_size)\n";
      min_gap = HexMesher::determine_min_gap_direct(mesh, [&](HexMesher::FaceIndex idx) {
        if(topo_dist[idx] == 0)
        {
          return max_diameter;
        }

        double topo_penalty = max_topo_dist / std::exp(topo_dist[idx]);
        double similarity_penalty = std::pow(1.0 - similarity[idx], 0.2) * max_mesh_delta;
        return diameter[idx] * (1.0 + topo_penalty) + similarity_penalty;
      }, std::string("f:gap_exp_similarity"));
      std::cout << "Min-gap (gap = diameter * (1 + max_topo_dist/e^topo_dist) + similarity * mesh_size) is " << min_gap << "\n";

      std::cout << "Computing possible min-gaps (weighted score)\n";
      min_gap = HexMesher::determine_min_gap_weighted(mesh, [&](HexMesher::FaceIndex idx) {
        double topo_factor(0.0);

        if(topo_dist[idx] != 0)
        {
          topo_factor = topo_dist[idx] / max_topo_dist;
        }

        return (1.0 / 2.0) * topo_factor + (1.0 / 2.0) * similarity[idx];
      }, std::string("f:MIS_diameter"), std::string("f:gap_score"));
      std::cout << "Min-gap (weighted score) is " << min_gap << "\n";

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
      HexMesher::PolygonWithHoles2D& component = union_components[i];

      int idx;
      int total_vertices = component.outer_boundary().size();
      for(const HexMesher::Polygon2D& hole : component.holes())
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

      HexMesher::Polygon2D simplified_boundary = HexMesher::simplify_by_normal(component.outer_boundary(), pred).first;
      HexMesher::Polygon2D grid_sampled_boundary = HexMesher::grid_sample(component.outer_boundary(), h);

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
    //const std::string filename("/home/user/mmuegge/nobackup/repos/feat/data/models/scalexa_gendie_simple.off");
    //const std::string filename("/home/user/mmuegge/nobackup/projects/hexmesher/meshes/surface_22630.off");
    const std::string filename("/home/user/mmuegge/nobackup/projects/hexmesher/meshes/impeller_KSB.obj");

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

    HexMesher::compute_vertex_normals(mesh);
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
