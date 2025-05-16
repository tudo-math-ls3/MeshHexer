#include <hexmesher.hpp>
#include <iostream>

int main(int argc, char* argv[])
{
  std::cout << "Running hexmesher\n";



  if(argc > 11)
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

    if(mode == "radial")
    {
      HexMesher::Point p(std::atof(argv[3]), std::atof(argv[4]), std::atof(argv[5]));
      HexMesher::Vector up(std::atof(argv[6]), std::atof(argv[7]), std::atof(argv[8]));
      HexMesher::Vector u(std::atof(argv[9]), std::atof(argv[10]), std::atof(argv[11]));

      int steps = std::atoi(argv[12]);

      HexMesher::RadialCrossSectionSampler sampler(steps, p, u, up);

      HexMesher::union_of_cross_sections(mesh, sampler);
    }

    if(mode == "line")
    {
      HexMesher::Point start(std::atof(argv[3]), std::atof(argv[4]), std::atof(argv[5]));
      HexMesher::Point end(std::atof(argv[6]), std::atof(argv[7]), std::atof(argv[8]));
      HexMesher::Vector normal(std::atof(argv[9]), std::atof(argv[10]), std::atof(argv[11]));
      HexMesher::Vector up(std::atof(argv[12]), std::atof(argv[13]), std::atof(argv[14]));

      int steps = std::atoi(argv[15]);

      HexMesher::LineCrossSectionSampler sampler(steps, start, end, normal, up);

      HexMesher::union_of_cross_sections(mesh, sampler);
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

    //HexMesher::Point start(0.0, 0.0, 80.0);
    //HexMesher::Point end(0.0, 0.0, 60.0);
    //HexMesher::Vector normal(0.0, 0.0, 1.0);
    //HexMesher::Vector up(0.0, 1.0, 0.0);

    HexMesher::Point origin(0.0, 0.0, 0.0);
    HexMesher::Vector up(0.0, 0.0, 1.0);
    HexMesher::Vector normal(1.0, 0.0, 0.0);

    //HexMesher::LineCrossSectionSampler sampler(50, start, end, normal, up);
    HexMesher::RadialCrossSectionSampler sampler(50, origin, normal, up);

    HexMesher::union_of_cross_sections(mesh, sampler);
  }

  return 0;
}
