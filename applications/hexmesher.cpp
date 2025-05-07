#include <hexmesher.hpp>
#include <iostream>

int main(int argc, char* argv[])
{
  std::cout << "Running hexmesher\n";

  if(argc > 10)
  {
    const std::string filename(argv[1]);

    HexMesher::Point p(std::atof(argv[2]), std::atof(argv[3]), std::atof(argv[4]));
    HexMesher::Vector up(std::atof(argv[5]), std::atof(argv[6]), std::atof(argv[7]));
    HexMesher::Vector u(std::atof(argv[8]), std::atof(argv[9]), std::atof(argv[10]));

    HexMesher::union_of_cross_sections(filename, p, up, u, 50);
  }
  else
  {
    const std::string filename("/home/user/mmuegge/nobackup/repos/feat/data/models/scalexa_gendie_small.off");
    HexMesher::Point p(0.0, 0.0, 0.0);
    HexMesher::Vector up(0.0, 0.0, 1.0);
    HexMesher::Vector u(1.0, 0.0, 0.0);

    HexMesher::union_of_cross_sections(filename, p, up, u, 50);
  }

  return 0;
}
