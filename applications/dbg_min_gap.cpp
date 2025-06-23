#include <hexmesher.hpp>

#include <iostream>

static void print_min_gap(const HexMesher::MinGap& min_gap)
{
  std::cout << "Min-gap of " << min_gap.gap << " between faces " << min_gap.origin << " and " << min_gap.limiting
            << "\n";
  std::cout << "Use `SelectIDs(IDs=[0, " << min_gap.origin << ", 0, " << min_gap.limiting
            << "], FieldType='CELL')` to select the chosen triangles in ParaView\n";
}

int main(int argc, char* argv[])
{
  if(argc > 1)
  {
    const std::string filename(argv[1]);

    HexMesher::Result<HexMesher::SurfaceMesh, std::string> result = HexMesher::load_from_file(filename);

    if(result.is_err())
    {
      std::cout << "Reading mesh failed with error: " << result.err_ref() << "\n";
      exit(1);
    }

    HexMesher::SurfaceMesh mesh = std::move(result).take_ok();

    HexMesher::MinGap first_percentile = mesh.min_gap_percentile(0.01);
    std::cout << "1st percentile min gap by score: ";
    print_min_gap(first_percentile);
    std::cout << "\n";

    HexMesher::MinGap fifth_percentile = mesh.min_gap_percentile(0.05);
    std::cout << "5th percentile min gap by score: ";
    print_min_gap(fifth_percentile);
    std::cout << "\n";

    HexMesher::MinGap tenth_percentile = mesh.min_gap_percentile(0.1);
    std::cout << "10th percentile min gap by score: ";
    print_min_gap(tenth_percentile);
    std::cout << "\n";

    HexMesher::MinGap twentieth_percentile = mesh.min_gap_percentile(0.20);
    std::cout << "20th percentile min gap by score: ";
    print_min_gap(twentieth_percentile);
    std::cout << "\n";

    HexMesher::MinGap fiftieth_percentile = mesh.min_gap_percentile(0.50);
    std::cout << "50th percentile min gap by score: ";
    print_min_gap(fiftieth_percentile);
    std::cout << "\n";

    HexMesher::MinGap fixed = mesh.min_gap();
    std::cout << "0th percentile min gap by score: ";
    print_min_gap(fixed);
    std::cout << "\n";
  }
}
