#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace HexMesher
{
  struct MinGap
  {
    /// Face at which the gap originates
    std::uint32_t origin = 0;
    /// Face which limits the gap
    std::uint32_t limiting = 0;

    /// Size of the gap
    double gap = 0.0;
  };

  struct SelfIntersectionWarning
  {
    static constexpr std::string_view name = "self-intersection";

    std::uint32_t tri_a;
    std::uint32_t tri_b;

    SelfIntersectionWarning(std::uint32_t a, std::uint32_t b) : tri_a(a), tri_b(b) {}
  };

  struct DegenerateTriangleWarning
  {
    static constexpr std::string_view name = "degenerate-triangle";

    std::uint32_t idx;

    DegenerateTriangleWarning(std::uint32_t i) : idx(i) {}
  };

  struct AnisotropicTriangleWarning
  {
    static constexpr std::string_view name = "anisotropic-triangle";

    std::uint32_t idx;

    AnisotropicTriangleWarning(std::uint32_t i) : idx(i) {}
  };

  struct MeshWarnings
  {
    std::vector<SelfIntersectionWarning> self_intersections;
    std::vector<DegenerateTriangleWarning> degenerate_triangles;
    std::vector<AnisotropicTriangleWarning> anisotropic_triangles;
  };
}
