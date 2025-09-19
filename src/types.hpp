#pragma once

#include <array>
#include <cstdint>
#include <iostream>
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

  struct Point
  {
    double x;
    double y;
    double z;
  };

  struct BoundingBox
  {
    Point min;
    Point max;
  };

  struct Slice
  {
    double coord = 0.0;
    std::uint64_t subdivision_level = 0;

    Slice() = default;

    Slice(double c, std::uint64_t lvl) : coord(c), subdivision_level(lvl)
    {
    }
  };

  class VolumeMesh
  {
  public:
    using Iterator = std::vector<Slice>::iterator;
    using ConstIterator = std::vector<Slice>::const_iterator;

  private:
    // Grid points, including start and end points.
    std::vector<Slice> _xs;
    std::vector<Slice> _ys;
    std::vector<Slice> _zs;

    std::array<std::size_t, 3> _num_slices;
    std::array<std::size_t, 4> _num_elements;

  public:
    VolumeMesh(std::size_t nx, std::size_t ny, std::size_t nz) :
      _xs(nx),
      _ys(ny),
      _zs(nz),
      _num_slices({nx - 1, ny - 1, nz - 1}),
      _num_elements()
    {
      _num_elements[0] = nx * ny * nz;
      _num_elements[1] = (nx - 1) * ny * nz + nx * (ny - 1) * nz + nx * ny * (nz - 1);
      _num_elements[2] = (nx - 1) * (ny - 1) * nz + nx * (ny - 1) * (nz - 1) + (nx - 1) * ny * (nz - 1);
      _num_elements[3] = (nx - 1) * (ny - 1) * (nz - 1);
    }

    template<typename SliceIter>
    VolumeMesh(
      SliceIter x_begin,
      SliceIter x_end,
      SliceIter y_begin,
      SliceIter y_end,
      SliceIter z_begin,
      SliceIter z_end) :
      _xs(x_begin, x_end),
      _ys(y_begin, y_end),
      _zs(z_begin, z_end),
      _num_slices({_xs.size() - 1, _ys.size() - 1, _zs.size() - 1}),
      _num_elements()
    {
      std::size_t nx = _xs.size();
      std::size_t ny = _ys.size();
      std::size_t nz = _zs.size();

      _num_elements[0] = nx * ny * nz;
      _num_elements[1] = (nx - 1) * ny * nz + nx * (ny - 1) * nz + nx * ny * (nz - 1);
      _num_elements[2] = (nx - 1) * (ny - 1) * nz + nx * (ny - 1) * (nz - 1) + (nx - 1) * ny * (nz - 1);
      _num_elements[3] = (nx - 1) * (ny - 1) * (nz - 1);
    }

    const std::size_t* num_slices() const
    {
      return _num_slices.data();
    }

    const std::size_t* num_elements() const
    {
      return _num_elements.data();
    }

    Iterator xs_begin()
    {
      return _xs.begin();
    }

    Iterator xs_end()
    {
      return _xs.end();
    }

    Iterator ys_begin()
    {
      return _ys.begin();
    }

    Iterator ys_end()
    {
      return _ys.end();
    }

    Iterator zs_begin()
    {
      return _zs.begin();
    }

    Iterator zs_end()
    {
      return _zs.end();
    }

    std::size_t num_vertices() const
    {
      return _xs.size() * _ys.size() * _zs.size();
    }

    std::size_t num_edges() const
    {
      return ((_xs.size() - 1) * _ys.size() * _zs.size()) + ((_ys.size() - 1) * _xs.size() * _zs.size()) +
             ((_zs.size() - 1) * _xs.size() * _ys.size());
    }

    std::size_t num_faces() const
    {
      return (_xs.size() * (_ys.size() - 1) * (_zs.size() - 1)) + (_ys.size() * (_xs.size() - 1) * (_zs.size() - 1)) +
             (_zs.size() * (_xs.size() - 1) * (_ys.size() - 1));
    }

    std::size_t num_cells() const
    {
      return (_xs.size() - 1) * (_ys.size() - 1) * (_zs.size() - 1);
    }

    HexMesher::Point vertex(std::size_t idx) const
    {
      std::size_t idx_z = idx / (_xs.size() * _ys.size());
      std::size_t idx_y = idx % (_xs.size() * _ys.size()) / _xs.size();
      std::size_t idx_x = idx % (_xs.size() * _ys.size()) % _xs.size();

      return HexMesher::Point{_xs[idx_x].coord, _ys[idx_y].coord, _zs[idx_z].coord};
    }

    std::size_t edge(std::size_t i, std::size_t j) const
    {
      // NOTE(mmuegge): Copied from FEAT3 (kernel/geometry/intern/struct_index_mapping.hpp)
      // FEATs structured mesh counts fences, rather than fence posts.
      // Hence the +- 1 everywhere

      const std::array<std::size_t, 3> num_slices = {_xs.size() - 1, _ys.size() - 1, _zs.size() - 1};

      // auxiliary index
      std::size_t pos = 0;

      // z-direction
      if(
        i >= num_slices[1] * (num_slices[0] + 1) * (num_slices[2] + 1) +
               (num_slices[2] + 1) * (num_slices[1] + 1) * num_slices[0])
      {
        pos = i - num_slices[1] * (num_slices[0] + 1) * (num_slices[2] + 1) -
              (num_slices[2] + 1) * (num_slices[1] + 1) * num_slices[0];
        return (pos / num_slices[2]) + ((pos % num_slices[2] + j) * (num_slices[0] + 1) * (num_slices[1] + 1));
      }
      // y-direction
      else if(i >= (num_slices[2] + 1) * (num_slices[1] + 1) * num_slices[0])
      {
        pos = i - (num_slices[2] + 1) * (num_slices[1] + 1) * num_slices[0];
        std::size_t x = (pos / num_slices[1]) % (num_slices[0] + 1);
        std::size_t y = pos % num_slices[1];
        std::size_t z = (pos / num_slices[1]) / (num_slices[0] + 1);
        return (z * (num_slices[0] + 1) * (num_slices[1] + 1)) + ((y + j) * (num_slices[0] + 1)) + x;
      }
      // x-direction
      else
      {
        return i + j + (i / num_slices[0]);
      }
    }

    std::size_t face(std::size_t i, std::size_t j) const
    {
      // NOTE(mmuegge): Copied from FEAT3 (kernel/geometry/intern/struct_index_mapping.hpp)
      // FEATs structured mesh counts fences, rather than fence posts.
      // Hence the +- 1 everywhere

      const std::array<std::size_t, 3> num_slices = {_xs.size() - 1, _ys.size() - 1, _zs.size() - 1};

      std::size_t pos = 0;
      if(i < num_slices[0] * num_slices[1] * (num_slices[2] + 1))
      {
        return i + (i / num_slices[0]) + ((i / (num_slices[0] * num_slices[1])) * (num_slices[0] + 1)) + (j % 2) +
               ((j / 2) * (num_slices[0] + 1));
      }
      else if(
        i >= num_slices[0] * num_slices[1] * (num_slices[2] + 1) &&
        i < num_slices[0] * (num_slices[1] * (num_slices[2] + 1) + num_slices[2] * (num_slices[1] + 1)))
      {
        pos = i - num_slices[0] * num_slices[1] * (num_slices[2] + 1);
        return (pos % num_slices[0]) + ((pos / (num_slices[0] * num_slices[2])) * (num_slices[0] + 1)) +
               ((pos % (num_slices[0] * num_slices[2])) / num_slices[0] * (num_slices[0] + 1) * (num_slices[1] + 1)) +
               (j % 2) + ((j / 2) * (num_slices[0] + 1) * (num_slices[1] + 1));
      }
      else
      {
        pos = i - num_slices[0] * (num_slices[1] * (num_slices[2] + 1) + num_slices[2] * (num_slices[1] + 1));
        return (pos / (num_slices[1] * num_slices[2])) + ((pos % num_slices[1]) * (num_slices[0] + 1)) +
               (((pos / num_slices[1]) % num_slices[2]) * (num_slices[0] + 1) * (num_slices[1] + 1)) +
               ((j % 2) * (num_slices[0] + 1)) + ((j / 2) * (num_slices[0] + 1) * (num_slices[1] + 1));
      }
    }

    std::size_t cell(std::size_t idx, std::size_t vert) const
    {
      const std::size_t z_layer = idx / ((_xs.size() - 1) * (_ys.size() - 1));
      const std::size_t y_layer = idx % ((_xs.size() - 1) * (_ys.size() - 1)) / (_xs.size() - 1);
      const std::size_t x_layer = idx % ((_xs.size() - 1) * (_ys.size() - 1)) % (_xs.size() - 1);

      const std::size_t base = (z_layer * (_xs.size() * _ys.size())) + (y_layer * _xs.size()) + x_layer;
      const std::size_t z_offset = _xs.size() * _ys.size();

      switch(vert)
      {
      case 0:
        return base;
      case 1:
        return base + 1;
      case 2:
        return base + _xs.size();
      case 3:
        return base + _xs.size() + 1;
      case 4:
        return base + z_offset;
      case 5:
        return base + 1 + z_offset;
      case 6:
        return base + _xs.size() + z_offset;
      case 7:
        return base + _xs.size() + 1 + z_offset;
      default:
        std::abort();
      }
    }

    std::uint64_t subdivision_level(std::size_t idx)
    {
      Point vtx = vertex(idx);

      std::uint64_t level = 0;

      for(std::size_t slice(0); slice < _num_slices[2]; slice++)
      {
        if(_zs[slice].coord <= vtx.z && vtx.z <= _zs[slice + 1].coord)
        {
          level = std::max(_zs[slice].subdivision_level, level);
        }
      }

      return level;
    }

    void write_feat_xml(std::ostream& stream) const
    {
      stream << "<FeatMeshFile version=\"1\" mesh=\"conformal:hypercube:3:3\">\n";
      stream << "<Mesh type=\"conformal:hypercube:3:3\" size=\""
      << num_vertices() << " "
      << num_edges() << " "
      << num_faces() << " "
      << num_cells() << "\">\n";
      stream << "<Vertices>\n";
      for(std::size_t i(0); i < num_vertices(); i++)
      {
        const Point point = vertex(i);
        stream << point.x << " " << point.y << " " << point.z << "\n";
      }
      stream << "</Vertices>\n";

      stream << "<Topology dim=\"1\">\n";
      for(std::size_t i(0); i < num_edges(); i++)
      {
        for(std::size_t j(0); j < 2; j++)
        {
          stream << edge(i, j);
          stream << (j == 1 ? "\n" : " ");
        }
      }
      stream << "</Topology>\n";

      stream << "<Topology dim=\"2\">\n";
      for(std::size_t i(0); i < num_faces(); i++)
      {
        for(std::size_t j(0); j < 4; j++)
        {
          stream << face(i, j);
          stream << (j == 3 ? "\n" : " ");
        }
      }
      stream << "</Topology>\n";

      stream << "<Topology dim=\"3\">\n";
      for(std::size_t i(0); i < num_cells(); i++)
      {
        for(std::size_t j(0); j < 8; j++)
        {
          stream << cell(i, j);
          stream << (j == 7 ? "\n" : " ");
        }
      }
      stream << "</Topology>\n";

      stream << "</Mesh>\n";
      stream << "</FeatMeshFile>\n";
    }
  };

  struct SelfIntersectionWarning
  {
    static constexpr std::string_view name = "self-intersection";

    std::uint32_t tri_a;
    std::uint32_t tri_b;

    SelfIntersectionWarning(std::uint32_t a, std::uint32_t b) : tri_a(a), tri_b(b)
    {
    }
  };

  struct DegenerateTriangleWarning
  {
    static constexpr std::string_view name = "degenerate-triangle";

    std::uint32_t idx;

    explicit DegenerateTriangleWarning(std::uint32_t i) : idx(i)
    {
    }
  };

  struct AnisotropicTriangleWarning
  {
    static constexpr std::string_view name = "anisotropic-triangle";

    std::uint32_t idx;

    explicit AnisotropicTriangleWarning(std::uint32_t i) : idx(i)
    {
    }
  };

  struct MeshWarnings
  {
    std::vector<SelfIntersectionWarning> self_intersections;
    std::vector<DegenerateTriangleWarning> degenerate_triangles;
    std::vector<AnisotropicTriangleWarning> anisotropic_triangles;
  };
} // namespace HexMesher
