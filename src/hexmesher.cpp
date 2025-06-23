#include <hexmesher.hpp>

#include <cgal_types.hpp>
#include <properties.hpp>
#include <types.hpp>
#include <warnings.hpp>

#include <algorithm>
#include <iostream>
#include <math.h>
#include <optional>
#include <unistd.h>
#include <variant>

#include <CGAL/Bbox_3.h>

#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/locate.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

namespace HexMesher
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  static bool ends_with(const std::string& string, const std::string& ending)
  {
    if(ending.size() > string.size())
    {
      return false;
    }
    return std::equal(ending.rbegin(), ending.rend(), string.rbegin());
  }

  class SurfaceMesh::SurfaceMeshImpl
  {
    /// Surface mesh
    Mesh _mesh;

    /// Cached AABBTree
    std::optional<AABBTree> _aabb_tree;

  public:
    SurfaceMeshImpl(Mesh&& m) : _mesh(std::move(m))
    {
    }

    AABBTree& aabb_tree()
    {
      if(!_aabb_tree.has_value())
      {
        _aabb_tree = AABBTree(_mesh.faces_begin(), _mesh.faces_end(), _mesh);
      }
      return _aabb_tree.value();
    }

    MinGap min_gap();
    MinGap min_gap_percentile(double percentile);

    Result<void, std::string> write_to_file(const std::string& filename);

    std::uint32_t num_vertices() const;
    std::uint32_t num_edges() const;
    std::uint32_t num_faces() const;

    bool is_closed() const;
    bool is_wound_consistently() const;
    bool is_outward_oriented() const;
    double minimal_aspect_ratio() const;
    double maximal_aspect_ratio() const;

    void warnings(MeshWarnings&) const;

  private:
    /**
     * \brief Ensures the given property exists on the mesh
     *
     * Does nothing if the property already exists.
     * Otherwise it calls the given function.
     */
    template<typename Index, typename T, typename Fn>
    void ensure_property(const std::string& property_name, Fn&& fn)
    {
      auto prop = _mesh.add_property_map<Index, T>(property_name, T{});

      if(prop.second)
      {
        fn();
      }
    }

    void prepare_for_min_gap();
  };

  void SurfaceMesh::SurfaceMeshImpl::prepare_for_min_gap()
  {
    // Ensure vertex normals are available
    ensure_property<VertexIndex, Vector3D>("v:normals", [&]() { compute_vertex_normals(_mesh); });

    // Ensure maximal dihedral angles are available
    ensure_property<FaceIndex, double>("f:dihedral_angle", [&]() { compute_max_dihedral_angle(_mesh); });

    // Ensure maximal inscribed spheres are available
    ensure_property<FaceIndex, double>("f:MIS_diameter", [&]() { maximal_inscribed_spheres(_mesh, aabb_tree()); });

    // Ensure topological distances are available
    ensure_property<FaceIndex, double>(
      "f:topological_distance",
      [&]() { topological_distances(_mesh, "f:MIS_id", "f:MIS_diameter"); });

    // Ensure gap scores are available
    ensure_property<FaceIndex, double>("f:gap_score", [&]() { score_gaps(_mesh); });
  }

  MinGap SurfaceMesh::SurfaceMeshImpl::min_gap()
  {
    prepare_for_min_gap();
    return HexMesher::min_gap(_mesh);
  }

  MinGap SurfaceMesh::SurfaceMeshImpl::min_gap_percentile(double percentile)
  {
    prepare_for_min_gap();
    return HexMesher::min_gap_percentile(_mesh, percentile);
  }

  Result<void, std::string> SurfaceMesh::SurfaceMeshImpl::write_to_file(const std::string& filename)
  {
    using ResultType = Result<void, std::string>;

    if(!ends_with(filename, ".ply"))
    {
      return ResultType::err("Can only write .ply files");
    }

    std::ofstream output(filename);

    if(!output)
    {
      return ResultType::err("Failed to open file for writing");
    }

    if(!CGAL::IO::write_PLY(output, _mesh))
    {
      return ResultType::err("Failed to write mesh to opened file");
    }

    return ResultType();
  }

  std::uint32_t SurfaceMesh::SurfaceMeshImpl::num_vertices() const
  {
    return _mesh.num_vertices();
  }

  std::uint32_t SurfaceMesh::SurfaceMeshImpl::num_edges() const
  {
    return _mesh.num_edges();
  }

  std::uint32_t SurfaceMesh::SurfaceMeshImpl::num_faces() const
  {
    return _mesh.num_faces();
  }

  bool SurfaceMesh::SurfaceMeshImpl::is_closed() const
  {
    return CGAL::is_closed(_mesh);
  }

  bool SurfaceMesh::SurfaceMeshImpl::is_wound_consistently() const
  {
    return HexMesher::is_wound_consistently(_mesh);
  }

  bool SurfaceMesh::SurfaceMeshImpl::is_outward_oriented() const
  {
    return CGAL::Polygon_mesh_processing::is_outward_oriented(_mesh);
  }

  double SurfaceMesh::SurfaceMeshImpl::minimal_aspect_ratio() const
  {
    double min_aspect_ratio = std::numeric_limits<double>::max();
    for(HexMesher::FaceIndex f : _mesh.faces())
    {
      double ratio = CGAL::Polygon_mesh_processing::face_aspect_ratio(f, _mesh);
      min_aspect_ratio = std::min(min_aspect_ratio, ratio);
    }

    return min_aspect_ratio;
  }

  double SurfaceMesh::SurfaceMeshImpl::maximal_aspect_ratio() const
  {
    double max_aspect_ratio = 0;
    for(HexMesher::FaceIndex f : _mesh.faces())
    {
      double ratio = CGAL::Polygon_mesh_processing::face_aspect_ratio(f, _mesh);
      max_aspect_ratio = std::max(max_aspect_ratio, ratio);
    }

    return max_aspect_ratio;
  }

  void SurfaceMesh::SurfaceMeshImpl::warnings(MeshWarnings& ws) const
  {
    ws.self_intersections.clear();
    ws.degenerate_triangles.clear();
    ws.anisotropic_triangles.clear();

    create_warnings(_mesh, ws);
  }

  SurfaceMesh::SurfaceMesh(std::unique_ptr<SurfaceMesh::SurfaceMeshImpl> ptr) : impl(std::move(ptr))
  {
  }

  SurfaceMesh::~SurfaceMesh() = default;
  SurfaceMesh::SurfaceMesh(SurfaceMesh&&) = default;
  SurfaceMesh& SurfaceMesh::operator=(SurfaceMesh&&) = default;

  MinGap SurfaceMesh::min_gap()
  {
    return impl->min_gap();
  }

  MinGap SurfaceMesh::min_gap_percentile(double percentile)
  {
    return impl->min_gap_percentile(percentile);
  }

  Result<void, std::string> SurfaceMesh::write_to_file(const std::string& filename)
  {
    return impl->write_to_file(filename);
  }

  std::uint32_t SurfaceMesh::num_vertices() const
  {
    return impl->num_vertices();
  }

  std::uint32_t SurfaceMesh::num_edges() const
  {
    return impl->num_edges();
  }

  std::uint32_t SurfaceMesh::num_faces() const
  {
    return impl->num_faces();
  }

  bool SurfaceMesh::is_closed() const
  {
    return impl->is_closed();
  }

  bool SurfaceMesh::is_wound_consistently() const
  {
    return impl->is_wound_consistently();
  }

  bool SurfaceMesh::is_outward_oriented() const
  {
    return impl->is_outward_oriented();
  }

  double SurfaceMesh::minimal_aspect_ratio() const
  {
    return impl->minimal_aspect_ratio();
  }

  double SurfaceMesh::maximal_aspect_ratio() const
  {
    return impl->maximal_aspect_ratio();
  }

  void SurfaceMesh::warnings(MeshWarnings& ws) const
  {
    impl->warnings(ws);
  }

  Result<SurfaceMesh, std::string> load_from_file(const std::string& filename, bool triangulate)
  {
    using ResultType = Result<SurfaceMesh, std::string>;

    HexMesher::Mesh mesh;
    if(ends_with(filename, ".ply"))
    {
      std::ifstream mesh_file(filename);
      std::string comment("");
      if(!CGAL::IO::read_PLY(mesh_file, mesh, comment, true))
      {
        return ResultType::err("Failed to read mesh " + filename);
      }
    }
    else if(!PMP::IO::read_polygon_mesh(filename, mesh))
    {
      return ResultType::err("Failed to read mesh " + filename);
    }

    if(CGAL::is_empty(mesh))
    {
      return ResultType::err("Mesh " + filename + " is empty.");
    }

    bool is_triangle_mesh = CGAL::is_triangle_mesh(mesh);
    if(!is_triangle_mesh && triangulate)
    {
      PMP::triangulate_faces(mesh);
    }
    else if(!is_triangle_mesh && !triangulate)
    {
      return ResultType::err("Mesh " + filename + " is not a triangle mesh.");
    }

    auto impl = std::make_unique<SurfaceMesh::SurfaceMeshImpl>(std::move(mesh));
    SurfaceMesh smesh(std::move(impl));

    return ResultType::ok(std::move(smesh));
  }
} // namespace HexMesher
