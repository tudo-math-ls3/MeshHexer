#pragma once

#include <types.hpp>

namespace HexMesher
{
  struct CuttingPlane
  {
    Plane3D plane;

    Point3D origin;
    Vector3D x_axis;
    Vector3D y_axis;

    Point2D project(const Point3D&) const;
  };

  class CrossSectionSampler
  {
  public:
    virtual Point3D origin() const = 0;
    virtual int num_planes() const = 0;
    virtual CuttingPlane get_plane(int) const = 0;
    virtual Point2D project(Point3D p) const = 0;
    virtual CuttingPlane get_plane_through_vertex(Point3D p) const = 0;
  };

  class RadialCrossSectionSampler : public CrossSectionSampler
  {
    int _num_planes;

    Point3D _origin;

    Vector3D _up;
    Vector3D _u;
    Vector3D _v;

  public:

    RadialCrossSectionSampler(int num_planes, const Point3D& o, const Vector3D& normal, const Vector3D& up) :
      _num_planes(num_planes),
      _origin(o),
      _up(up),
      _u(normal),
      _v(CGAL::cross_product(normal, up))
    {
      if(CGAL::scalar_product(normal, up) > 1e-8)
      {
        throw std::invalid_argument("RadialCrossSectionSampler: normal and up must be orthogonal!");
      }

      _up = _up / CGAL::approximate_sqrt(_up.squared_length());
      _u = _u / CGAL::approximate_sqrt(_u.squared_length());
      _v = _v / CGAL::approximate_sqrt(_v.squared_length());
    }

    Point3D origin() const override;
    int num_planes() const override;
    CuttingPlane get_plane(int idx) const override;

    Point2D project(Point3D p) const override;
    CuttingPlane get_plane_through_vertex(Point3D p) const override;
  };

  class LineCrossSectionSampler : public CrossSectionSampler
  {
    int _num_planes;

    Point3D _start;
    Point3D _end;

    Vector3D _normal;
    Vector3D _x_axis;
    Vector3D _y_axis;

  public:

    LineCrossSectionSampler(int num_planes, const Point3D& start, const Point3D& end, const Vector3D& normal, const Vector3D& up) :
      _num_planes(num_planes),
      _start(start),
      _end(end),
      _normal(normal),
      _x_axis(CGAL::cross_product(normal, up)),
      _y_axis(up)
    {
      if(CGAL::scalar_product(normal, up) > 1e-8)
      {
        throw std::invalid_argument("LineCrossSectionSampler: normal and up must be orthogonal!");
      }

      _normal = _normal / CGAL::approximate_sqrt(_normal.squared_length());
      _x_axis = _x_axis / CGAL::approximate_sqrt(_x_axis.squared_length());
      _y_axis = _y_axis / CGAL::approximate_sqrt(_y_axis.squared_length());
    }

    Point3D origin() const override;
    int num_planes() const override;
    CuttingPlane get_plane(int idx) const override;

    Point2D project(Point3D p) const override;
    CuttingPlane get_plane_through_vertex(Point3D p) const override;
  };

  double angle(const Vector2D& a, const Vector2D& b);

  std::vector<PolygonWithHoles2D> union_of_cross_sections(const Mesh& mesh, const CrossSectionSampler& sampler);

  std::pair<Polygon2D, std::vector<int>> simplify_by_normal(
    const Polygon2D& polygon,
    const std::function<bool(const std::vector<Vector2D>&, const Vector2D&)>& continue_pred);

  Polygon2D grid_sample(const Polygon2D& polygon, Real min_dist);

  struct MinGap
  {
    /// Face at which the gap originates
    FaceIndex origin = FaceIndex(0);
    /// Face which limits the gap
    FaceIndex limiting = FaceIndex(0);

    /// Size of the gap
    Real gap = Real(0.0);
  };

  /**
   * \brief Compute a measure of mesh thickness at each face of a mesh.
   *
   * Determines the maximum inscribed sphere at each face of the mesh,
   * i.e. the largest possible sphere that touches the centroid of the mesh and any other point of the mesh,
   * without intersecting the mesh.
   *
   * The diameter of each sphere is made available in a "f:MIS_diameter" mesh property of type double.
   * The id of the _other_ primitive that is touched by the sphere is made available in a "f:MIS_id" mesh property.
   *
   * See the following for details on the algorithm:
   * Shrinking sphere: A parallel algorithm for computing the thickness of 3D objects
   * Masatomo Inui, Nobuyuki Umezu, Ryohei Shimane
   * COMPUTER-AIDED DESIGN & APPLICATIONS, 2016, VOL. 13, NO. 2, 199–207
   * http://dx.doi.org/10.1080/16864360.2015.1084186
   */
  void compute_mesh_thickness(Mesh& mesh);

  /**
   * \brief Calculates topological distances on the mesh
   *
   * \param mesh Mesh to calculate distances on
   * \param property Name of a mesh property of with type FaceIndex
   *
   * Calculates the smallest distance along the edges of the mesh between any face f
   * and the corresponding face property[f].
   */
  void topological_distances(Mesh& mesh, const std::string& property, double max_distance = 0.0);
  void topological_distances(Mesh& mesh, const std::string& targets_property, const std::string& max_distance_property);

  MinGap determine_min_gap_weighted(
    Mesh& mesh,
    std::function<double(FaceIndex)> weighting,
    const std::string& diameter_property,
    const std::string& id_property,
    const std::string& property = std::string("f:gap"));

  MinGap determine_min_gap_direct(
    Mesh& mesh,
    std::function<double(FaceIndex)> gap_calc,
    const std::string& id_property,
    const std::string& property = std::string("f:gap"));

  /**
   * \brief Computes vertex normals and makes them available as v:normals
   */
  void compute_vertex_normals(Mesh& mesh);

  void compute_curvature(Mesh& mesh);

  void compute_max_dihedral_angle(Mesh& mesh);

  /**
   * \brief Returns true if face normals calculated via a cross product point towards the "unbounded" side of the mesh
   *
   * We write "unbounded" because we do not assume that the mesh is watertight. There thus might not be a "bounded" side.angle
   * This function will work either way.
   */
  bool do_normals_point_outside(const Mesh& mesh);

  bool is_wound_consistently(const Mesh& mesh);

  Vector3D surface_normal(Mesh& mesh, FaceIndex f, Point3D point);

  void score_gaps(Mesh& mesh);
  MinGap select_min_gap(Mesh& mesh, double percentile);
  MinGap select_min_gap2(Mesh& mesh);
}
