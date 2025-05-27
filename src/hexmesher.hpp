#pragma once

#include <types.hpp>

namespace HexMesher
{
  using CGALKernel = CGAL::Exact_predicates_exact_constructions_kernel;

  using Mesh = CGAL::Surface_mesh<CGALKernel::Point_3>;

  using VertexIndex = Mesh::Vertex_index;
  using EdgeIndex = Mesh::Edge_index;
  using FaceIndex = Mesh::Face_index;

  using Real = CGALKernel::RT;

  using Point = CGALKernel::Point_3;
  using Vector = CGALKernel::Vector_3;
  using Plane = CGALKernel::Plane_3;

  using Polyline = std::vector<Point>;
  using Polylines = std::list<Polyline>;

  using Point2D = CGALKernel::Point_2;
  using Vector2D = CGALKernel::Vector_2;

  using Polyline2D = std::vector<Point2D>;
  using Polylines2D = std::list<Polyline2D>;

  using Polygon = CGAL::Polygon_2<CGALKernel>;
  using PolygonWithHoles = CGAL::Polygon_with_holes_2<CGALKernel>;

  struct CuttingPlane
  {
    Plane plane;

    Point origin;
    Vector x_axis;
    Vector y_axis;

    Point2D project(const Point&) const;
  };

  class CrossSectionSampler
  {
  public:
    virtual Point origin() const = 0;
    virtual int num_planes() const = 0;
    virtual CuttingPlane get_plane(int) const = 0;
    virtual Point2D project(Point p) const = 0;
    virtual CuttingPlane get_plane_through_vertex(Point p) const = 0;
  };

  class RadialCrossSectionSampler : public CrossSectionSampler
  {
    int _num_planes;

    Point _origin;

    Vector _up;
    Vector _u;
    Vector _v;

  public:

    RadialCrossSectionSampler(int num_planes, const Point& o, const Vector& normal, const Vector& up) :
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

    Point origin() const override;
    int num_planes() const override;
    CuttingPlane get_plane(int idx) const override;

    Point2D project(Point p) const override;
    CuttingPlane get_plane_through_vertex(Point p) const override;
  };

  class LineCrossSectionSampler : public CrossSectionSampler
  {
    int _num_planes;

    Point _start;
    Point _end;

    Vector _normal;
    Vector _x_axis;
    Vector _y_axis;

  public:

    LineCrossSectionSampler(int num_planes, const Point& start, const Point& end, const Vector& normal, const Vector& up) :
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

    Point origin() const override;
    int num_planes() const override;
    CuttingPlane get_plane(int idx) const override;

    Point2D project(Point p) const override;
    CuttingPlane get_plane_through_vertex(Point p) const override;
  };

  double angle(const Vector2D& a, const Vector2D& b);

  std::vector<PolygonWithHoles> union_of_cross_sections(const Mesh& mesh, const CrossSectionSampler& sampler);

  std::pair<Polygon, std::vector<int>> simplify_by_normal(
    const Polygon& polygon,
    const std::function<bool(const std::vector<Vector2D>&, const Vector2D&)>& continue_pred);

  Polygon grid_sample(const Polygon& polygon, Real min_dist);

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
   * \brief Returns the smallest distance along the edges of the mesh between any two vertices of the given faces
   */
  double topological_distance(FaceIndex a, FaceIndex b, const Mesh& mesh);

  /**
   * \brief Calculates topological distances on the mesh
   * 
   * \param mesh Mesh to calculate distances on
   * \param property Name of a mesh property of with type FaceIndex
   * 
   * Calculates the smallest distance along the edges of the mesh between any face f
   * and the corresponding face property[f].
   */
  void topological_distances(Mesh& mesh, const std::string& property);

  double determine_min_gap_weighted(Mesh& mesh, std::function<double(FaceIndex)> weighting, const std::string& diameter_property, const std::string& property = std::string("f:gap"));
  double determine_min_gap_direct(Mesh& mesh, std::function<double(FaceIndex)> gap_calc, const std::string& property = std::string("f:gap"));

  /**
   * \brief Computes vertex normals and makes them available as v:normals
   */
  void compute_vertex_normals(Mesh& mesh);

  Vector surface_normal(Mesh& mesh, FaceIndex f, Point point);
}
