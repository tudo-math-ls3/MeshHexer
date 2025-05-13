#pragma once

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_slicer.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>

#include <string>

namespace HexMesher
{
  using CGALKernel = CGAL::Exact_predicates_exact_constructions_kernel;
  using Mesh = CGAL::Surface_mesh<CGALKernel::Point_3>;

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

  void union_of_cross_sections(const Mesh& mesh, const CrossSectionSampler& sampler);
}
