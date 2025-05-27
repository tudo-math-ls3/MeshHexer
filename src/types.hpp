#pragma once

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>

namespace HexMesher
{
  // Exact types
  using ExactKernel = CGAL::Exact_predicates_exact_constructions_kernel;

  // Inexact types

  using Kernel = CGAL::Simple_cartesian<double>;

  using Point3D = Kernel::Point_3;
  using Vector3D = Kernel::Vector_3;
  using Plane3D = Kernel::Plane_3;

  using Point2D = Kernel::Point_2;
  using Vector2D = Kernel::Vector_2;

  using Mesh = CGAL::Surface_mesh<Point3D>;
  using VertexIndex = Mesh::Vertex_index;
  using EdgeIndex = Mesh::Edge_index;
  using FaceIndex = Mesh::Face_index;
  using HalfedgeIndex = Mesh::Halfedge_index;

  using Real = Kernel::RT;

  using Polyline3D = std::vector<Point3D>;
  using Polylines3D = std::list<Polyline3D>;

  using Polyline2D = std::vector<Point2D>;
  using Polylines2D = std::list<Polyline2D>;

  using Polygon2D = CGAL::Polygon_2<Kernel>;
  using PolygonWithHoles2D = CGAL::Polygon_with_holes_2<Kernel>;
}
