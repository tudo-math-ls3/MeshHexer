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

  using Polyline = std::vector<Point>;
  using Polylines = std::list<Polyline>;

  using Point2D = CGALKernel::Point_2;
  using Vector2D = CGALKernel::Vector_2;

  using Polyline2D = std::vector<Point2D>;
  using Polylines2D = std::list<Polyline2D>;

  using Polygon = CGAL::Polygon_2<CGALKernel>;
  using PolygonWithHoles = CGAL::Polygon_with_holes_2<CGALKernel>;

  void union_of_cross_sections(const std::string& filename, Point p, Vector up, Vector u, int steps);
}
