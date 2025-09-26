#pragma once

#include <cgal_types.hpp>
#include <meshhexer/types.hpp>

namespace MeshHexer
{
  void write_geo_compound_2d(const std::string& filename, const Polygon2D& poly);
  void write_geo(const std::string& filename, const Polygon2D& poly);
  void write_polygon_brep(const std::string& filename, const Polygon2D& poly);
  void write_polygon(const std::string& filename, const Polygon2D& poly);
  void write_polygon(const std::string& filename, const PolygonWithHoles2D& poly);
  void write_polylines(const std::string& filename, const Polylines3D& polylines);
  void write_polylines(const std::string& filename, const Polylines2D& polylines);

  void write_feat_xml(std::ostream& stream, const VolumeMesh& vmesh);
} // namespace MeshHexer
