#pragma once

#include <types.hpp>

namespace HexMesher
{
  void write_geo_compound_2d(const std::string& filename, const Polygon& poly);
  void write_geo(const std::string& filename, const Polygon& poly);
  void write_polygon_brep(const std::string& filename, const Polygon& poly);
  void write_polygon(const std::string& filename, const Polygon& poly);
  void write_polygon(const std::string& filename, const PolygonWithHoles& poly);
  void write_polylines(const std::string& filename, const Polylines& polylines);
  void write_polylines(const std::string& filename, const Polylines2D& polylines);
}
