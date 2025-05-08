#include <hexmesher.hpp>

#include <iostream>
#include <math.h>

#include <CGAL/Boolean_set_operations_2.h>

namespace HexMesher
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  PolygonWithHoles make_polygon(Polylines2D& polylines)
  {
    // Create polygons from the polylines.
    // This gives us inside-out tests without
    // without having to think about the orientation
    // of the polyline
    std::vector<Polygon> polygons;

    for(auto it = polylines.begin(); it != polylines.end(); it++)
    {
      Polygon poly(it->begin(), it->end() - 1);
      if(poly.is_clockwise_oriented())
      {
        polygons.emplace_back(std::reverse_iterator(it->end() - 1), std::reverse_iterator(it->begin()));
      }
      else
      {
        polygons.emplace_back(std::move(poly));
      }
    }

    // We assume that there is a single closed polyline that contains all others
    // Find it and save it as the boundary polyline.
    Polygon* outer = &polygons[0];
    int outer_idx = 0;

    for(int i(1); i < polygons.size(); i++)
    {
      Polygon& candidate = polygons[i];
      for(auto it = candidate.vertices_begin(); it != candidate.vertices_end(); it++)
      {
        if(outer->bounded_side(*it) == CGAL::ON_UNBOUNDED_SIDE)
        {
          // Candidate polygon has at least one vertex outside the current outer polygon
          outer = &polygons[i];
          outer_idx = i;
          break;
        }
      }
    }

    std::vector<Polygon> holes;

    auto it = polygons.begin();
    for(int i(0); i < polygons.size(); i++)
    {
      if(i != outer_idx)
      {
        holes.push_back(*it);
      }
      it++;
    }

    return PolygonWithHoles(*outer, holes.begin(), holes.end());
  }

  PolygonWithHoles merge(PolygonWithHoles& a, PolygonWithHoles& b)
  {
    PolygonWithHoles new_boundary;
    CGAL::join(a.outer_boundary(), b.outer_boundary(), new_boundary);

    std::vector<PolygonWithHoles> intersections;

    for(auto& hole_a : a.holes())
    {
      for(auto& hole_b : b.holes())
      {
        CGAL::intersection(hole_a, hole_b, std::back_inserter(intersections));
      }
    }

    std::vector<Polygon> new_holes;

    for(auto& intersection : intersections)
    {
      new_holes.push_back(intersection.outer_boundary());
    }

    return PolygonWithHoles(new_boundary.outer_boundary(), new_holes.begin(), new_holes.end());
  }

  void write_polygon(const std::string& filename, const Polygon& poly)
  {
    std::ofstream output(filename);

    if (output) {
      output<< "<VTKFile type=\"PolyData\">\n";
      output << "<PolyData>\n";
      output << "<Piece NumberOfPoints=\"" << poly.size() << "\" NumberOfVerts=\"0\" NumberOfLines=\"1\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n";
      output << "<PointData></PointData>\n";
      output << "<CellData></CellData>\n";
      output << "<Points>\n";
      output << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">\n";

      for(const Point2D& p : poly.vertices())
      {
        output << p.x() << " " << p.y() << " 0\n";
      }

      output << "</DataArray>\n";
      output << "</Points>\n";
      output << "<Verts></Verts>\n";
      output << "<Lines>\n";
      output << "<DataArray type=\"UInt32\" Name=\"connectivity\">\n";

      for(int i(0); i < poly.size(); i++)
      {
        output << i << " ";
      }
      output << "0\n";
      output << "</DataArray>\n";

      output << "<DataArray type=\"UInt32\" Name=\"offsets\">\n";
      output << poly.size() << "\n";
      output << "</DataArray>\n";

      output << "</Lines>\n";
      output << "<Strips></Strips>\n";
      output << "<Polys></Polys>\n";
      output << "</Piece>\n";
      output << "</PolyData>\n";
        output<< "</VTKFile>\n";
      }
      else {
        std::cerr << "Failed to write polygon to " << filename << "!\n";
      }
  }

  void write_polygon(const std::string& filename, const PolygonWithHoles& poly)
  {
    std::ofstream output(filename);

    if (output) {
      int num_polys = 1 + poly.holes().size();
      int total_points = 0;
      std::vector<int> starting_points;

      starting_points.push_back(total_points);
      total_points += poly.outer_boundary().size();
      for(const Polygon& hole : poly.holes())
      {
        starting_points.push_back(total_points);
        total_points += hole.size();
      }

      output<< "<VTKFile type=\"PolyData\">\n";
      output << "<PolyData>\n";
      output << "<Piece NumberOfPoints=\"" << total_points << "\" NumberOfVerts=\"0\" NumberOfLines=\"" << num_polys << "\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n";
      output << "<PointData></PointData>\n";
      output << "<CellData></CellData>\n";
      output << "<Points>\n";
      output << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">\n";

      for(const Point2D& p : poly.outer_boundary().vertices())
      {
        output << p.x() << " " << p.y() << " 0\n";
      }

      for(const Polygon& hole : poly.holes())
      {
        for(const Point2D& p : hole.vertices())
        {
          output << p.x() << " " << p.y() << " 0\n";
        }
      }

      output << "</DataArray>\n";
      output << "</Points>\n";
      output << "<Verts></Verts>\n";
      output << "<Lines>\n";
      output << "<DataArray type=\"UInt32\" Name=\"connectivity\">\n";

      if(poly.outer_boundary().size() > 0)
      {
        int idx = 0;
        for(int i(0); i < poly.outer_boundary().size() - 1; i++)
        {
          output << i << " ";
        }
        output << starting_points[idx] << "\n";
        idx++;

        for(const Polygon& hole : poly.holes())
        {
          int starting_index = starting_points[idx];
          for(int i(0); i < hole.size() - 1; i++)
          {
            output << starting_index + i << "\n";
          }
          output << starting_index << "\n";
          idx++;
        }
      }

      output << "</DataArray>\n";

      output << "<DataArray type=\"UInt32\" Name=\"offsets\">\n";
      for(int i(1); i < starting_points.size(); i++)
      {
        output << starting_points[i] << "\n";
      }
      output << total_points << "\n";
      output << "</DataArray>\n";

      output << "</Lines>\n";
      output << "<Strips></Strips>\n";
      output << "<Polys></Polys>\n";
      output << "</Piece>\n";
      output << "</PolyData>\n";
      output<< "</VTKFile>\n";
    }
    else {
      std::cerr << "Failed to write polygon to " << filename << "!\n";
    }
  }

  void write_polylines(const std::string& filename, const Polylines& polylines)
  {
    std::ofstream output(filename);

    if (output) {
      int num_lines = polylines.size();
      int total_points = 0;
      std::vector<int> starting_points;

      for(const Polyline& line : polylines)
      {
        starting_points.push_back(total_points);
        total_points += line.size();
      }

      output<< "<VTKFile type=\"PolyData\">\n";
      output << "<PolyData>\n";
      output << "<Piece NumberOfPoints=\"" << total_points << "\" NumberOfVerts=\"0\" NumberOfLines=\"" << num_lines << "\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n";
      output << "<PointData></PointData>\n";
      output << "<CellData></CellData>\n";
      output << "<Points>\n";
      output << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">\n";

      for(const Polyline& polyline : polylines)
      {
        for(const Point& p : polyline)
        {
          output << p.x() << " " << p.y() << " " << p.z() << "\n";
        }
      }

      output << "</DataArray>\n";
      output << "</Points>\n";
      output << "<Verts></Verts>\n";
      output << "<Lines>\n";
      output << "<DataArray type=\"UInt32\" Name=\"connectivity\">\n";

      auto it = polylines.begin();
      for(int i(0); i < polylines.size(); i++)
      {
        int starting_idx = starting_points[i];

        for(int j(0); j < it->size() - 1; j++)
        {
          output << starting_idx + j << "\n";
        }
        output << starting_idx << "\n";
        it++;
      }
      output << "</DataArray>\n";

      output << "<DataArray type=\"UInt32\" Name=\"offsets\">\n";
      for(int i(1); i < starting_points.size(); i++)
      {
        output << starting_points[i] << "\n";
      }
      output << total_points << "\n";
      output << "</DataArray>\n";

      output << "</Lines>\n";
      output << "<Strips></Strips>\n";
      output << "<Polys></Polys>\n";
      output << "</Piece>\n";
      output << "</PolyData>\n";
      output<< "</VTKFile>\n";
    }
    else
    {
      std::cerr << "Failed to write polyline to " << filename << "!\n";
    }
  }

  void write_polylines(const std::string& filename, const Polylines2D& polylines)
  {
    std::ofstream output(filename);

    if (output) {
      int num_lines = polylines.size();
      int total_points = 0;
      std::vector<int> starting_points;

      for(const Polyline2D& line : polylines)
      {
        starting_points.push_back(total_points);
        total_points += line.size();
      }

      output<< "<VTKFile type=\"PolyData\">\n";
      output << "<PolyData>\n";
      output << "<Piece NumberOfPoints=\"" << total_points << "\" NumberOfVerts=\"0\" NumberOfLines=\"" << num_lines << "\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n";
      output << "<PointData></PointData>\n";
      output << "<CellData></CellData>\n";
      output << "<Points>\n";
      output << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">\n";

      for(const Polyline2D& polyline : polylines)
      {
        for(const Point2D& p : polyline)
        {
          output << p.x() << " " << p.y() << " 0\n";
        }
      }

      output << "</DataArray>\n";
      output << "</Points>\n";
      output << "<Verts></Verts>\n";
      output << "<Lines>\n";
      output << "<DataArray type=\"UInt32\" Name=\"connectivity\">\n";

      auto it = polylines.begin();
      for(int i(0); i < polylines.size(); i++)
      {
        int starting_idx = starting_points[i];

        for(int j(0); j < it->size() - 1; j++)
        {
          output << starting_idx + j << "\n";
        }
        output << starting_idx << "\n";
        it++;
      }
      output << "</DataArray>\n";

      output << "<DataArray type=\"UInt32\" Name=\"offsets\">\n";
      for(int i(1); i < starting_points.size(); i++)
      {
        output << starting_points[i] << "\n";
      }
      output << total_points << "\n";
      output << "</DataArray>\n";

      output << "</Lines>\n";
      output << "<Strips></Strips>\n";
      output << "<Polys></Polys>\n";
      output << "</Piece>\n";
      output << "</PolyData>\n";
      output<< "</VTKFile>\n";
    }
    else {
      std::cerr << "Failed to write polyline to " << filename << "!\n";
    }
  }

  Point2D CuttingPlane::project(const Point& point) const
  {
    auto x = CGAL::scalar_product(Vector(origin, point), x_axis);
    auto y = CGAL::scalar_product(Vector(origin, point), y_axis);

    return Point2D(x, y);
  }

  int RadialCrossSectionSampler::num_planes() const
  {
    return _num_planes;
  }

  CuttingPlane RadialCrossSectionSampler::get_plane(int idx) const
  {
    idx = std::max(0, std::min(idx, _num_planes - 1));

    const double delta_angle = (2.0 * M_PI) / double(_num_planes - 1);
    const double angle = idx * delta_angle;

    Vector plane_normal = std::cos(angle) * _u + std::sin(angle) * _v;
    CGALKernel::Plane_3 plane(_origin, plane_normal);

    Vector x_axis = _up;
    Vector y_axis = CGAL::cross_product(plane_normal, x_axis);

    return CuttingPlane {
      plane,
      _origin,
      x_axis,
      y_axis
    };
  }

  int LineCrossSectionSampler::num_planes() const
  {
    return _num_planes;
  }

  CuttingPlane LineCrossSectionSampler::get_plane(int idx) const
  {
    idx = std::max(0, std::min(idx, _num_planes - 1));

    const double delta = 1.0 / double(_num_planes - 1);

    const Point origin = _start + double(idx) * delta * Vector(_start, _end);

    Plane plane(origin, _normal);

    return CuttingPlane
    {
      plane,
      origin,
      _x_axis,
      _y_axis
    };
  }

  /**
   * \brief Computes the union of cross sections of a surface mesh
   *
   * Cross sections are generated by rotating a plane.
   *
   * The axis of rotation is given by the point \c p and the vector \c up.
   * The initial cross section is given by the plane going through \c p with normal \c u.
   * The plane is then rotated 360 degrees around the axis of rotation to produce \c cross_sections cross sections.
   *
   * \returns The union of all generated cross sections.
   *
   * \param filename Path of the surface mesh file
   * \param p Point on cutting plane
   * \param up Up-direction. Cutting plane will run through axis described by p and up
   * \param u Vector orthogonal to up. Initial normal of cutting plane
   * \param cross_sections Number of cross sections
   */
  void union_of_cross_sections(const Mesh& mesh, const CrossSectionSampler& sampler)
  {
    // Slicer constructor from the mesh
    CGAL::Polygon_mesh_slicer<Mesh, CGALKernel> slicer(mesh);

    Polylines polylines;
    Polylines2D projected_polylines;

    PolygonWithHoles shadow;

    for(int i(0); i < sampler.num_planes(); i++)
    {
      // Create slicing plane
      CuttingPlane cutting_plane = sampler.get_plane(i);

      std::cout << "Plane: " << cutting_plane.plane << "\n";

      // Find polylines
      polylines.clear();
      projected_polylines.clear();
      slicer(cutting_plane.plane, std::back_inserter(polylines));

      if(polylines.empty())
      {
        std::cerr << "Empty intersection for plane " << i << "!\n";
        continue;
      }

      /*auto it = polylines.begin();
      for(int j(0); j < polylines.size(); j++)
      {
        write_polyline("intersection_" + std::to_string(i) + "_polyline_" + std::to_string(j) + ".vtp", *it);
        it++;
      }*/
      write_polylines("intersection_" + std::to_string(i) + ".vtp", polylines);

      // The found polylines are 3d objects. We are actually only interested in the projection onto the current cutting plane.
      for(const auto& polyline : polylines)
      {
        Polyline2D projected;
        for(const Point& point : polyline)
        {
          projected.push_back(cutting_plane.project(point));
        }
        projected_polylines.push_back(projected);
      }

      write_polylines("projected_" + std::to_string(i) + ".vtp", projected_polylines);

      PolygonWithHoles next_slice = make_polygon(projected_polylines);

      if(next_slice.outer_boundary().is_counterclockwise_oriented())
      {
        if(i > 0)
        {
          shadow = merge(shadow, next_slice);
        }
        else
        {
          shadow = next_slice;
        }
      }
      else {
        std::cerr << "Outer boundary with wrong orientation. Skipping plane!\n";
      }
    }

    std::cout << "Done!\n";

    std::cout << "Found Union of cross sections with:\n";
    std::cout << "  Boundary Vertices: " << shadow.outer_boundary().size() << "\n";
    std::cout << "  Number of holes: " << shadow.holes().size() << "\n";

    int idx;
    int total_vertices = shadow.outer_boundary().size();
    for(const Polygon& hole : shadow.holes())
    {
      std::cout << "  Hole " << idx << " with " << hole.size() << " vertices\n";
      total_vertices += hole.size();
    }
    std::cout << "  Total vertices: " << total_vertices << "\n";

    // Output the found shadow
    write_polygon("shadow.vtp", shadow);
  }
}
