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
        if(outer->oriented_side(*it) == CGAL::ON_POSITIVE_SIDE)
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

  void write_polyline(const std::string& filename, const Polyline& polyline)
  {
    std::ofstream output(filename);

    if (output) {
      output<< "<VTKFile type=\"PolyData\">\n";
      output << "<PolyData>\n";
      output << "<Piece NumberOfPoints=\"" << polyline.size() << "\" NumberOfVerts=\"0\" NumberOfLines=\"1\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n";
      output << "<PointData>\n";
      //output << "<DataArray type=\"UInt32\" NumberOfComponents=\"1\">\n";

      //for(int i(0); i < polyline.size(); i++)
      //{
        //output << i << " ";
      //}

      //output << "\n</DataArray>\n";
      output << "</PointData>\n";
      output << "<CellData></CellData>\n";
      output << "<Points>\n";
      output << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">\n";

      for(const Point& p : polyline)
      {
        output << p.x() << " " << p.y() << " " << p.z() << "\n";
      }

      output << "</DataArray>\n";
      output << "</Points>\n";
      output << "<Verts></Verts>\n";
      output << "<Lines>\n";
      output << "<DataArray type=\"UInt32\" Name=\"connectivity\">\n";

      for(int i(0); i < polyline.size() - 1; i++)
      {
        output << i << " ";
      }
      output << "0\n";
      output << "</DataArray>\n";

      output << "<DataArray type=\"UInt32\" Name=\"offsets\">\n";
      output << polyline.size() << "\n";
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

  void write_polyline(const std::string& filename, const Polyline2D& polyline)
  {
    std::ofstream output(filename);

    if (output) {
      output<< "<VTKFile type=\"PolyData\">\n";
      output << "<PolyData>\n";
      output << "<Piece NumberOfPoints=\"" << polyline.size() << "\" NumberOfVerts=\"0\" NumberOfLines=\"1\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n";
      output << "<PointData>\n";
      //output << "<DataArray type=\"UInt32\" NumberOfComponents=\"1\">\n";

      //for(int i(0); i < polyline.size(); i++)
      //{
        //output << i << " ";
      //}

      //output << "\n</DataArray>\n";
      output << "</PointData>\n";
      output << "<CellData></CellData>\n";
      output << "<Points>\n";
      output << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">\n";

      for(const Point2D& p : polyline)
      {
        output << "0 " << p.y() << " " << p.x() << "\n";
      }

      output << "</DataArray>\n";
      output << "</Points>\n";
      output << "<Verts></Verts>\n";
      output << "<Lines>\n";
      output << "<DataArray type=\"UInt32\" Name=\"connectivity\">\n";

      for(int i(0); i < polyline.size() - 1; i++)
      {
        output << i << " ";
      }
      output << "0\n";
      output << "</DataArray>\n";

      output << "<DataArray type=\"UInt32\" Name=\"offsets\">\n";
      output << polyline.size() << "\n";
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
  void union_of_cross_sections(const std::string& filename, Point p, Vector up, Vector u, int cross_sections)
  {
    const double eps = 1e-8;

    // Normalize up and u
    // We later need to project onto the various cutting planes and need distances to be compatible
    u = u / CGAL::approximate_sqrt(u.squared_length());
    up = up / CGAL::approximate_sqrt(up.squared_length());

    if(CGAL::scalar_product(u, up) > eps)
    {
      std::cerr << "max_shadow: u and up must be orthogonal!\n";
      return;
    }

    Mesh mesh;
    if(!PMP::IO::read_polygon_mesh(filename, mesh))
    {
      std::cerr << "Could not read mesh!\n";
      return;
    }

    if(CGAL::is_empty(mesh))
    {
      std::cerr << "Input mesh is empty!\n";
      return;
    }

    if(!CGAL::is_triangle_mesh(mesh))
    {
      std::cerr << "Input mesh is not a triangle mesh\n";
      return;
    }

    Vector v = CGAL::cross_product(u, up);

    std::cout << "p: " << p << "\n";
    std::cout << "u: " << u << "\n";
    std::cout << "v: " << v << "\n";
    std::cout << "up: " << up << "\n";

    double angle = 0.0;
    double delta_angle = (2.0 * M_PI) / double(cross_sections - 1);

    // Slicer constructor from the mesh
    CGAL::Polygon_mesh_slicer<Mesh, CGALKernel> slicer(mesh);

    Polylines polylines;
    Polylines2D projected_polylines;

    PolygonWithHoles shadow;

    for(int i(0); i < cross_sections; i++)
    {
      // Create slicing plane
      Vector plane_normal = std::cos(angle) * u + std::sin(angle) * v;
      CGALKernel::Plane_3 plane(p, plane_normal);

      Vector x_axis = up;
      Vector y_axis = CGAL::cross_product(plane_normal, x_axis);

      // Find polylines
      polylines.clear();
      projected_polylines.clear();
      slicer(plane, std::back_inserter(polylines));

      /*auto it = polylines.begin();
      for(int j(0); j < polylines.size(); j++)
      {
        write_polyline("intersection_" + std::to_string(i) + "_polyline_" + std::to_string(j) + ".vtp", *it);
        it++;
      }*/
      write_polyline("intersection_" + std::to_string(i) + ".vtp", polylines.front());

      // The found polylines are 3d objects. We are actually only interested in the projection onto the current cutting plane.
      for(const auto& polyline : polylines)
      {
        Polyline2D projected;
        for(const Point& point : polyline)
        {
          auto x = CGAL::scalar_product(Vector(p, point), x_axis);
          auto y = CGAL::scalar_product(Vector(p, point), y_axis);

          Point2D projected_point(x, y);
          projected.push_back(projected_point);
        }
        projected_polylines.push_back(projected);
      }

      write_polyline("projected_" + std::to_string(i) + ".vtp", projected_polylines.front());

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

      // Update angle
      angle += delta_angle;
    }

    std::cout << "Done!\n";

    // Output the found shadow
    write_polygon("shadow.vtp", shadow.outer_boundary());
  }
}
