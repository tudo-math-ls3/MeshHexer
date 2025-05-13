#include <hexmesher.hpp>

#include <iostream>
#include <math.h>

#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Polyline_simplification_2/simplify.h>

namespace HexMesher
{
  namespace PMP = CGAL::Polygon_mesh_processing;
  namespace PS = CGAL::Polyline_simplification_2;

  /**
   * \brief Returns true if all vertices of \c b are contained in \c a.
   */
  bool contains(const Polygon& a, const Polygon& b)
  {
    for(const Point2D& p : b)
    {
      if(a.bounded_side(p) == CGAL::ON_UNBOUNDED_SIDE)
      {
        return false;
      }
    }
    return true;
  }

  std::vector<PolygonWithHoles> make_polygon(Polylines2D& polylines)
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

    // Simplyify polygons
    int total_vertices_pre = 0;
    int total_vertices_post = 0;
    for(Polygon& poly : polygons)
    {
      total_vertices_pre += poly.size();
      poly = PS::simplify(poly, PS::Squared_distance_cost(), PS::Stop_above_cost_threshold(1e-8));
      total_vertices_post += poly.size();
    }

    std::cout << "Simplified polygons: " << total_vertices_pre << " -> " << total_vertices_post << "\n";

    // Sort polygons in descending order by area, i.e. largest polygon comes first.
    std::sort(
      polygons.begin(),
      polygons.end(),
      [](const Polygon& a, const Polygon& b) { return a.area() > b.area(); }
    );

    // Create inclusion tree for polygons.
    const std::size_t num_polygons = polygons.size();
    std::vector<int> depths(num_polygons, 0);
    std::vector<int> parents(num_polygons, -1);

    std::size_t current_idx = 0;
    for(const Polygon& polygon : polygons)
    {
      int parent_candidate = -1;
      int parent_depth = -1;
      // Find the deepest, i.e. smallest, already handled polygon
      // that still contains the current polygon
      for(int j(0); j < current_idx; j++)
      {
        if(contains(polygons[j], polygon) && depths[j] > parent_depth)
        {
          parent_depth = depths[j];
          parent_candidate = j;
        }
      }

      if(parent_candidate != -1)
      {
        // Polygon is contained in another polygon. Set parent
        parents[current_idx] = parent_candidate;
        depths[current_idx] = parent_depth + 1;
      }
      else
      {
        // Polygon is not contained in another larger polygon. Add it as new root
        parents[current_idx] = current_idx;
        depths[current_idx] = 0;
      }
      current_idx++;
    }

    std::vector<std::size_t> root_mapping(num_polygons, -1);
    std::vector<PolygonWithHoles> result;

    std::size_t handled_polygons = 0;

    // Collect all roots
    for(int i(0); i < num_polygons; i++)
    {
      if(depths[i] == 0)
      {
        root_mapping[i] = result.size();
        result.emplace_back(polygons[i]);
        handled_polygons++;
      }
    }

    // Collect holes at depth 1
    for(int i(0); i < num_polygons; i++)
    {
      if(depths[i] == 1)
      {
        std::size_t parent = parents[i];
        result[root_mapping[parent]].add_hole(polygons[i]);
        handled_polygons++;
      }
    }

    std::size_t unhandled_polygons = num_polygons - handled_polygons;
    if(unhandled_polygons > 0)
    {
      std::cerr << "Warning: " << unhandled_polygons << " unhandled polygons in HexMesher::make_polygon!";
    }

    return result;
  }

  bool merge_polygons_with_holes(const PolygonWithHoles& a, const PolygonWithHoles& b, PolygonWithHoles& res)
  {
    if(CGAL::join(a.outer_boundary(), b.outer_boundary(), res))
    {
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

      res.holes().clear();
      res.holes().insert(res.holes().end(), new_holes.begin(), new_holes.end());
      return true;
    }
    else
    {
      return false;
    }
  }

  std::vector<PolygonWithHoles> merge(std::vector<PolygonWithHoles>& cross_sections)
  {
    std::vector<PolygonWithHoles> result;
    while(!cross_sections.empty())
    {
      PolygonWithHoles merged = cross_sections.back();
      cross_sections.pop_back();

      bool progress = false;
      do
      {
        progress = false;
        for(auto it = cross_sections.rbegin(); it != cross_sections.rend();)
        {
          if(it->outer_boundary().is_clockwise_oriented())
          {
            it = decltype(it){cross_sections.erase(std::next(it).base())};
            continue;
          }

          PolygonWithHoles new_shadow;
          if(merge_polygons_with_holes(merged, *it, new_shadow))
          {
            it = decltype(it){cross_sections.erase(std::next(it).base())};

            merged = new_shadow;
            progress = true;
            continue;
          }
          it++;
        }
      } while(progress);

      merged = PS::simplify(merged, PS::Squared_distance_cost(), PS::Stop_above_cost_threshold(1e-8));
      result.push_back(merged);
    }

    return result;
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

  Point RadialCrossSectionSampler::origin() const
  {
    return _origin;
  }

  int RadialCrossSectionSampler::num_planes() const
  {
    return _num_planes;
  }

  CuttingPlane RadialCrossSectionSampler::get_plane(int idx) const
  {
    idx = std::max(0, std::min(idx, _num_planes - 1));

    const double delta_angle = 
      _num_planes > 1 ? (2.0 * M_PI) / double(_num_planes - 1) : 2.0 * M_PI;
    const double angle = idx * delta_angle;

    Vector plane_normal = std::cos(angle) * _u + std::sin(angle) * _v;
    CGALKernel::Plane_3 plane(_origin, plane_normal);

    Vector y_axis = _up;
    Vector x_axis = CGAL::cross_product(plane_normal, y_axis);

    return CuttingPlane {
      plane,
      _origin,
      x_axis,
      y_axis
    };
  }

  Point2D RadialCrossSectionSampler::project(Point p) const
  {
    // We need a radial projection onto the 0-th cutting plane
    CuttingPlane plane = get_plane(0);

    // Figure out height
    const auto height = CGAL::scalar_product(Vector(plane.origin, p), plane.y_axis);

    Point ref = plane.origin + height * plane.y_axis;

    // p, ref, and the projected point now all lie on a plane orthogonal to the axis
    // Target point is then the point on the 0-th cutting plane with the same distance
    // as the vector p - ref.

    const auto dist = CGAL::approximate_sqrt(Vector(ref, p).squared_length());

    return Point2D(dist, height);
  }

  CuttingPlane RadialCrossSectionSampler::get_plane_through_vertex(Point p) const
  {
    Vector y_axis = _up;
    Vector x_axis = Vector(_origin, p) - y_axis * CGAL::scalar_product(Vector(_origin, p), y_axis);
    x_axis = x_axis / CGAL::approximate_sqrt(x_axis.squared_length());
    Vector normal = CGAL::cross_product(x_axis, y_axis);

    return CuttingPlane {
      CGALKernel::Plane_3(_origin, normal),
      _origin,
      x_axis,
      y_axis
    };
  }

  Point LineCrossSectionSampler::origin() const
  {
    return _start;
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

  Point2D LineCrossSectionSampler::project(Point p) const
  {
    // Orthogonal projection onto 0-th cutting plane

    CuttingPlane plane = get_plane(0);

    const auto x = CGAL::scalar_product(Vector(plane.origin, p), plane.x_axis);
    const auto y = CGAL::scalar_product(Vector(plane.origin, p), plane.y_axis);

    return Point2D(x, y);
  }

  CuttingPlane LineCrossSectionSampler::get_plane_through_vertex(Point p) const
  {
    const auto distance = CGAL::scalar_product(_normal, Vector(_start, p));
    const Point origin = _start + distance * _normal;

    Plane plane(origin, _normal);

    return CuttingPlane
    {
      plane,
      origin,
      _x_axis,
      _y_axis
    };
  }

  template<typename OutputIterator>
  void find_cross_section(CGAL::Polygon_mesh_slicer<Mesh, CGALKernel>& slicer, const CuttingPlane& plane, OutputIterator out)
  {
    Polylines polylines;
    Polylines2D projected_polylines;

    // Find polylines
    polylines.clear();
    projected_polylines.clear();
    slicer(plane.plane, std::back_inserter(polylines));

    if(polylines.empty())
    {
      std::cerr << "Empty intersection.\n";
      return;
    }

    /*auto it = polylines.begin();
    for(int j(0); j < polylines.size(); j++)
    {
      write_polyline("intersection_" + std::to_string(i) + "_polyline_" + std::to_string(j) + ".vtp", *it);
      it++;
    }*/
    //write_polylines("intersection_" + std::to_string(i) + ".vtp", polylines);

    // The found polylines are 3d objects. We are actually only interested in the projection onto the current cutting plane.
    for(const auto& polyline : polylines)
    {
      Polyline2D projected;
      for(const Point& point : polyline)
      {
        projected.push_back(plane.project(point));
      }
      projected_polylines.push_back(projected);
    }

    //write_polylines("projected_" + std::to_string(i) + ".vtp", projected_polylines);

    std::vector<PolygonWithHoles> new_cross_sections = make_polygon(projected_polylines);

    for(PolygonWithHoles& poly : new_cross_sections)
    {
      *out = poly;
    }
  }

  template<typename OutputIterator>
  void find_cross_sections(CGAL::Polygon_mesh_slicer<Mesh, CGALKernel>& slicer, const CrossSectionSampler& sampler, OutputIterator out)
  {
    for(int i(0); i < sampler.num_planes(); i++)
    {
      CuttingPlane cutting_plane = sampler.get_plane(i);
      find_cross_section(slicer, cutting_plane, out);
    }
  }

  template<typename CuttingPlaneIterator, typename OutputIterator>
  void find_cross_sections(CGAL::Polygon_mesh_slicer<Mesh, CGALKernel>& slicer, const CuttingPlaneIterator start, const CuttingPlaneIterator end, OutputIterator out)
  {
    for(CuttingPlaneIterator it = start; it != end; it++)
    {
      find_cross_section(slicer, *it, out);
    }
  }

  bool is_outside_union(Point2D p, std::vector<PolygonWithHoles>& union_components)
  {
    bool outside_boundary = true;
    for(const PolygonWithHoles& poly : union_components)
    {
      outside_boundary = outside_boundary && poly.outer_boundary().bounded_side(p) == CGAL::ON_UNBOUNDED_SIDE;
      for(const Polygon& hole : poly.holes())
      {
        if(hole.bounded_side(p) == CGAL::ON_BOUNDED_SIDE)
        {
          return true;
        }
      }
    }

    return outside_boundary;
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

    std::vector<PolygonWithHoles> all_cross_sections;

    // Running initial sampling of mesh. Collect evenly spaced cross sections and merge them
    find_cross_sections(slicer, sampler, std::back_inserter(all_cross_sections));
    std::vector<PolygonWithHoles> union_components = merge(all_cross_sections);

    // The inital samples should have captured the major features of the input mesh,
    // but there is no guarantee that all features have been captured.
    // As a next step the vertices of the input mesh onto the components of the union.
    // If a vertex gets projected outside of the outer boundary or inside a hole,
    // we have missed a feature of the mesh.
    // In that case we place additional cutting planes through these vertices.

    std::vector<Point> outside_vertices;

    for(const Point& vertex : mesh.points())
    {
      Point2D projected = sampler.project(vertex);
      if(is_outside_union(projected, union_components))
      {
        outside_vertices.push_back(vertex);
      }
    }

    // Merge additional feature planes into union
    std::cout << "Potentially checking an additional " << outside_vertices.size() << " feature planes.\n";

    // Sort vertices by distance to sampler origin.
    // For radial sampler this should save some time,
    // because more distant vertices will cover more other vertices
    // when added to the union.
    std::sort(outside_vertices.begin(), outside_vertices.end(), [&](Point& a, Point& b) {
      return Vector(sampler.origin(), a).squared_length() > Vector(sampler.origin(), b).squared_length();
    });

    std::cout << "Done sorting\n";

    int feature_planes_checked = 0;
    for(Point& vertex : outside_vertices)
    {
      Point2D projected = sampler.project(vertex);
      if(is_outside_union(projected, union_components))
      {
        feature_planes_checked++;
        CuttingPlane plane = sampler.get_plane_through_vertex(vertex);
        find_cross_section(slicer, plane, std::back_inserter(union_components));
        union_components = merge(union_components);
      }
    }

    std::cout << "Checked " << feature_planes_checked << " of " << outside_vertices.size() << " potential feature planes.\n";

    std::cout << "Done!\n";
    std::cout << "Found union of cross sections with " << union_components.size() << " components.\n";

    for(int i(0); i < union_components.size(); i++)
    {
      PolygonWithHoles& component = union_components[i];

      int idx;
      int total_vertices = component.outer_boundary().size();
      for(const Polygon& hole : component.holes())
      {
        std::cout << "  Hole " << idx << " with " << hole.size() << " vertices\n";
        total_vertices += hole.size();
      }
      std::cout << "Component " << i << ":\n";
      std::cout << "  Boundary Vertices: " << component.outer_boundary().size() << "\n";
      std::cout << "  Number of holes: " << component.holes().size() << "\n";
      std::cout << "  Total vertices: " << total_vertices << "\n";

      // Output the found shadow
      write_polygon("shadow_" + std::to_string(i) + ".vtp", component);
    }
  }
}
