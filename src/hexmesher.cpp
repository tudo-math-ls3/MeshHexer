#include <hexmesher.hpp>

#include <iostream>
#include <math.h>
#include <unistd.h>
#include <optional>
#include <variant>

#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Polyline_simplification_2/simplify.h>
#include <CGAL/Polygon_mesh_slicer.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

#include <CGAL/Polygon_mesh_processing/compute_normal.h>

namespace HexMesher
{
  namespace PMP = CGAL::Polygon_mesh_processing;
  namespace PS = CGAL::Polyline_simplification_2;

  void status(const std::string& message)
  {
    if(isatty(fileno(stdout)))
    {
      std::cout << "\33[2K\r"; // Clear line and reset cursor to beginning
      std::cout << message;
      std::cout.flush();
    }
    else
    {
      std::cout << message;
    }
  }

  enum class WindingOrder
  {
    Clockwise,
    CounterClockwise,
  };

  WindingOrder winding_order(const Polyline2D& polyline)
  {
    Real sum(0);
    for(int i(0); i < polyline.size(); i++)
    {
      const Point2D& a = polyline[i];
      const Point2D& b = polyline[(i + 1) % polyline.size()];

      sum += (b.x() - a.x()) * (b.y() + b.y());
    }

    return sum < 0 ? WindingOrder::CounterClockwise : WindingOrder::Clockwise;
  }

  /**
   * \brief Compute left normal of line segment from a to b. Result has unit length.
   */
  Vector2D left_normal(const Point2D& a, const Point2D& b)
  {
    Vector2D delta = b - a;
    Vector2D normal = Vector2D(-delta.y(), delta.x());
    return normal / CGAL::approximate_sqrt(normal.squared_length());
  }

  /**
   * \brief Compute right normal of line segment from a to b. Result has unit length.
   */
  Vector2D right_normal(const Point2D& a, const Point2D& b)
  {
    return -left_normal(a, b);
  }

  Vector2D outside_normal(const Point2D& a, const Point2D& b, WindingOrder order)
  {
    return order == WindingOrder::CounterClockwise ? right_normal(a, b) : left_normal(a, b);
  }

  Vector2D inside_normal(const Point2D& a, const Point2D& b, WindingOrder order)
  {
    return -outside_normal(a, b, order);
  }

  double angle(const Vector2D& a, const Vector2D& b)
  {
    return CGAL::to_double(CGAL::approximate_angle(Vector(a.x(), a.y(), 0), Vector(b.x(), b.y(), 0)));
  }

  /**
   * \brief Simplifies a polyline by merging runs of segments with similar normals
   * 
   * Normals whose angles differ by at most \c threshold are considered similar.
   * 
   * \param polyline The polyline to simplify
   * \param threshold Similarity threshold for normals, in degrees.
   */
  std::pair<Polygon, std::vector<int>> simplify_by_normal(
    const Polygon& polygon, 
    const std::function<bool(const std::vector<Vector2D>&, const Vector2D&)>& continue_pred)
  {
    Polyline2D result;
    std::vector<int> chosen_points;

    // Find vertex of polyline with smallest inside angle as a safe starting point
    const std::size_t size = polygon.size();
    int starting_vertex = 0;
    Real minimal_angle = 360.0;
    for(int i(0); i < polygon.size(); i++)
    {
      const Point2D& a = polygon[i];
      const Point2D& b = polygon[(i + 1) % size];
      const Point2D& c = polygon[(i + 2) % size];

      Vector2D ba = a - b;
      ba = ba / CGAL::approximate_sqrt(ba.squared_length());
      Vector2D bc = c - b;
      bc = bc / CGAL::approximate_sqrt(bc.squared_length());

      const Real theta = angle(ba, bc);

      if(theta < minimal_angle)
      {
        starting_vertex = i + 1;
        minimal_angle = theta;
      }
    }

    // Build simplified result by finding sections with similar normals
    WindingOrder order = polygon.is_clockwise_oriented() ? WindingOrder::Clockwise : WindingOrder::CounterClockwise;

    int section_start = starting_vertex;
    int section_end = (starting_vertex + 1) % size;
    int candidate = (starting_vertex + 2) % size;

    Vector2D reference_normal = outside_normal(polygon[section_start], polygon[section_end], order);
    Vector2D link_normal = outside_normal(polygon[section_end], polygon[candidate], order);

    result.push_back(polygon[section_start]);
    chosen_points.push_back(section_start);

    std::vector<Vector2D> section_normals;
    section_normals.push_back(reference_normal);
    
    while(true)
    {
      while(continue_pred(section_normals, link_normal) && section_end != starting_vertex)
      {
        section_end = candidate;
        candidate = (candidate + 1) % size;

        const Point2D& a = polygon[section_end];
        const Point2D& b = polygon[candidate];

        section_normals.push_back(link_normal);
        link_normal = outside_normal(a, b, order);
      }

      // Latest link_normal is outside threshold
      // Merge the links from [section_start, section_end]
      // and continue next section at section_end if not all links are covered yet.

      if(section_end == starting_vertex)
      {
        // With this section we have covered all links of the original polyline.
        // The last vertex we pushed thus just gets connected to the inital vertex.

        std::cout << "Polygon simplification done. Reduced from " << polygon.size() << " vertices to " << result.size() << " vertices.\n";
        return {Polygon(result.begin(), result.end()), chosen_points};
      }

      // We have not yet covered all links of the original polyline.
      // Merge the current section into a single segment and then continue.
      result.push_back(polygon[section_end]);
      chosen_points.push_back(section_end);

      section_start = section_end;
      section_end = candidate;
      candidate = (candidate + 1) % size;

      reference_normal = outside_normal(polygon[section_start], polygon[section_end], order);
      link_normal = outside_normal(polygon[section_end], polygon[candidate], order);

      section_normals.clear();
      section_normals.push_back(reference_normal);

      status("[SimplifyByNormal]: " + std::to_string(section_end) + " / " + std::to_string(starting_vertex));
    }

    return {Polygon(), std::vector<int>()};
  }

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
      poly = PS::simplify(poly, PS::Squared_distance_cost(), PS::Stop_above_cost_threshold(1e-4));
      total_vertices_post += poly.size();
    }

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
      status("Finding cross section " + std::to_string(i) + " of " + std::to_string(sampler.num_planes()));
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
  std::vector<PolygonWithHoles> union_of_cross_sections(const Mesh& mesh, const CrossSectionSampler& sampler)
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

    /*
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
        CuttingPlane plane = sampler.get_plane_through_vertex(vertex);
        find_cross_section(slicer, plane, std::back_inserter(union_components));
        union_components = merge(union_components);
      }
      feature_planes_checked++;
      status("Checked " + std::to_string(feature_planes_checked + 1) + " of " + std::to_string(outside_vertices.size()) + " feature planes");
    }

    std::cout << "Checked " << feature_planes_checked << " of " << outside_vertices.size() << " potential feature planes.\n";
    */


    std::cout << "Done!\n";
    std::cout << "Found union of cross sections with " << union_components.size() << " components.\n";

    return union_components;
  }

  std::vector<Vector2D> laplace(const Polygon& polygon)
  {
    std::vector<Vector2D> result(polygon.size());

    const std::size_t size = polygon.size();
    for(int i(0); i < polygon.size(); i++)
    {
      const Point2D& a = polygon[i];
      const Point2D& b = polygon[(i + 1) % size];
      const Point2D& c = polygon[(i + 2) % size];

      const Real length_ab = CGAL::approximate_sqrt((b - a).squared_length());
      const Real length_bc = CGAL::approximate_sqrt((c - b).squared_length());
      const Real total_length = length_ab + length_bc;

      const Vector2D delta = (length_bc * Vector2D(b, a) + length_ab * Vector2D(b, c)) / total_length;

      result[(i + 1) % size] = delta;
    }

    return result;
  }

  struct Ray
  {
    Point2D origin;
    Vector2D direction;
  };

  struct Segment
  {

    Segment() = default;

    Segment(const Point2D& a_, const Point2D& b_, const Vector2D& normal_) :
      a(a_),
      b(b_),
      normal(normal_)
    {
    }

    Point2D a;
    Point2D b;

    Vector2D normal;
  };

  struct Intersection
  {
    Real distance;
    Point2D point;
    Vector2D normal;
    bool outside;

    int index;
  };

  /**
   * Returns the intersection between the given ray and the given segment, if they intersect.
   */
  std::optional<Intersection> ray_segment_intersection(const Ray& ray, const Segment& segment)
  {
    // We need to solve
    // ray.origin + t * ray.direction = segment.a + u * (segment.b - segment.a)
    // or in short
    // o + td = a + u(b - a)
    // for t, u.
    // If t >= 0 and 0 <= u <= 1, the lines intersect.

    // Rearrange to
    // td - u(b - a) = a - o
    // and use s = a - b, c = a - o
    // td + us = c

    // Segment direction
    const Vector2D s = segment.a - segment.b;
    // Rhs of system of equations
    const Vector2D c = segment.a - ray.origin;
    const Vector2D d = ray.direction;

    // Solution via Cramers rule
    const Real denom = d.x() * s.y() - s.x() * d.y();

    if(denom == 0)
    {
      // Line and segment are either colinear or parallel
      // There is technically a intersection in the colinear case
      // but we don't handle that right know
      return std::nullopt;
    }

    const Real t = (c.x() * s.y() - s.x() * c.y()) / denom;
    const Real u = (d.x() * c.y() - c.x() * d.y()) / denom;

    if(t > 0 && 0 < u && u < 1)
    {
      // Intersection between line and segment

      const Point2D point = ray.origin + t * ray.direction;
      const Real distance = CGAL::approximate_sqrt((t * ray.direction).squared_length());

      Vector2D normal = segment.normal;
      bool outside = true;
      if(CGAL::scalar_product(normal, ray.direction) > 0)
      {
        normal = -normal;
        outside = false;
      }

      std::cout << "Intersection:\n";
      std::cout << "Ray { origin: " << ray.origin << ", direction: " << ray.direction << "}\n";
      std::cout << "Segment {a: " << segment.a << ", " << segment.b << " }\n";
      std::cout << "Intersection point: " << point << "\n";
      std::cout << "t: " << t << "\n";
      std::cout << "u: " << u << "\n";

      return Intersection{distance, point, normal, outside};
    }

    return std::nullopt;
  };

  template<typename SegmentIter>
  std::optional<Intersection> closest_ray_segment_intersection(const Ray& ray, SegmentIter begin, SegmentIter end)
  {
    std::optional<Intersection> result(std::nullopt);

    int idx = 0;
    for(SegmentIter it = begin; it != end; it++)
    {
      std::optional<Intersection> intersection = ray_segment_intersection(ray, *it);

      if(!result && intersection)
      {
        result = intersection;
        result.value().index = idx;
      }

      if(result && intersection && result.value().distance > intersection.value().distance)
      {
        result = intersection;
        result.value().index = idx;
      }

      idx++;
    }

    return result;
  }

  Polygon grid_sample(const Polygon& polygon, Real min_dist)
  {
    // 1. Determine characteristic points of the polygon
    // Characteristic points are those points we definitely want to keep in the
    // final mesh, because they belong to geometric features we want to keep

    // Choosing characteristic points is difficult.
    // As a proof of concept we choose
    // - points with interior angles <= 90deg
    // - points with interior angles >= 270deg, i.e. inner and outer 90deg corners
    // - points chosen by a simplify_by_normals with an allowed normal difference of 20deg.
    // - points with large curvature, top 5%

    // We have to be careful to balance capturing all geometric features
    // with the number of chosen points.
    // Gmsh does (per default) not discard boundary points,
    // which means choosing too many points causes too small elements.

    std::set<int> characteristic_points;

    // Choose points by interior angle
    const std::size_t size = polygon.size();
    for(int i(0); i < polygon.size(); i++)
    {
      const Point2D& a = polygon[i];
      const Point2D& b = polygon[(i + 1) % size];
      const Point2D& c = polygon[(i + 2) % size];

      Vector2D ba = a - b;
      ba = ba / CGAL::approximate_sqrt(ba.squared_length());
      Vector2D bc = c - b;
      bc = bc / CGAL::approximate_sqrt(bc.squared_length());

      const Real theta = angle(ba, bc);

      if(theta <= 90)
      {
        characteristic_points.insert(i + 1);
      }

      if(theta >= 270)
      {
        characteristic_points.insert(i + 1);
      }
    }

    // Choose points of simplified polygon
    std::vector<int> simplified_points = simplify_by_normal(polygon, [](auto& normals, auto& next)
    {
      return angle(normals.front(), next) < 20.0;
    }).second;;

    characteristic_points.insert(simplified_points.begin(), simplified_points.end());

    // Choose top points based on curvature

    // Determine laplace for every vertex
    std::vector<Vector2D> offsets = laplace(polygon);

    // Sort indices based on laplace, in descending order
    std::vector<int> indices(polygon.size());
    for(int i(0); i < polygon.size(); i++)
    {
      indices[i] = i;
    }
    std::sort(indices.begin(), indices.end(), [&](int a, int b) {
      return CGAL::approximate_sqrt(offsets[a].squared_length()) > CGAL::approximate_sqrt(offsets[b].squared_length());
    });

    // Pick top 5% of vertices by curvature

    for(int i(0); i < (double)polygon.size() * 0.05; i++)
    {
      characteristic_points.insert(indices[i]);
    }

    // 2. Create a grid from the chosen characteristic points
    // Each characteristic point describes an intersection of grid lines.
    // The grid only extends into the interior of the polygon,
    // i.e. a inner corner will stop the grid.
    // This shields parts of the polygon from vertices, which shouldn't
    // influence the meshing somewhere else.
    // Imagine a lopsided y-pipe. Characteristic pipes in one of the
    // pipes should not influence the meshing of the other pipe.
    // Resample the polyline by choosing all characteristic points
    // and collecting all segment-grid intersections.
    // The resulting sampling lies on a tensor-grid,
    // which should hopefully allow the mesher to create nicer meshes
    // than when sampling the polyline at random.

    // Build polyline from points chosen so far
    std::vector<Point2D> polyline;
    for(int idx : characteristic_points)
    {
      std::cout << "Adding vertex " << idx << " to new polyline\n";
      polyline.push_back(polygon[idx]);
    }

    // Build list of segments
    const WindingOrder order = polygon.is_clockwise_oriented() ? WindingOrder::Clockwise : WindingOrder::CounterClockwise;
    std::vector<Segment> segments;
    for(int i(0); i < polyline.size(); i++)
    {
      const Point2D a = polyline[i];
      const Point2D b = polyline[(i + 1) % polyline.size()];
      segments.push_back(Segment(a, b, outside_normal(a, b, order)));
    }

    std::array<Vector2D, 4> directions = {Vector2D(-1.0, 0.0), Vector2D(1.0, 0.0), Vector2D(0.0, -1.0), Vector2D(0.0, 1.0)};
    for(const Point2D& point : polyline)
    {
      for(const Vector2D& dir : directions)
      {
        Ray ray{point, dir};
        std::optional<Intersection> intersection = closest_ray_segment_intersection(ray, segments.begin(), segments.end());

        if(intersection && !intersection.value().outside)
        {
          const int idx = intersection.value().index;
          const Point2D intersection_point = intersection.value().point;

          // Update segment list
          // Erase the hit segment and insert the two new segments at its place
          auto segment_iter = segments.begin() + idx;
          Segment hit = *segment_iter;

          segment_iter = segments.erase(segment_iter);
          segment_iter = segments.insert(segment_iter, Segment(intersection_point, hit.b, outside_normal(intersection_point, hit.b, order)));
          segments.insert(segment_iter, Segment(hit.a, intersection_point, outside_normal(hit.a, intersection_point, order)));
        }
      }
    }

    {
      for(int i(0); i < segments.size() - 1; i++)
      {
        const Segment& a = segments[i];
        const Segment& b = segments[i + 1];

        if(a.b != b.a)
        {
          std::cerr << "Segments not contiguous!\n";
          exit(1);
        }
      }
      const Segment& start = segments[0];
      const Segment& end = segments[segments.size() - 1];

      if(start.a != end.b)
      {
        std::cerr << "Segments not contiguous!\n";
        exit(1);
      }

    }

    // Filter out too small segments
    for(auto it = segments.begin() + 1; it <= segments.end() - 1;)
    {
      const Segment& segment = *it;
      const auto length = CGAL::approximate_sqrt((segment.b - segment.a).squared_length());

      const Point2D midpoint = segment.a + (length / 2.0) * (segment.b - segment.a);

      if(length < min_dist)
      {
        Segment& prev = *std::prev(it);
        Segment& next = *std::next(it);

        prev.b = midpoint;
        next.a = midpoint;

        it = segments.erase(it);
      }
      else
      {
        it++;
      }
    }

    std::vector<Point2D> grid_sampled_polyline;

    for(const Segment& segment : segments)
    {
      grid_sampled_polyline.push_back(segment.a);
    }

    std::cout << "Finished grid sampling. Reduced from " << polygon.size() << " vertices to " << grid_sampled_polyline.size() << " vertices.\n";
    return Polygon(grid_sampled_polyline.begin(), grid_sampled_polyline.end());
  }

  void compute_mesh_thickness(Mesh& mesh)
  {
    if(mesh.num_faces() == 0)
    {
      // No work to be done
      return;
    }

    // Create AABB tree for intersection and distance queries
    using Primitive = CGAL::AABB_face_graph_triangle_primitive<Mesh>;
    using AABBTraits = CGAL::AABB_traits_3<CGALKernel, Primitive>;
    using AABBTree = CGAL::AABB_tree<AABBTraits>;

    using FaceIndex = Mesh::Face_index;

    using Ray3 = CGALKernel::Ray_3;
    using Segment3 = CGALKernel::Segment_3;

    using RayIntersection = std::optional<AABBTree::Intersection_and_primitive_id<Ray3>::Type>;

    AABBTree aabb_tree(mesh.faces_begin(), mesh.faces_end(), mesh);

    // Create mesh property for storing thicknesses
    Mesh::Property_map<FaceIndex, double> diameter_property = 
      mesh.add_property_map<FaceIndex, double>("f:MIS_diameter", 0).first;

    Mesh::Property_map<FaceIndex, FaceIndex> id_property = 
      mesh.add_property_map<FaceIndex, FaceIndex>("f:MIS_id", FaceIndex(0)).first;

    // Determine bounding box of mesh
    Real max_radius = std::min(
      {
        aabb_tree.bbox().xmax() - aabb_tree.bbox().xmin(),
        aabb_tree.bbox().ymax() - aabb_tree.bbox().ymin(),
        aabb_tree.bbox().zmax() - aabb_tree.bbox().zmin(),
      }
    ) / Real(2.0);

    for(FaceIndex face_index : mesh.faces())
    {
      // Determine centroid of face
      Point centroid(0.0, 0.0, 0.0);
      auto circulators = CGAL::vertices_around_face(mesh.halfedge(face_index), mesh);
      for(auto vcirc = circulators.first; vcirc != circulators.second; vcirc++)
      {
        centroid += Real(1.0 / 3.0) * Vector(Point(CGAL::Origin()), mesh.point(*vcirc));
      }

      // Determine initial sphere radius 
      // We need an initial radius that is large enough for our initial sphere
      // to intersect at least one object other than the current face.
      // We cast a ray from the centroid towards the face normal.
      // If we hit something, we have an initial guess for the radius.
      // Otherwise the current face is pointing towards a hole in the surface mesh.
      // In that case we use the smallest length of the axis aligned bounding box as the initial diameter.
      // This is the largest value any MIS can be.

      Real initial_radius(max_radius);
      FaceIndex initial_id(0);

      // Determine (inward) normal
      Vector normal = -CGAL::Polygon_mesh_processing::compute_face_normal(face_index, mesh);

      // Cast ray
      Ray3 ray(centroid, normal);
      auto skip = [=](FaceIndex idx) { return idx == face_index; };
      RayIntersection intersection = aabb_tree.first_intersection(ray, skip);

      if(intersection && std::holds_alternative<Point>(intersection.value().first))
      {
        const Point& p = std::get<Point>(intersection.value().first);
        Real distance = CGAL::approximate_sqrt(Vector(centroid, p).squared_length());
        initial_radius = std::min(initial_radius, distance / 2.0);
        initial_id = intersection.value().second;
      }

      if(intersection && std::holds_alternative<Segment3>(intersection.value().first))
      {
        const Segment3& segment = std::get<Segment3>(intersection.value().first);
        Real distance = CGAL::approximate_sqrt(Vector(centroid, segment.source()).squared_length());
        initial_radius = std::min(initial_radius, distance / 2.0);
        initial_id = intersection.value().second;
      }

      Real prev_radius = initial_radius + 1;
      Real radius = initial_radius;
      FaceIndex id = initial_id;

      // Shrink sphere until change in radius becomes too small
      while(prev_radius - radius > 1e-6)
      {
        Point center = centroid + radius * normal;
        auto closest_point_data = aabb_tree.closest_point_and_primitive(center);

        Point closest_point = closest_point_data.first;
        Real distance_to_closest = CGAL::approximate_sqrt(Vector(center, closest_point).squared_length());

        if(distance_to_closest >= radius - 1e-6)
        {
          // Closest point is not in sphere. Current sphere is correct
          break;
        }

        // We now need to find the radius r' of the new sphere that touches both the
        // centroid of our face and the closest point we just found.
        // Let the centroid be p, the closest point p' and let c' be the center of the new sphere.
        // p, p', and c' prime form a isosceles triangle with two sides of length r', a side of length d.
        //
        //     r'
        // p'------c'
        // \       |
        //  \      |
        //   \     |
        //  d \    | r'
        //     \   |
        //      \  |
        //       \a|
        //        p
        //
        // The angle a is the angle between the normal of the face and the vector pp'.
        // We can then determine r' via the law of sines as r' = d * sin(a) * sin(180 - 2a)

        double alpha = CGAL::to_double(CGAL::approximate_angle(normal, Vector(centroid, closest_point))) / 180.0 * M_PI;

        if(alpha > 1.3)
        {
          // Closest point is at a very shallow angle.
          // Likely a meshing artifact. Ignore it.
          break;
        }

        double d = CGAL::to_double(CGAL::approximate_sqrt((closest_point - centroid).squared_length()));
        Real new_radius = Real(d * std::sin(alpha) / std::sin(M_PI - 2.0 * alpha));

        prev_radius = radius;
        radius = std::min(new_radius, radius);
        id = closest_point_data.second;
      }

      // Record thickness
      diameter_property[face_index] = CGAL::to_double(2.0 * radius);
      id_property[face_index] = id;
    }
  }
}
