#include <CGAL/enum.h>
#include <meshing.hpp>
#include <properties.hpp>

#include <limits>

#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Polygon_mesh_slicer.h>
#include <CGAL/Polyline_simplification_2/simplify.h>

namespace HexMesher
{
  namespace PS = CGAL::Polyline_simplification_2;

  namespace
  {
    void status(const std::string& message)
    {
      if(isatty(fileno(stdout)) != 0)
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

    enum class WindingOrder : std::uint8_t
    {
      Clockwise,
      CounterClockwise,
    };

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

    double angle(const Vector2D& a, const Vector2D& b)
    {
      return CGAL::to_double(CGAL::approximate_angle(Vector3D(a.x(), a.y(), 0), Vector3D(b.x(), b.y(), 0)));
    }

    /**
    * \brief Returns true if all vertices of \c b are contained in \c a.
    */
    bool contains(const Polygon2D& a, const Polygon2D& b)
    {
      return std::all_of(b.begin(), b.end(), [&a](const Point2D& p) {
        return a.bounded_side(p) != CGAL::ON_UNBOUNDED_SIDE;
      });
    }

    std::vector<PolygonWithHoles2D> make_polygon(Polylines2D& polylines)
    {
      // Create polygons from the polylines.
      // This gives us inside-out tests without
      // without having to think about the orientation
      // of the polyline
      std::vector<Polygon2D> polygons;

      for(auto & polyline : polylines)
      {
        Polygon2D poly(polyline.begin(), polyline.end() - 1);
        if(poly.is_clockwise_oriented())
        {
          polygons.emplace_back(std::reverse_iterator(polyline.end() - 1), std::reverse_iterator(polyline.begin()));
        }
        else
        {
          polygons.emplace_back(std::move(poly));
        }
      }

      // Simplyify polygons
      for(Polygon2D& poly : polygons)
      {
        poly = PS::simplify(poly, PS::Squared_distance_cost(), PS::Stop_above_cost_threshold(1e-4));
      }

      // Sort polygons in descending order by area, i.e. largest polygon comes first.
      std::sort(
        polygons.begin(),
        polygons.end(),
        [](const Polygon2D& a, const Polygon2D& b) { return a.area() > b.area(); });

      // Create inclusion tree for polygons.
      const std::size_t num_polygons = polygons.size();
      std::vector<int> depths(num_polygons, 0);
      std::vector<std::size_t> parents(num_polygons, -1);

      std::size_t current_idx = 0;
      for(const Polygon2D& polygon : polygons)
      {
        std::size_t parent_candidate = ~std::size_t(0);
        int parent_depth = -1;
        // Find the deepest, i.e. smallest, already handled polygon
        // that still contains the current polygon
        for(std::size_t j(0); j < current_idx; j++)
        {
          if(contains(polygons[j], polygon) && depths[j] > parent_depth)
          {
            parent_depth = depths[j];
            parent_candidate = j;
          }
        }

        if(parent_candidate != ~std::size_t(0))
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
      std::vector<PolygonWithHoles2D> result;

      std::size_t handled_polygons = 0;

      // Collect all roots
      for(std::size_t i(0); i < num_polygons; i++)
      {
        if(depths[i] == 0)
        {
          root_mapping[i] = result.size();
          result.emplace_back(polygons[i]);
          handled_polygons++;
        }
      }

      // Collect holes at depth 1
      for(std::size_t i(0); i < num_polygons; i++)
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

    bool merge_polygons_with_holes(const PolygonWithHoles2D& a, const PolygonWithHoles2D& b, PolygonWithHoles2D& res)
    {
      if(CGAL::join(a.outer_boundary(), b.outer_boundary(), res))
      {
        std::vector<PolygonWithHoles2D> intersections;

        for(const auto& hole_a : a.holes())
        {
          for(const auto& hole_b : b.holes())
          {
            CGAL::intersection(hole_a, hole_b, std::back_inserter(intersections));
          }
        }

        std::vector<Polygon2D> new_holes;

        new_holes.reserve(intersections.size());
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

    std::vector<PolygonWithHoles2D> merge(std::vector<PolygonWithHoles2D>& cross_sections)
    {
      std::vector<PolygonWithHoles2D> result;
      while(!cross_sections.empty())
      {
        PolygonWithHoles2D merged = cross_sections.back();
        cross_sections.pop_back();

        bool progress = false;
        do //NOLINT
        {
          progress = false;
          for(auto it = cross_sections.rbegin(); it != cross_sections.rend();)
          {
            if(it->outer_boundary().is_clockwise_oriented())
            {
              it = decltype(it){cross_sections.erase(std::next(it).base())};
              continue;
            }

            PolygonWithHoles2D new_shadow;
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

    template<typename OutputIterator>
    void
    find_cross_section(CGAL::Polygon_mesh_slicer<Mesh, Kernel>& slicer, const CuttingPlane& plane, OutputIterator out)
    {
      Polylines3D polylines;
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

      // The found polylines are 3d objects. We are actually only interested in the projection onto the current cutting
      // plane.
      for(const auto& polyline : polylines)
      {
        Polyline2D projected;
        for(const Point3D& point : polyline)
        {
          projected.push_back(plane.project(point));
        }
        projected_polylines.push_back(projected);
      }

      std::vector<PolygonWithHoles2D> new_cross_sections = make_polygon(projected_polylines);

      for(PolygonWithHoles2D& poly : new_cross_sections)
      {
        *out = poly;
      }
    }

    template<typename OutputIterator>
    void find_cross_sections(
      CGAL::Polygon_mesh_slicer<Mesh, Kernel>& slicer,
      const CrossSectionSampler& sampler,
      OutputIterator out)
    {
      for(int i(0); i < sampler.num_planes(); i++)
      {
        status("Finding cross section " + std::to_string(i) + " of " + std::to_string(sampler.num_planes()));
        CuttingPlane cutting_plane = sampler.get_plane(i);
        find_cross_section(slicer, cutting_plane, out);
      }
    }

    template<typename CuttingPlaneIterator, typename OutputIterator>
    void find_cross_sections(
      CGAL::Polygon_mesh_slicer<Mesh, Kernel>& slicer,
      const CuttingPlaneIterator start,
      const CuttingPlaneIterator end,
      OutputIterator out)
    {
      for(CuttingPlaneIterator it = start; it != end; it++)
      {
        find_cross_section(slicer, *it, out);
      }
    }

    template<typename PointType, typename PointIterator>
    std::vector<PointType> k_means_clustering(PointIterator begin, PointIterator end, std::size_t num_clusters)
    {
      if(std::distance(begin, end) < num_clusters)
      {
        num_clusters = std::distance(begin, end);
      }

      std::vector<PointType> cluster_centers(num_clusters);
      std::vector<std::vector<PointIterator>> clusters(num_clusters);

      const auto distance = [](const PointType& a, const PointType& b) { return (a - b).squared_length(); };

      PointIterator init = begin;
      for(std::size_t i = 0; i < num_clusters; i++)
      {
        cluster_centers[i] = *init++;
      }

      bool converged = false;
      while(!converged)
      {
        // Reset clusters
        for(auto& cluster : clusters)
        {
          cluster.clear();
        }

        // Assign each gap to closest cluster
        PointIterator next_point = begin;
        while(next_point != end)
        {
          std::size_t closest = 0;
          double closest_distance = std::numeric_limits<double>::max();

          for(std::size_t j(0); j < num_clusters; j++)
          {
            double d = distance(*next_point, cluster_centers[j]);
            if(d < closest_distance)
            {
              closest = j;
              closest_distance = d;
            }
          }

          clusters[closest].push_back(next_point);
          next_point++;
        }

        // Determine new cluster centers
        std::vector<PointType> new_centers;
        for(std::size_t cluster(0); cluster < num_clusters; cluster++)
        {
          if(clusters[cluster].size() > 0)
          {
            PointType center{CGAL::NULL_VECTOR};
            for(const auto& datum : clusters[cluster])
            {
              center += *datum;
            }
            center /= Real(clusters[cluster].size());

            new_centers.push_back(center);
          }
          else
          {
            new_centers.push_back(*begin);
          }
        }

        // Check for convergence
        converged = true;
        for(std::size_t cluster(0); cluster < num_clusters; cluster++)
        {
          if(distance(cluster_centers[cluster], new_centers[cluster]) > 0.01)
          {
            converged = false;
          }
        }

        cluster_centers = new_centers;
      }

      return cluster_centers;
    }

    enum class Axis : std::uint8_t
    {
      X = 0,
      Y = 1,
      Z = 2,
    };

    std::array<Axis, 2> other_axes(Axis axis)
    {
      switch(axis)
      {
        case Axis::X: return {Axis::Y, Axis::Z};
        case Axis::Y: return {Axis::X, Axis::Z};
        default: return {Axis::X, Axis::Y};
      }
    }
  }

  Point2D CuttingPlane::project(const Point3D& point) const
  {
    auto x = CGAL::scalar_product(Vector3D(origin, point), x_axis);
    auto y = CGAL::scalar_product(Vector3D(origin, point), y_axis);

    return {x, y};
  }

  Point3D RadialCrossSectionSampler::origin() const
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

    const double delta_angle = _num_planes > 1 ? (2.0 * M_PI) / double(_num_planes - 1) : 2.0 * M_PI;
    const double angle = idx * delta_angle;

    Vector3D plane_normal = std::cos(angle) * _u + std::sin(angle) * _v;
    Plane3D plane(_origin, plane_normal);

    Vector3D y_axis = _up;
    Vector3D x_axis = CGAL::cross_product(plane_normal, y_axis);

    return CuttingPlane{plane, _origin, x_axis, y_axis};
  }

  Point2D RadialCrossSectionSampler::project(Point3D p) const
  {
    // We need a radial projection onto the 0-th cutting plane
    CuttingPlane plane = get_plane(0);

    // Figure out height
    const auto height = CGAL::scalar_product(Vector3D(plane.origin, p), plane.y_axis);

    Point3D ref = plane.origin + height * plane.y_axis;

    // p, ref, and the projected point now all lie on a plane orthogonal to the axis
    // Target point is then the point on the 0-th cutting plane with the same distance
    // as the vector p - ref.

    const auto dist = CGAL::approximate_sqrt(Vector3D(ref, p).squared_length());

    return {dist, height};
  }

  CuttingPlane RadialCrossSectionSampler::get_plane_through_vertex(Point3D p) const
  {
    Vector3D y_axis = _up;
    Vector3D x_axis = Vector3D(_origin, p) - y_axis * CGAL::scalar_product(Vector3D(_origin, p), y_axis);
    x_axis = x_axis / CGAL::approximate_sqrt(x_axis.squared_length());
    Vector3D normal = CGAL::cross_product(x_axis, y_axis);

    return CuttingPlane{Plane3D(_origin, normal), _origin, x_axis, y_axis};
  }

  Point3D LineCrossSectionSampler::origin() const
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

    const Point3D origin = _start + double(idx) * delta * Vector3D(_start, _end);

    Plane3D plane(origin, _normal);

    return CuttingPlane{plane, origin, _x_axis, _y_axis};
  }

  Point2D LineCrossSectionSampler::project(Point3D p) const
  {
    // Orthogonal projection onto 0-th cutting plane

    CuttingPlane plane = get_plane(0);

    const auto x = CGAL::scalar_product(Vector3D(plane.origin, p), plane.x_axis);
    const auto y = CGAL::scalar_product(Vector3D(plane.origin, p), plane.y_axis);

    return {x, y};
  }

  CuttingPlane LineCrossSectionSampler::get_plane_through_vertex(Point3D p) const
  {
    const auto distance = CGAL::scalar_product(_normal, Vector3D(_start, p));
    const Point3D origin = _start + distance * _normal;

    Plane3D plane(origin, _normal);

    return CuttingPlane{plane, origin, _x_axis, _y_axis};
  }

  /**
   * \brief Simplifies a polyline by merging runs of segments with similar normals
   *
   * Normals whose angles differ by at most \c threshold are considered similar.
   *
   * \param polyline The polyline to simplify
   * \param threshold Similarity threshold for normals, in degrees.
   */
  std::pair<Polygon2D, std::vector<std::size_t>> simplify_by_normal(
    const Polygon2D& polygon,
    const std::function<bool(const std::vector<Vector2D>&, const Vector2D&)>& continue_pred)
  {
    Polyline2D result;
    std::vector<std::size_t> chosen_points;

    // Find vertex of polyline with smallest inside angle as a safe starting point
    const std::size_t size = polygon.size();
    std::size_t starting_vertex = 0;
    Real minimal_angle = 360.0;
    for(std::size_t i(0); i < polygon.size(); i++)
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

    std::size_t section_start = starting_vertex;
    std::size_t section_end = (starting_vertex + 1) % size;
    std::size_t candidate = (starting_vertex + 2) % size;

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

        std::cout << "Polygon simplification done. Reduced from " << polygon.size() << " vertices to " << result.size()
                  << " vertices.\n";
        return {Polygon2D(result.begin(), result.end()), chosen_points};
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

    return {Polygon2D(), std::vector<std::size_t>()};
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
  std::vector<PolygonWithHoles2D> union_of_cross_sections(const Mesh& mesh, const CrossSectionSampler& sampler)
  {
    // Slicer constructor from the mesh
    CGAL::Polygon_mesh_slicer<Mesh, Kernel> slicer(mesh);

    std::vector<PolygonWithHoles2D> all_cross_sections;

    // Running initial sampling of mesh. Collect evenly spaced cross sections and merge them
    find_cross_sections(slicer, sampler, std::back_inserter(all_cross_sections));
    std::vector<PolygonWithHoles2D> union_components = merge(all_cross_sections);

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
      status("Checked " + std::to_string(feature_planes_checked + 1) + " of " + std::to_string(outside_vertices.size())
    + " feature planes");
    }

    std::cout << "Checked " << feature_planes_checked << " of " << outside_vertices.size() << " potential feature
    planes.\n";
    */

    std::cout << "Done!\n";
    std::cout << "Found union of cross sections with " << union_components.size() << " components.\n";

    return union_components;
  }

  class VolumeMeshBuilder
  {
  public:
    using Iterator = std::vector<Slice>::iterator;
    using ConstIterator = std::vector<Slice>::const_iterator;

  private:
    // Grid points, including start and end points.
    std::array<std::vector<Slice>, 3> slices_per_axis;

  public:
    VolumeMeshBuilder() = default;

    void add_slice(Axis axis, double coord)
    {
      add_slice(axis, Slice(coord, 0));
    }

    VolumeMesh build(Mesh& mesh, std::uint64_t levels)
    {
      for(int i = 0; i < 6; i++)
      {
        fix_aspect_ratios(Axis::X);
        fix_aspect_ratios(Axis::Y);
        fix_aspect_ratios(Axis::Z);
      }

      /*
      while(slices_per_axis[0].size() > 3 && slices_per_axis[1].size() > 3 && slices_per_axis[2].size() > 3 && num_cells() > target_cells)
      {
        coarsen(Axis::X);
        coarsen(Axis::Y);
        coarsen(Axis::Z);
      }
      */

      choose_subdivision_levels(mesh, levels);

      return {
        slices_per_axis[0].begin(),
        slices_per_axis[0].end(),
        slices_per_axis[1].begin(),
        slices_per_axis[1].end(),
        slices_per_axis[2].begin(),
        slices_per_axis[2].end()};
    }

  private:
    void add_slice(Axis axis, Slice slice)
    {
      std::vector<Slice>& slices = slices_per_axis.at(int(axis));

      slices.insert(
        std::upper_bound(
          slices.begin(),
          slices.end(),
          slice,
          [](const Slice& a, const Slice& b) { return a.coord < b.coord; }),
        slice);
    }

    std::size_t num_cells() const
    {
      return (slices_per_axis[0].size() - 1) * (slices_per_axis[1].size() - 1) * (slices_per_axis[2].size() - 1);
    }

    std::pair<double, double> slice_extrema(Axis axis)
    {
      double max = 0;
      double min = std::numeric_limits<double>::max();

      const auto& slices = slices_per_axis.at(int(axis));
      for(std::size_t slice(0); slice < slices.size() - 1; slice++)
      {
        const double diff = slices[slice + 1].coord - slices[slice].coord;
        min = std::min(min, diff);
        max = std::max(max, diff);
      }
      return std::make_pair(min, max);
    }

    BoundingBox bounding_box()
    {
      BoundingBox result{};

      result.min.x = slices_per_axis[0].front().coord;
      result.min.y = slices_per_axis[1].front().coord;
      result.min.z = slices_per_axis[2].front().coord;

      result.max.x = slices_per_axis[0].back().coord;
      result.max.y = slices_per_axis[1].back().coord;
      result.max.z = slices_per_axis[2].back().coord;

      return result;
    }

    /// \brief Insert extra slices to fix segments with too large aspect ratios in z-direction
    void fix_aspect_ratios(Axis axis)
    {
      // Determine subdivision levels for all z-slices
      std::vector<Slice>& slices = slices_per_axis.at(int(axis));

      std::pair<double, double> x_extrema = slice_extrema(other_axes(axis)[0]);
      std::pair<double, double> y_extrema = slice_extrema(other_axes(axis)[1]);

      std::vector<double> new_slices;

      for(std::size_t slice(0); slice < slices.size() - 1; slice++)
      {
        const double width = slices[slice + 1].coord - slices[slice].coord;

        if(width < x_extrema.first && width < y_extrema.first)
        {
          continue;
        }

        const double aspect_ratio = std::max({
          width / x_extrema.first,
          x_extrema.first / width,
          width / y_extrema.first,
          y_extrema.first / width});
        const auto new_segments = static_cast<std::size_t>(aspect_ratio / 1.5);
        const double new_width = width / static_cast<double>(new_segments);

        if(new_segments > 0)
        {
          for(std::size_t i(0); i < new_segments - 1; i++)
          {
            new_slices.push_back(slices[slice].coord + (static_cast<double>(i + 1) * new_width));
          }
        }
      }

      for(double slice : new_slices)
      {
        add_slice(axis, slice);
      }
    }

    /// Delete every second slice
    void coarsen(Axis axis)
    {
      std::vector<Slice>& slices = slices_per_axis.at(int(axis));
      auto iter = slices.begin() + 2;

      // Don't delete last slice. Keep bounding box the same
      while(iter < slices.end() - 1)
      {
        iter = slices.erase(iter);
      }
    }

    void choose_subdivision_levels(Mesh& mesh, std::uint64_t levels)
    {
      // Find largest x and y lengths of any cell
      const double max_x_edge = slice_extrema(Axis::X).second;
      const double max_y_edge = slice_extrema(Axis::Y).second;

      // Determine subdivision levels for all z-slices
      std::vector<Slice>& z_slices = slices_per_axis[int(Axis::Z)];

      for(std::size_t slice(0); slice < z_slices.size() - 1; slice++)
      {
        BoundingBox slice_bb = bounding_box();
        slice_bb.min.z = z_slices[slice].coord;
        slice_bb.max.z = z_slices[slice + 1].coord;

        std::vector<std::pair<HexMesher::Point, double>> slice_gaps = gaps(mesh, slice_bb);

        double slice_min_gap = std::numeric_limits<double>::max();
        for(auto gap : slice_gaps)
        {
          slice_min_gap = std::min(gap.second, slice_min_gap);
        }

        double max_edge_length = std::max({slice_bb.max.z - slice_bb.min.z, max_x_edge, max_y_edge});

        // Each intended refinement level with divide the maximum edge length in two
        double max_refined_edge_length = max_edge_length / std::pow(2, levels);

        // Each adaptive refinement will divide edge length by three
        // max_edge_length * (1/3)^(levels + n) = slice_min_gap
        // log_1/3()

        std::uint64_t ref_level =
          std::uint64_t(std::max(0.0, std::ceil(std::log(max_refined_edge_length / slice_min_gap) / log(3.0))));
        ref_level = std::clamp<std::uint64_t>(ref_level, 0, 2);
        z_slices[slice].subdivision_level = std::max(ref_level, z_slices[slice].subdivision_level);
        z_slices[slice + 1].subdivision_level = ref_level;
      }

      // Insert buffer slices to keep transition elements contained

      std::vector<Slice> unbuffered_slices(z_slices.begin(), z_slices.end());

      // We want the transition elements to be as nice as possible, i.e.
      // have an aspect ratio of 1.
      // We thus match the buffer width to the x/y widths
      const double buffer_width = std::max(max_x_edge, max_y_edge);

      for(std::size_t slice(0); slice < unbuffered_slices.size() - 1; slice++)
      {
        const Slice& left_slice = unbuffered_slices[slice];
        const Slice& right_slice = unbuffered_slices[slice + 1];
        const std::uint64_t level_left = left_slice.subdivision_level;
        const std::uint64_t level_right = right_slice.subdivision_level;

        /*
        std::cout << "Slice: " << unsigned(slice) << "\n";
        std::cout << "[" << left_slice.coord << ", " << right_slice.coord << "]\n";
        std::cout << "Level left: " << unsigned(level_left) << "\n";
        std::cout << "Level right: " << unsigned(level_right) << "\n";
        */

        if(level_left == level_right)
        {
          // Same subdivision levels. No transition layer needed
          continue;
        }

        const bool rising = level_left < level_right;
        double start = rising ? right_slice.coord : left_slice.coord;
        double dir = rising ? -1.0 : 1.0;
        std::uint64_t start_level = rising ? level_right : level_left;

        // Transition from lower to higher subdivision level in this segment
        // Insert transition slices at end of segment

        const std::uint64_t required_transition_slices = rising ? level_right - level_left : level_left - level_right;
        const double available_width = (right_slice.coord - left_slice.coord);
        const auto available_transition_slices = std::uint64_t(std::floor(available_width / buffer_width));

        const std::uint64_t transition_slices = std::min(required_transition_slices, available_transition_slices);

        for(std::uint64_t i(0); i < transition_slices; i++)
        {
          add_slice(Axis::Z, Slice(start + (dir * static_cast<double>(i + 1) * buffer_width), start_level - i - 1));
        }
      }

      // We have now inserted all possible transition layers
      // NOTE(mmuegge): This assumes that we can't move slices.
      // That is a possible extension
      // As a last step smooth out the subdivision levels

      bool converged = false;

      while(!converged)
      {
        converged = true;
        for(std::size_t slice(0); slice < z_slices.size() - 1; slice++)
        {
          Slice& left_slice = z_slices[slice];
          Slice& right_slice = z_slices[slice + 1];
          const std::uint64_t level_left = left_slice.subdivision_level;
          const std::uint64_t level_right = right_slice.subdivision_level;

          const std::uint64_t diff = std::max(level_left, level_right) - std::min(level_left, level_right);

          if(diff > 1)
          {
            converged = false;
            if(level_left < level_right)
            {
              left_slice.subdivision_level = level_right - 1;
            }
            if(level_left > level_right)
            {
              right_slice.subdivision_level = level_left - 1;
            }
          }
        }
      }
    }
  };

  VolumeMesh fbm_mesh(Mesh& mesh, AABBTree& aabb_tree, std::uint64_t levels)
  {
    BoundingBox bb = bounding_box(mesh);
    const double mesh_size = std::max({bb.max.x - bb.min.x, bb.max.y - bb.min.y, bb.max.z - bb.min.z});
    std::vector<std::pair<HexMesher::Point, double>> gap_vec = gaps(mesh);
    std::vector<std::pair<HexMesher::Point2D, double>> depths = z_depths(mesh, aabb_tree);
    MinGap mg = min_gap(mesh);

    // Partition the mesh into areas that need the _same_ amount of refinements to hit the min-gap.
    // Because we are creating a coarse mesh, we do not want to create any cells smaller than the min-gap.
    // This allows us to discretize the gap information into global-min-gap sized buckets.

    const auto update_bucket = [](std::optional<double>& bucket, double value, const auto op)
    {
      if(!bucket.has_value())
      {
        bucket = value;
      }
      else
      {
        bucket = op(bucket.value(), value);
      }
    };

    const auto num_buckets_x = std::size_t(std::ceil((bb.max.x - bb.min.x) / mg.gap));
    const auto num_buckets_y = std::size_t(std::ceil((bb.max.y - bb.min.y) / mg.gap));
    const auto num_buckets_z = std::size_t(std::ceil((bb.max.z - bb.min.z) / mg.gap));

    std::vector<std::optional<double>> buckets_x(num_buckets_x);
    std::vector<std::optional<double>> buckets_y(num_buckets_y);
    std::vector<std::optional<double>> buckets_z(num_buckets_z);

    for(const std::pair<HexMesher::Point, double>& gap : gap_vec)
    {
      const std::size_t idx_z = std::clamp(std::size_t(std::floor((gap.first.z - bb.min.z) / mg.gap)), std::size_t(0), num_buckets_z - 1);

      std::optional<double>& bucket_z = buckets_z.at(idx_z);

      update_bucket(bucket_z, gap.second, [](auto a, auto b) { return std::min(a, b); });

      const std::size_t idx_x = std::clamp(std::size_t(std::floor((gap.first.x - bb.min.x) / mg.gap)), std::size_t(0), num_buckets_x - 1);
      const std::size_t idx_y = std::clamp(std::size_t(std::floor((gap.first.y - bb.min.y) / mg.gap)), std::size_t(0), num_buckets_y - 1);

      std::optional<double>& bucket_x = buckets_x.at(idx_x);
      std::optional<double>& bucket_y = buckets_y.at(idx_y);

      update_bucket(bucket_x, gap.second, [](auto a, auto b) { return std::min(a, b); });
      update_bucket(bucket_y, gap.second, [](auto a, auto b) { return std::min(a, b); });
    }

    // We can then partition those buckets into intervals, such that each interval needs a consistent
    // amount of refinements to hit the min-gap in the final mesh hierarchy.
    // The partitioning criterium is the ratio between the minimum and maximum min-gaps in the interval.
    // The easy solution is to keep that factor below 2, which works in any case.
    // A more complex solution considers regular and adaptive refinements, which divide the mesh size into 2 and 3 respectively.
    // We can determine the required number of refinements for an interval from the length of the interval and the minimum min-gap of the interval.
    // We also know how many regular refinements are intended.
    // We could thus keep the partitioning criterium as < 2 until we have surpassed the number of regular refinements.
    // From that point on we know that we adpatively refine at least once anyway and we can allow the minimum and maximum min-gap to differ by up to a factor of 3.
    // There are bucket for which we have no reliable min-gap information. We can distribute these buckets between adjacent intervals as we see fit.

    using BucketIterator = std::vector<std::optional<double>>::const_iterator;

    struct Interval
    {
      /// @brief First bucket of interval
      BucketIterator begin;
      /// @brief One after last bucket of interval
      BucketIterator end;

      std::optional<std::pair<double, double>> min_max = std::nullopt;

      bool can_expand(const std::optional<double>& bucket)
      {
        if(!min_max.has_value() && !bucket.has_value())
        {
          // Can expand "empty" intervals with empty buckets
          return true;
        }

        if(min_max.has_value() && bucket.has_value())
        {
          const double new_min = std::min(min_max.value().first, bucket.value());
          const double new_max = std::max(min_max.value().second, bucket.value());

          if(new_max / new_min < 2.0)
          {
            // Can expand with buckets that dont require more or less refinement
            return true;
          }
        }

        // Otherwise interval can't be expanded
        return false;
      }

      bool can_merge(const Interval& other)
      {
        if(!min_max.has_value() && !other.min_max.has_value())
        {
          return true;
        }

        if(min_max.has_value() && other.min_max.has_value())
        {
          const double new_min = std::min(min_max.value().first, other.min_max.value().first);
          const double new_max = std::max(min_max.value().second, other.min_max.value().second);

          if(new_max / new_min < 2.0)
          {
            return true;
          }
        }

        return false;
      }

      BucketIterator expand(BucketIterator b, BucketIterator e)
      {
        auto next = b;
        while(next != e && can_expand(*next))
        {
          end = next + 1;
          if(min_max.has_value())
          {
            const double new_min = std::min(min_max.value().first, next->value());
            const double new_max = std::max(min_max.value().second, next->value());

            min_max = std::make_pair(new_min, new_max);
          }
          next++;
        }

        return next;
      }
    };

    const auto make_intervals = [&](const std::vector<std::optional<double>>& buckets)
    {
      std::vector<Interval> intervals;

      auto current = buckets.begin();

      while(current != buckets.end())
      {
        Interval interval{current, current + 1};
        if(current->has_value())
        {
          interval.min_max = std::make_pair(current->value(), current->value());
        }

        current = interval.expand(current + 1, buckets.end());

        intervals.push_back(interval);
      }

      return intervals;
    };

    std::vector<Interval> intervals_x = make_intervals(buckets_x);
    std::vector<Interval> intervals_y = make_intervals(buckets_y);
    std::vector<Interval> intervals_z = make_intervals(buckets_z);


    // Distribute "empty" intervals. Merge adjacent intervals across empty if compatible, assign to smaller adjacent interval otherwise

    const auto merge_empty = [&](std::vector<Interval>& intervals)
    {
      auto iter = intervals.begin();

      while(iter != intervals.end())
      {
        if(iter->min_max.has_value() && std::distance(iter->begin, iter->end) > 1)
        {
          // Valid interval. Keep it
          iter++;
        }
        else
        {
          // "Empty" interval. Add to smaller of adjacent intervals
          auto prev = std::prev(iter);
          auto next = std::next(iter);

          bool prev_valid = iter > intervals.begin();
          bool next_valid = iter < std::prev(intervals.end());

          if(prev_valid && next_valid)
          {
            if(prev->can_merge(*next))
            {
              prev->end = next->end;
              intervals.erase(next);
            }
            else if(std::distance(prev->begin, prev->end) < std::distance(next->begin, next->end))
            {
              prev->end = iter->end;
            }
            else
            {
              next->begin = iter->begin;
            }
            iter = intervals.erase(iter);
          }
          else if(next_valid)
          {
            // Empty interval at start
            next->begin = iter->begin;
            iter = intervals.erase(iter);
          }
          else if(prev_valid)
          {
            // Empty interval at end
            prev->end = iter->end;
            iter = intervals.erase(iter);
          }
          else
          {
            // Empty interval is only interval. Keep it
            iter++;
          }
        }
      }
    };

    merge_empty(intervals_x);
    merge_empty(intervals_y);
    merge_empty(intervals_z);

    const auto merge_too_small = [&](std::vector<Interval>& intervals)
    {
      auto iter = intervals.begin();

      while(iter != intervals.end())
      {
        const Interval& interval = *iter;

        auto size = std::distance(interval.begin, interval.end);

        if(static_cast<double>(size) * mg.gap < 0.01 * mesh_size || size < 3)
        {
          // Too small. Merge into smaller of adjacent intervals
          auto prev = std::prev(iter);
          auto next = std::next(iter);

          bool prev_valid = iter > intervals.begin();
          bool next_valid = iter < std::prev(intervals.end());

          if(prev_valid && next_valid)
          {
            if(std::distance(prev->begin, prev->end) < std::distance(next->begin, next->end))
            {
              prev->end = iter->end;
            }
            else
            {
              next->begin = iter->begin;
            }
            iter = intervals.erase(iter);
          }
          else if(next_valid)
          {
            // Empty interval at start
            next->begin = iter->begin;
            iter = intervals.erase(iter);
          }
          else if(prev_valid)
          {
            // Empty interval at end
            prev->end = iter->end;
            iter = intervals.erase(iter);
          }
          else
          {
            // Empty interval is only interval. Keep it
            iter++;
          }
        }
        else
        {
          iter++;
        }
      }
    };

    merge_too_small(intervals_x);
    merge_too_small(intervals_y);
    merge_too_small(intervals_z);

    VolumeMeshBuilder builder;
    for(const Interval& interval : intervals_x)
    {
      builder.add_slice(Axis::X, bb.min.x + (mg.gap * static_cast<double>(std::distance(buckets_x.cbegin(), interval.begin))));
    }
    builder.add_slice(Axis::X, bb.max.x);

    for(const Interval& interval : intervals_y)
    {
      builder.add_slice(Axis::Y, bb.min.y + (mg.gap * static_cast<double>(std::distance(buckets_y.cbegin(), interval.begin))));
    }
    builder.add_slice(Axis::Y, bb.max.y);

    for(const Interval& interval : intervals_z)
    {
      builder.add_slice(Axis::Z, bb.min.z + (mg.gap * static_cast<double>(std::distance(buckets_z.cbegin(), interval.begin))));
    }
    builder.add_slice(Axis::Z, bb.max.z);

    return builder.build(mesh, levels);
  }
} // namespace HexMesher
