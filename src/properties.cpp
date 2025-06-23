#include <properties.hpp>

#include <queue>

#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_tree.h>

#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/locate.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>

namespace HexMesher::Intern
{
  Vector3D surface_normal(Mesh& mesh, FaceIndex f, Point3D point)
  {
    auto maybe_normals_map = mesh.property_map<VertexIndex, Vector3D>("v:normals");
    if(!maybe_normals_map.has_value())
    {
      // TODO: Should be an assert
      std::cout << "surface_normals failed: missing property v:normals\n";
      return Vector3D(0.0, 0.0, 0.0);
    }
    Mesh::Property_map<VertexIndex, Vector3D> vertex_normals = maybe_normals_map.value();

    // Get barycentric coordinates of point
    auto location = CGAL::Polygon_mesh_processing::locate_in_face(point, f, mesh);

    std::array<Vector3D, 3> normals;

    // Get vertex normals
    int i(0);
    for(VertexIndex v : mesh.vertices_around_face(mesh.halfedge(f)))
    {
      normals[i++] = vertex_normals[v];
    }

    // Interpolate
    Vector3D result =
      location.second[0] * normals[0] + location.second[1] * normals[1] + location.second[2] * normals[2];
    return result;
  }
} // namespace HexMesher::Intern

namespace HexMesher
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  bool is_wound_consistently(const Mesh& mesh)
  {
    for(HalfedgeIndex idx : mesh.halfedges())
    {
      VertexIndex source = mesh.source(idx);
      VertexIndex target = mesh.target(idx);

      HalfedgeIndex opposite = mesh.opposite(idx);

      if(opposite != mesh.null_halfedge())
      {
        // If the mesh is wound consistently, then opposite halfedges need
        // to visit the same vertices in the opposite order.
        // Otherwise the two faces are wound differently.
        if(source != mesh.target(opposite) || target != mesh.source(opposite))
        {
          return false;
        }
      }
    }
    return true;
  }

  // Compute the maximum dihedral angle for each face in radians
  void compute_max_dihedral_angle(Mesh& mesh)
  {
    Mesh::Property_map<FaceIndex, double> dihedral_angles =
      mesh.add_property_map<FaceIndex, double>("f:dihedral_angle", 0).first;

    for(FaceIndex current : mesh.faces())
    {
      double max_angle = 0;
      Vector3D normal = PMP::compute_face_normal(current, mesh);

      for(VertexIndex v : mesh.vertices_around_face(mesh.halfedge(current)))
      {
        for(FaceIndex neighbor : mesh.faces_around_target(mesh.halfedge(v)))
        {
          if(neighbor > mesh.num_faces())
          {
            // NOTE(mmuegge): The inner loop produces invalid indices for the 10062 mesh
            continue;
          }
          Vector3D n = PMP::compute_face_normal(neighbor, mesh);

          double cosine = CGAL::scalar_product(normal, n);
          cosine = std::max(std::min(cosine, 1.0), 0.0);

          max_angle = std::max(max_angle, std::acos(cosine));
        }
      }

      dihedral_angles.begin();

      dihedral_angles[current] = max_angle;
    }
  }

  void compute_vertex_normals(Mesh& mesh)
  {
    Mesh::Property_map<VertexIndex, Vector3D> normals =
      mesh.add_property_map<VertexIndex, Vector3D>("v:normals", Vector3D(0.0, 0.0, 0.0)).first;

    CGAL::Polygon_mesh_processing::compute_vertex_normals(mesh, normals);
  }

  void maximal_inscribed_spheres(Mesh& mesh, const AABBTree& aabb_tree)
  {
    namespace PMP = CGAL::Polygon_mesh_processing;

    if(mesh.num_faces() == 0)
    {
      // No work to be done
      return;
    }

    // Create mesh property for storing thicknesses
    Mesh::Property_map<FaceIndex, double> diameter_property =
      mesh.add_property_map<FaceIndex, double>("f:MIS_diameter", 0).first;

    Mesh::Property_map<FaceIndex, std::uint32_t> id_property =
      mesh.add_property_map<FaceIndex, std::uint32_t>("f:MIS_id", 0).first;

    Mesh::Property_map<FaceIndex, double> similarity =
      mesh.add_property_map<FaceIndex, double>("f:similarity_of_normals", 0).first;

    Mesh::Property_map<FaceIndex, int> iters = mesh.add_property_map<FaceIndex, int>("f:MIS_iters", 0).first;

    // Determine bounding box of mesh
    Real max_radius = std::min({
                        aabb_tree.bbox().xmax() - aabb_tree.bbox().xmin(),
                        aabb_tree.bbox().ymax() - aabb_tree.bbox().ymin(),
                        aabb_tree.bbox().zmax() - aabb_tree.bbox().zmin(),
                      }) /
                      Real(2.0);

    // Determine maximum edge length of mesh
    Real max_edge_length(0);
    Real average_edge_length(0);
    for(auto e : mesh.edges())
    {
      max_edge_length = std::max(max_edge_length, PMP::edge_length(e, mesh));
      average_edge_length += PMP::edge_length(e, mesh);
    }
    average_edge_length /= Real(mesh.num_edges());

    Real normal_direction = PMP::is_outward_oriented(mesh) ? Real(-1.0) : Real(1.0);

    for(FaceIndex face_index : mesh.faces())
    {
      // Determine centroid of face
      Point3D centroid(0.0, 0.0, 0.0);
      for(VertexIndex v : mesh.vertices_around_face(mesh.halfedge(face_index)))
      {
        centroid += Real(1.0 / 3.0) * Vector3D(Point3D(CGAL::Origin()), mesh.point(v));
      }

      // Determine average edge length around the face.
      // We later use this to judge if a new closest point is
      // actually a closer point or just a slightly closer
      // vertex of the discretization of the same inscribed sphere

      Real local_edge_length(0);
      Real local_max_edge_length(0);
      for(HalfedgeIndex e : mesh.halfedges_around_face(mesh.halfedge(face_index)))
      {
        local_edge_length += PMP::edge_length(e, mesh);
        local_max_edge_length = std::max(local_max_edge_length, PMP::edge_length(e, mesh));
      }
      local_edge_length /= Real(3.0);

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
      // Optimization(mmuegge): We know the barycentric coordinates of the centroid already.
      // No reason for surface_normal to recalculate them
      // Either offer a version of surface_normal that accepts a point in barycentric coordinates
      // or a centroid_normal function.
      Vector3D normal(Intern::surface_normal(mesh, face_index, centroid));
      normal /= std::sqrt(normal.squared_length());
      Vector3D inward_normal = normal_direction * normal;

      if(std::abs(CGAL::approximate_sqrt(inward_normal.squared_length()) - 1.0) > 1e-2)
      {
        std::cout << "Non-unit normal at face " << face_index << "\n";
        std::cout << "Vector: " << normal << "\n";
        std::cout << "Length: " << CGAL::approximate_sqrt(normal.squared_length()) << "\n";
        std::cout << "Aborting\n";
        exit(1);
      }

      // Cast ray
      Ray3D ray(centroid, inward_normal);
      auto skip = [=](FaceIndex idx) { return idx == face_index; };
      RayIntersection intersection = aabb_tree.first_intersection(ray, skip);

      if(intersection && std::holds_alternative<Point3D>(intersection.value().first))
      {
        const Point3D& p = std::get<Point3D>(intersection.value().first);
        Real distance = CGAL::approximate_sqrt(Vector3D(centroid, p).squared_length());

        // NOTE(mmuegge): In meshes with self-intersections it can occur that a raycast hits
        // a triangle distinct from the current face at distance 0.
        // In that case we want to keep the initial radius based on the bounding box.
        if(distance > 0.0)
        {
          initial_radius = std::min(initial_radius, distance / 2.0);
          initial_id = intersection.value().second;
        }
      }

      if(intersection && std::holds_alternative<Segment3D>(intersection.value().first))
      {
        const Segment3D& segment = std::get<Segment3D>(intersection.value().first);
        Real distance = CGAL::approximate_sqrt(Vector3D(centroid, segment.source()).squared_length());

        if(distance > 0.0)
        {
          initial_radius = std::min(initial_radius, distance / 2.0);
          initial_id = intersection.value().second;
        }
      }

      Real prev_radius = initial_radius + 1;
      Real radius = initial_radius;
      FaceIndex id = initial_id;
      Point3D closest_point;

      // Shrink sphere until change in radius becomes too small
      while(prev_radius - radius > 1e-6)
      {
        Point3D center = centroid + radius * inward_normal;
        auto closest_point_data = aabb_tree.closest_point_and_primitive(center);

        if(closest_point_data.second == face_index)
        {
          // Disregard closest points on self
          break;
        }

        closest_point = closest_point_data.first;

        Vector3D to_closest = Vector3D(centroid, closest_point);
        if(CGAL::scalar_product(inward_normal, to_closest) <= Real(0))
        {
          // Next closest point is behind the current face.
          // Disregard it
          break;
        }

        Real distance_to_closest = CGAL::approximate_sqrt(Vector3D(center, closest_point).squared_length());
        Vector3D closest_normal = -Intern::surface_normal(mesh, closest_point_data.second, closest_point_data.first);

        Real discretization_error(0);
        // NOTE(mmuegge): Must be local_max_edge_length < 2 * radius to be a valid triangle
        // Tests show 1.5 works better though
        if(local_max_edge_length < 1.5 * radius)
        {
          // Largest distance between sphere and an edge of maximum length, assuming the vertices of that edge lie on
          // the sphere Any closest points within this distance belong to the same radius, they are just slightly offset
          // due to the surface discretization.
          discretization_error =
            radius - CGAL::approximate_sqrt(radius * radius - std::pow(local_max_edge_length / 2.0, 2));
        }
        else
        {
          discretization_error = 0.02 * radius;
        }

        if(distance_to_closest >= radius - 2 * discretization_error)
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

        double d = CGAL::to_double(CGAL::approximate_sqrt((closest_point - centroid).squared_length()));
        double alpha =
          CGAL::to_double(CGAL::approximate_angle(inward_normal, Vector3D(centroid, closest_point))) / 180.0 * M_PI;

        Real new_radius;
        if(alpha != 0)
        {
          new_radius = Real(d * std::sin(alpha) / std::sin(M_PI - 2.0 * alpha));
        }
        else
        {
          new_radius = d;
        }

        prev_radius = radius;
        radius = std::min(new_radius, radius);
        id = closest_point_data.second;

        iters[face_index]++;
      }

      // Record thickness
      diameter_property[face_index] = CGAL::to_double(2.0 * radius);
      id_property[face_index] = id;

      // Record abs(cos(theta))
      // similarity[face_index] = std::abs(CGAL::to_double(CGAL::scalar_product(inward_normal, surface_normal(mesh, id,
      // closest_point))));
      similarity[face_index] = CGAL::scalar_product(inward_normal, -PMP::compute_face_normal(id, mesh));
    }
  }

  // Set up priority queue for choosing next vertex to explore
  struct FrontierEntry
  {
    VertexIndex idx;
    double priority;

    FrontierEntry() : idx(0), priority(0)
    {
    }

    FrontierEntry(VertexIndex v, double p) : idx(v), priority(p)
    {
    }
  };

  // The priority_queue returns the last element of the defined ordering first.
  // We want the last element to be the one with the lowest estimated distance.
  // Thus we sort in descending order.
  struct Comparator
  {
    bool operator()(FrontierEntry a, FrontierEntry b)
    {
      return a.priority > b.priority;
    }
  };

  double topological_distance(
    FaceIndex a,
    FaceIndex b,
    const Mesh& mesh,
    const std::vector<double>& edge_lengths,
    std::vector<double>& distances,
    double max_distance = 0.0)
  {
    // Prepare scratch memory. -1.0 is sentinel value for an unvisited node.
    for(auto& entry : distances)
    {
      entry = -1.0;
    }

    // Allocate frontier
    std::priority_queue<FrontierEntry, std::vector<FrontierEntry>, Comparator> frontier;

    // Initialise starting points.
    for(VertexIndex v : mesh.vertices_around_face(mesh.halfedge(a)))
    {
      // Add start point to frontier
      frontier.push(FrontierEntry(v, 0));

      // Start point has best known distance 0.0
      distances[v] = 0.0;
    }

    // Determine end points
    std::vector<VertexIndex> goal_vertices;
    for(VertexIndex v : mesh.vertices_around_face(mesh.halfedge(b)))
    {
      goal_vertices.push_back(v);
    }

    // Set up heuristic
    // Heuristic is shortest direct distance to any of the goal vertices
    auto heuristic = [&](VertexIndex v)
    {
      double result = std::numeric_limits<double>::max();
      for(VertexIndex g : goal_vertices)
      {
        const Point3D& a = mesh.point(v);
        const Point3D& b = mesh.point(g);
        double distance = CGAL::to_double(CGAL::approximate_sqrt(Vector3D(a, b).squared_length()));
        result = std::min(result, distance);
      }
      return result;
    };

    while(!frontier.empty())
    {
      // Remove next best element from the frontier
      VertexIndex current = frontier.top().idx;
      frontier.pop();

      if(std::find(goal_vertices.begin(), goal_vertices.end(), current) != goal_vertices.end())
      {
        // We found one of the target vertices.
        return distances[current];
      }

      for(HalfedgeIndex hedge : mesh.halfedges_around_target(mesh.halfedge(current)))
      {
        VertexIndex neighbor = mesh.source(hedge);

        double edge_length = edge_lengths[mesh.edge(hedge)];
        double new_cost = distances[current] + edge_length;

        if(max_distance != 0.0 && new_cost > max_distance)
        {
          continue;
        }

        if(distances[neighbor] == -1.0 || new_cost < distances[neighbor])
        {
          distances[neighbor] = new_cost;

          double priority = new_cost + heuristic(neighbor);
          frontier.push(FrontierEntry(neighbor, priority));
        }
      }
    }

    return -1.0;
  }

  void topological_distances(Mesh& mesh, const std::string& targets_property, double max_distance)
  {
    // Get targets
    auto maybe_target_map = mesh.property_map<FaceIndex, std::uint32_t>(targets_property);
    if(!maybe_target_map.has_value())
    {
      return;
    }
    Mesh::Property_map<FaceIndex, std::uint32_t> targets = maybe_target_map.value();

    // Get output property
    Mesh::Property_map<FaceIndex, double> topo_distance =
      mesh.add_property_map<FaceIndex, double>("f:topological_distance", 0).first;

    // Pre-calculate edge lengths
    std::vector<double> edge_lengths;
    for(EdgeIndex idx : mesh.edges())
    {
      Point3D a = mesh.point(mesh.target(mesh.halfedge(idx)));
      Point3D b = mesh.point(mesh.source(mesh.halfedge(idx)));
      edge_lengths.push_back(std::sqrt((a - b).squared_length()));
    }

    // Allocate scratch memory for distances
    std::vector<double> scratch_distances(mesh.num_vertices(), -1.0);

    for(FaceIndex f : mesh.faces())
    {
      topo_distance[f] =
        topological_distance(f, FaceIndex(targets[f]), mesh, edge_lengths, scratch_distances, max_distance);
    }
  }

  void topological_distances(Mesh& mesh, const std::string& targets_property, const std::string& max_distance_property)
  {
    // Get targets
    auto maybe_target_map = mesh.property_map<FaceIndex, std::uint32_t>(targets_property);
    if(!maybe_target_map.has_value())
    {
      return;
    }
    Mesh::Property_map<FaceIndex, std::uint32_t> targets = maybe_target_map.value();

    // Get maximum search distances
    auto maybe_max_distances_map = mesh.property_map<FaceIndex, double>(max_distance_property);
    if(!maybe_max_distances_map.has_value())
    {
      return;
    }
    Mesh::Property_map<FaceIndex, double> max_distances = maybe_max_distances_map.value();

    // Get output property
    Mesh::Property_map<FaceIndex, double> topo_distance =
      mesh.add_property_map<FaceIndex, double>("f:topological_distance", 0).first;

    // Pre-calculate edge lengths
    std::vector<double> edge_lengths;
    for(EdgeIndex idx : mesh.edges())
    {
      Point3D a = mesh.point(mesh.target(mesh.halfedge(idx)));
      Point3D b = mesh.point(mesh.source(mesh.halfedge(idx)));
      edge_lengths.push_back(std::sqrt((a - b).squared_length()));
    }

    double min_x;
    double min_y;
    double min_z;

    double max_x;
    double max_y;
    double max_z;

    for(const HexMesher::Point3D& p : mesh.points())
    {
      min_x = std::min(CGAL::to_double(p.x()), min_x);
      min_y = std::min(CGAL::to_double(p.y()), min_y);
      min_z = std::min(CGAL::to_double(p.z()), min_z);

      max_x = std::max(CGAL::to_double(p.x()), max_x);
      max_y = std::max(CGAL::to_double(p.y()), max_y);
      max_z = std::max(CGAL::to_double(p.z()), max_z);
    }

    double max_mesh_delta = std::max({max_x - min_x, max_y - min_y, max_z - min_z});

    // Allocate scratch memory for distances
    std::vector<double> scratch_distances(mesh.num_vertices(), -1.0);

    for(FaceIndex f : mesh.faces())
    {
      double max_distance = M_PI * max_distances[f];
      max_distance = std::max(max_distance, 0.001 * max_mesh_delta);
      topo_distance[f] =
        topological_distance(f, FaceIndex(targets[f]), mesh, edge_lengths, scratch_distances, max_distance);
    }
  }

  void score_gaps(Mesh& mesh)
  {
    HexMesher::Mesh::Property_map<HexMesher::FaceIndex, double> diameters =
      mesh.property_map<HexMesher::FaceIndex, double>("f:MIS_diameter").value();

    HexMesher::Mesh::Property_map<HexMesher::FaceIndex, double> topo_dists =
      mesh.property_map<HexMesher::FaceIndex, double>("f:topological_distance").value();

    HexMesher::Mesh::Property_map<HexMesher::FaceIndex, double> similarity_of_normals =
      mesh.property_map<HexMesher::FaceIndex, double>("f:similarity_of_normals").value();

    HexMesher::Mesh::Property_map<HexMesher::FaceIndex, double> dihedral_angles =
      mesh.property_map<HexMesher::FaceIndex, double>("f:dihedral_angle").value();

    HexMesher::Mesh::Property_map<HexMesher::FaceIndex, std::uint32_t> ids =
      mesh.property_map<HexMesher::FaceIndex, std::uint32_t>("f:MIS_id").value();

    // NOTE(mmuegge): Tried to store these as HexMesher::Points,
    // but CGAL produced a stack overflow on cleanup,
    // because it reuses values of the points in the exact kernel
    // and cleans them up recursively.
    // Try again once we support the inexact kernel.
    double min_x;
    double min_y;
    double min_z;

    double max_x;
    double max_y;
    double max_z;

    for(const HexMesher::Point3D& p : mesh.points())
    {
      min_x = std::min(CGAL::to_double(p.x()), min_x);
      min_y = std::min(CGAL::to_double(p.y()), min_y);
      min_z = std::min(CGAL::to_double(p.z()), min_z);

      max_x = std::max(CGAL::to_double(p.x()), max_x);
      max_y = std::max(CGAL::to_double(p.y()), max_y);
      max_z = std::max(CGAL::to_double(p.z()), max_z);
    }

    double max_diameter = std::min({max_x - min_x, max_y - min_y, max_z - min_z});

    auto topo_dist_score = [&](FaceIndex f)
    {
      double dist = topo_dists[f];

      // If the topological_distance is zero, this is not a valid gap
      if(dist == 0.0)
        return 0.0;

      // If the faces are either unconnected or farther apart than the search radius,
      // then this is an ideal gap
      if(dist == -1.0)
        return 1.0;

      // Divide by maximum search range for this face
      double rel = topo_dists[f] / std::max(0.001 * max_diameter, M_PI * diameters[f]);

      // At a relative distance of 0.5 the limiting face is (at best) at the opposite end of the MIS.
      // Penalize all distances below that.
      if(rel < 0.5)
        return 4 * rel * rel;

      // Above a relative distance of 0.5 there must be some amount of space between the MIS and the mesh.
      // All these gaps are ok.
      return 1.0;
    };

    auto normalized_diameter_score = [&](FaceIndex f)
    {
      double normalized = diameters[f] / max_diameter;

      // Penalize any gaps smaller than a tenth of a percent of the mesh size
      if(normalized < 1e-3)
      {
        return 1e3 * normalized;
      }

      // Full score for anything else
      return 1.0;
    };

    auto similarity_of_normals_score = [&](FaceIndex f)
    {
      if(-similarity_of_normals[f] < 0.0)
      {
        return 0.0;
      }

      return std::pow(-similarity_of_normals[f], 1.0 / 3.0);
    };

    auto dihedral_angle_score = [&](FaceIndex f)
    {
      if(dihedral_angles[f] < 0.3)
      {
        return 1.0;
      }
      return std::sqrt(1.0 - ((dihedral_angles[f] + 0.3) / M_PI));
    };

    auto max_edge_length = [&](FaceIndex f)
    {
      double result = 0.0;
      for(HalfedgeIndex idx : mesh.halfedges_around_face(mesh.halfedge(f)))
      {
        result = std::max(result, PMP::edge_length(idx, mesh));
      }
      return result;
    };

    std::vector<std::pair<HexMesher::FaceIndex, HexMesher::FaceIndex>> self_intersections;
    CGAL::Polygon_mesh_processing::self_intersections(mesh, std::back_inserter(self_intersections));

    Mesh::Property_map<FaceIndex, double> gap_score = mesh.add_property_map<FaceIndex, double>("f:gap_score", 0).first;

    for(FaceIndex f : mesh.faces())
    {
      FaceIndex limiting = FaceIndex(ids[f]);

      double edge_ratio = max_edge_length(f) / max_edge_length(limiting);
      if(edge_ratio < 1.0)
      {
        edge_ratio = 1.0 / edge_ratio;
      }

      if(
        PMP::face_aspect_ratio(f, mesh) < 5.0 && PMP::face_aspect_ratio(limiting, mesh) < 5.0 &&
        std::find(self_intersections.begin(), self_intersections.end(), std::make_pair(f, limiting)) ==
          self_intersections.end() &&
        std::find(self_intersections.begin(), self_intersections.end(), std::make_pair(limiting, f)) ==
          self_intersections.end() &&
        diameters[f] / max_diameter > 1e-4 && edge_ratio < 5.0)
      {
        if(topo_dists[f] == -1.0 && topo_dists[FaceIndex(ids[f])] == -1.0)
        {
          // Perfect gap if both triangles are unconnected.
          // NOTE(mmuegge): There are situations where a large triangle
          // hits a very small triangle across a corner of the mesh.
          // In that case the small triangle might not find its limiting triangle
          // within the search radius for topological distances, but the large
          // triangle might, because its vertices are further from its centroid,
          // giving an initial distance boost.
          // We thus check that neither triangle found its limiting triangle,
          // to determine if we should score the gap more exactly.
          // Note that the above condition does not imply f = ids[f].
          // TODO(mmuegge): Can we fix this by considering the distances
          // from the centroid to a starting vertex in the topological distance search?
          gap_score[f] = 1.0;
        }
        else
        {
          gap_score[f] = (1.0 / 3.0) * (topo_dist_score(f) + similarity_of_normals_score(f) + dihedral_angle_score(f));
        }
      }
      else
      {
        gap_score[f] = 0;
      }
    }
  }

  MinGap min_gap_percentile(Mesh& mesh, double percentile)
  {
    HexMesher::Mesh::Property_map<HexMesher::FaceIndex, double> diameters =
      mesh.property_map<HexMesher::FaceIndex, double>("f:MIS_diameter").value();

    HexMesher::Mesh::Property_map<HexMesher::FaceIndex, double> scores =
      mesh.property_map<HexMesher::FaceIndex, double>("f:gap_score").value();

    HexMesher::Mesh::Property_map<HexMesher::FaceIndex, std::uint32_t> ids =
      mesh.property_map<HexMesher::FaceIndex, std::uint32_t>("f:MIS_id").value();

    std::vector<double> percentiles(scores.begin(), scores.end());

    // Sort scores in descending order
    std::sort(percentiles.begin(), percentiles.end(), [](double a, double b) { return a > b; });

    /*
    auto last = std::unique(percentiles.begin(), percentiles.end(), [](double a, double b)
    {
      return std::abs(a -b) < 1e-6;
    });
    percentiles.erase(last, percentiles.end());
    */

    std::size_t percentile_size = std::ceil(Real(percentiles.size()) * Real(percentile));
    double cutoff_score = *(percentiles.begin() + percentile_size);

    std::cout << "Score cutoff is " << cutoff_score << "\n";

    MinGap result;
    result.gap = std::numeric_limits<double>::max();

    for(FaceIndex f : mesh.faces())
    {
      if(scores[f] < cutoff_score)
      {
        continue;
      }
      if(diameters[f] < result.gap)
      {
        result.gap = diameters[f];
        result.origin = f;
        result.limiting = FaceIndex(ids[f]);
      }
    }

    return result;
  }

  MinGap min_gap(Mesh& mesh)
  {
    HexMesher::Mesh::Property_map<HexMesher::FaceIndex, double> diameters =
      mesh.property_map<HexMesher::FaceIndex, double>("f:MIS_diameter").value();

    HexMesher::Mesh::Property_map<HexMesher::FaceIndex, double> scores =
      mesh.property_map<HexMesher::FaceIndex, double>("f:gap_score").value();

    HexMesher::Mesh::Property_map<HexMesher::FaceIndex, std::uint32_t> ids =
      mesh.property_map<HexMesher::FaceIndex, std::uint32_t>("f:MIS_id").value();

    double threshold = 0.95;

    std::vector<FaceIndex> candidates;
    std::copy_if(
      mesh.faces_begin(),
      mesh.faces_end(),
      std::back_inserter(candidates),
      [&](FaceIndex f) { return scores[f] > threshold; });

    if(candidates.empty())
    {
      // No gaps with score above 0.95. Take 10th percentile instead.

      std::vector<double> percentiles(scores.begin(), scores.end());
      // Sort scores in descending order
      std::sort(percentiles.begin(), percentiles.end(), [](double a, double b) { return a > b; });

      std::size_t percentile_size = std::ceil(Real(percentiles.size()) * Real(0.1));
      double cutoff_score = *(percentiles.begin() + percentile_size);

      std::copy_if(
        mesh.faces_begin(),
        mesh.faces_end(),
        std::back_inserter(candidates),
        [&](FaceIndex f) { return scores[f] > cutoff_score; });
    }

    MinGap result;
    result.gap = std::numeric_limits<double>::max();

    for(FaceIndex f : candidates)
    {
      if(diameters[f] < result.gap)
      {
        result.gap = diameters[f];
        result.origin = f;
        result.limiting = FaceIndex(ids[f]);
      }
    }

    return result;
  }
} // namespace HexMesher
