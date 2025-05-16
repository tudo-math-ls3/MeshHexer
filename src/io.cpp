#include <io.hpp>

namespace HexMesher
{
  void write_geo_compound_2d(const std::string& filename, const Polygon& poly)
  {
    std::ofstream output(filename);

    if(!output)
    {
      std::cerr << "Could not open file " << filename << " for writing.\n";
      return;
    }

    // Set options
    output << "Mesh.Algorithm = 8;\n"; // Frontal-Delaunay for Quads
    output << "Mesh.RecombinationAlgorithm = 3;\n"; // Blossom Full-Quad
    output << "Mesh.Format = 16;\n"; // vtk
    output << "Mesh.RecombineAll = 1;\n"; // Always produce quad meshes

    int next_tag = 1;

    for(const Point2D& p : poly)
    {
      output << "Point(" << next_tag++ << ") = {" << p.x() << ", " << p.y() << ", 0, 5.0};\n";
    }


    next_tag = 1;
    for(int i(0); i < poly.size() - 1; i++)
    {
      output << "Curve(" << next_tag++ << ") = {" << i + 1 << ", " << i + 2 << "};\n";
    }
    output << "Curve(" << next_tag++ << ") = {" << poly.size() << ", 1};\n";

    int line_loop_tag = next_tag;
    output << "Curve Loop(" << next_tag++ << ") = {";
    for(int i(0); i < poly.size(); i++)
    {
      output << (i + 1);
      if(i + 1 != poly.size())
      {
        output << ", ";
      }
    }
    output << "};\n";

    int surface_tag = next_tag;
    output << "Plane Surface(" << next_tag++ << ") = {" << line_loop_tag << "};\n"; 

    output << "Compound Curve{";
    for(int i(0); i < poly.size(); i++)
    {
      output << (i + 1);
      if(i + 1 != poly.size())
      {
        output << ", ";
      }
    }
    output << "};\n";
    output << "Compound Surface{" << surface_tag << "};\n";
    output << "Mesh 2\n";
  }


  void write_geo(const std::string& filename, const Polygon& poly)
  {
    std::ofstream output(filename);

    if(!output)
    {
      std::cerr << "Could not open file " << filename << " for writing.\n";
      return;
    }

    // Set options
    output << "Mesh.Algorithm = 8;\n"; // Frontal-Delaunay for Quads
    output << "Mesh.RecombinationAlgorithm = 3;\n"; // Blossom Full-Quad
    output << "Mesh.Format = 16;\n"; // vtk
    output << "Mesh.RecombineAll = 1;\n"; // Always produce quad meshes

    int next_tag = 1;

    for(const Point2D& p : poly)
    {
      output << "Point(" << next_tag++ << ") = {" << p.x() << ", " << p.y() << ", 0, 5.0};\n";
    }


    next_tag = 1;
    for(int i(0); i < poly.size() - 1; i++)
    {
      output << "Line(" << next_tag++ << ") = {" << i + 1 << ", " << i + 2 << "};\n";
    }
    output << "Line(" << next_tag++ << ") = {" << poly.size() << ", 1};\n";

    int line_loop_tag = next_tag;
    output << "Line Loop(" << next_tag++ << ") = {";
    for(int i(0); i < poly.size(); i++)
    {
      output << (i + 1);
      if(i + 1 != poly.size())
      {
        output << ", ";
      }
    }
    output << "};\n";

    int surface_tag = next_tag;
    output << "Plane Surface(" << next_tag++ << ") = {" << line_loop_tag << "};\n"; 

    //output << "Extrude {0, 0, 10} {Surface{" << surface_tag << "}; Layers {1}; }\n";
    output << "Mesh 2\n";
  }

  void write_polygon_brep(const std::string& filename, const Polygon& poly)
  {
    std::ofstream output(filename);

    if(!output)
    {
      std::cerr << "Could not open file " << filename << " for writing.\n";
      return;
    }

    output << "DBRep_DrawableShape\n\n";
    output << "CASCADE Topology V1, (c)  Matra-Datavision\n";
    output << "Locations 0\n";
    output << "Curve2ds 0\n";
    output << "Curves 0\n";
    output << "Polygon3D 1\n"; // One 3D polygon
    output << std::to_string(poly.size()) << " 0\n"; // Number of points and parameter presence. 0 = no parameters.
    output << "0.1\n"; // deflection

    // Points of the polygon, all on a single line, in triplets
    for(const Point2D& p: poly)
    {
      output << p.x() << " " << p.y() << " 0 ";
    }
    output << "\n";

    output << "PolygonOnTriangulations 0\n";
    output << "Surfaces 0\n";
    output << "Triangulations 0\n";
    output << "TShapes 0\n";
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
      output << poly.size() + 1 << "\n";
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
}
