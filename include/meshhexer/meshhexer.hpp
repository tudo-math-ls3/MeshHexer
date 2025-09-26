#pragma once

#include <meshhexer/types.hpp>

#include <memory>
#include <optional>
#include <variant>

namespace MeshHexer
{
  /**
   * \brief Tagged union for error-handling. Contains either a success value of type T or an error value of type E
   */
  template<typename T, typename E>
  class Result
  {
    /**
     * \brief Type-tag for constructing a success-result.
     */
    struct TagOK
    {
    };

    /**
     * \brief Type-tag for constructing an error-result.
     */
    struct TagErr
    {
    };

    /// Result content
    std::variant<T, E> content;


    /// Ccopy constructor for success-results
    Result(TagOK /*tag*/, const T& value) : content(std::in_place_index<0>, value)
    {
    }

    /// Move constructor for success-results
    Result(TagOK /*tag*/, T&& value) : content(std::in_place_index<0>, std::move(value))
    {
    }

    /// Copy constructor for error-results
    Result(TagErr /*tag*/, const E& value) : content(std::in_place_index<1>, value)
    {
    }

    /// Move constructor for error-results
    Result(TagErr /*tag*/, E&& value) :
      content(std::in_place_index<1>, std::move(value))
    {
    }

  public:

    // Rule of five
    Result() = delete;

    /// Deleted copy constructor
    Result(const Result&) = delete;

    /// Deleted copy assignment
    Result& operator=(const Result&) = delete;

    /// Move constructor
    Result(Result&&) = default;

    /// Move-assignment
    Result& operator=(Result&&) = default;

    /// Destructor
    ~Result() = default;

    // Static factory methods

    /**
     * \brief Create a success-result from the given value
     *
     * \param[in] value Value of success
     */
    static Result ok(const T& value)
    {
      return Result(TagOK{}, value);
    }

    /**
     * \brief Create a success-result from the given value
     *
     * \param[in] value Value of success
     */
    static Result ok(T&& value)
    {
      return Result(TagOK{}, std::move(value));
    }

    /**
     * \brief Create a error-result from the given value
     *
     * \param[in] value Value of error
     */
    static Result err(const E& value)
    {
      return Result(TagErr{}, value);
    }

    /**
     * \brief Create a error-result from the given value
     *
     * \param[in] value Value of error
     */
    static Result err(E&& value)
    {
      return Result(TagErr{}, std::move(value));
    }

    /// Returns true if this is a success-result
    bool is_ok() const
    {
      return content.index() == 0;
    }

    /// Returns true if this is an error-result
    bool is_err() const
    {
      return content.index() == 1;
    }

    /**
     * \brief Success value accessor
     *
     * \returns A copy of the stored value. Errors if this is not a success-result.
     */
    T ok_value() const
    {
      return std::get<0>(content);
    }

    /**
     * \brief Error value accessor
     *
     * \returns A copy of the stored value. Errors if this is not an error-result.
     */
    E err_value() const
    {
      return std::get<1>(content);
    }

    /**
     * \brief Success reference accessor
     *
     * \returns A reference to the stored value. Errors if this is not a success-result.
     */
    T& ok_ref()
    {
      return std::get<0>(content);
    }

    /**
     * \brief Error reference accessor
     *
     * \returns A reference to the stored value. Errors if this is not an error-result.
     */
    E& err_ref()
    {
      return std::get<1>(content);
    }

    /**
     * \brief Success reference accessor
     *
     * \returns A reference to the stored value. Errors if this is not a success-result.
     */
    const T& ok_ref() const
    {
      return std::get<0>(content);
    }

    /**
     * \brief Error reference accessor
     *
     * \returns A reference to the stored value. Errors if this is not an error-result.
     */
    const E& err_ref() const
    {
      return std::get<1>(content);
    }

    /**
     * \brief Success move accessor
     *
     * \returns A rvalue-reference to the stored value. Errors if this is not a success-result.
     */
    T&& take_ok() &&
    {
      return std::get<0>(std::move(content));
    };

    /**
     * \brief Error move accessor
     *
     * \returns A rvalue-reference to the stored value. Errors if this is not an error-result.
     */
    E&& take_err() &&
    {
      return std::get<1>(std::move(content));
    };

    /**
     * \brief Operator bool overload
     *
     * \returns True if parsing was succesful
     */
    explicit operator bool() const
    {
      return is_ok();
    }
  };

  /**
   * \brief Tagged union for error-handling. Contains either a success value of type T or an error value of type E
   *
   * This is a specialization for results containing either nothing or an error.
   */
  template<typename E>
  class Result<void, E>
  {
    /// Optional error value
    std::optional<E> error;

  public:
    // Rule of five
    Result() : error(std::nullopt)
    {
    }

    /// Constructor
    explicit Result(const E& value) : error(value)
    {
    }

    /// Constructor
    explicit Result(E&& value) : error(std::move(value))
    {
    }

    /// Copy constructor
    Result(const Result&) = delete;

    /// Copy-assignment
    Result& operator=(const Result&) = delete;

    /// Move-constructor
    Result(Result&&) = default;

    /// Move-assignment
    Result& operator=(Result&&) = default;

    /// Destructor
    ~Result() = default;

    // Static factory methods

    /**
     * \brief Create a success-result
     */
    static Result ok()
    {
      return Result();
    }

    /**
     * \brief Create an error-result from the given value
     */
    static Result err(const E& value)
    {
      return Result(value);
    }

    /**
     * \brief Create an error-result by moving the given value
     */
    static Result err(E&& value)
    {
      return Result(std::move(value));
    }

    /// Returns true if this is a success-result
    bool is_ok() const
    {
      return !error.has_value();
    }

    /// Returns true if this is an error-result
    bool is_err() const
    {
      return error.has_value();
    }

    /**
     * \brief Error value accessor
     *
     * \returns A copy of the stored value. Errors if this is not an error-result.
     */
    E err_value() const
    {
      return error.value();
    }

    /**
     * \brief Error reference accessor
     *
     * \returns A reference to the stored value. Errors if this is not an error-result.
     */
    E& err_ref()
    {
      return error.value();
    }

    /**
     * \brief Error reference accessor
     *
     * \returns A reference to the stored value. Errors if this is not an error-result.
     */
    const E& err_ref() const
    {
      return error.value();
    }

    /**
     * \brief Error move accessor
     *
     * \returns A rvalue-reference to the stored value. Errors if this is not an error-result.
     */
    E&& take_err() &&
    {
      return std::move(error).value();
    };

    /**
     * \brief Operator bool overload
     *
     * \returns True if parsing was succesful
     */
    explicit operator bool() const
    {
      return is_ok();
    }
  };

  /**
   * \brief MeshHexer surface mesh class
   *
   * Main API-entry-point for MeshHexer.
   * Contains a triangulated surface and can be used to either calculate properties of the surface mesh
   * or create volume meshes from the surface (non-aligned meshes only, for now).
   *
   * To create a SurfaceMesh instance, load a triangle mesh from disk using the MeshHexer::load_from_file() method.
   */
  class SurfaceMesh
  {
  public:
    // Forward declaration for PIMPL
    class SurfaceMeshImpl;

  private:

    /// Pointer to implementation
    std::unique_ptr<SurfaceMeshImpl> impl;

  public:
    /// Constructor
    explicit SurfaceMesh(std::unique_ptr<SurfaceMeshImpl> ptr);

    /// Destructor
    ~SurfaceMesh();

    /// Move-constructor
    SurfaceMesh(SurfaceMesh&&) noexcept;

    /// Copy-constructor
    SurfaceMesh(const SurfaceMesh&) = delete;

    /// Move-assign operator
    SurfaceMesh& operator=(SurfaceMesh&&) noexcept;

    /// Copy-assign operator
    SurfaceMesh& operator=(const SurfaceMesh&) = delete;

    /// standard ctor
    SurfaceMesh() = default;


    /**
     * \brief Axis-aligned bounding box of the mesh
     *
     * \returns The axis-aligned bounding box of the mesh
     */
    BoundingBox bounding_box() const;

    /// Returns the number of vertices in the surface mesh
    std::uint32_t num_vertices() const;
    /// Returns the number of edges in the surface mesh
    std::uint32_t num_edges() const;
    /// Returns the number of faces in the surface mesh
    std::uint32_t num_faces() const;

    /**
     * \brief Check if the surface mesh describes a closed surface
     *
     * \returns True, if there are no border edges in the surface mesh
     */
    bool is_closed() const;

    /**
     * \brief Check if the winding order of the faces of the surface mesh is consistent
     *
     * \return True, if either all faces of the surface mesh are wound clockwise or all faces of the surface mesh are wound counter-clockwise
     */
    bool is_wound_consistently() const;

    /**
     * \brief Check if the surface normals of all faces of the surface mesh point towards the unbounded side
     *
     * \pre The surface mesh is closed, call SurfaceMesh::is_closed() to check
     * \pre The surface mesh is wound consistently, call SurfaceMesh::is_wound_consistently() to check
     *
     * \returns True, if the surface normals of all faces of the surface mesh point towards the unbounded side, i.e. outside.
     */
    bool is_outward_oriented() const;

    /**
     * \brief Compute the smallest aspect ratio of any face of the surface mesh
     *
     * \returns The minimal aspect ratio of any face of the surface mesh
     */
    double minimal_aspect_ratio() const;

    /**
     * \brief Compute the largest aspect ratio of any face of the surface mesh
     *
     * \returns The maximum aspect ratio of any face of the surface mesh
     */
    double maximal_aspect_ratio() const;

    /**
     * \brief Compute Gap candidates of the surface mesh.
     *
     * Gaps try to capture the diameter of the volume described by the surface mesh.
     * The diameter is measured from the barycentric centers of all faces of the surface mesh.
     * For each face of the surface mesh the largest sphere that just touches the center of
     * the face and one other point of the surface mesh is computed.
     * We call this sphere the maximal inscribed sphere (MIS) and use the diameter of that sphere
     * as the diameter at that face.
     *
     * \verbatim
     * │                 │
     * │     ∙∙∙∙∙∙∙     │
     * │  ∙∙∙       ∙∙∙  │
     * │ ∙             ∙ │
     * │∙    Diameter   ∙│
     * │∙---------------∙│
     * │∙               ∙│
     * │ ∙             ∙ │
     * │  ∙∙∙       ∙∙∙  │
     * │     ∙∙∙∙∙∙∙     │
     * │                 │
     * \endverbatim
     *
     * Note that not all these MIS span gaps that correspond to the intuitive diameter of the
     * surface mesh at that point. Badly reconstructed meshes for example might contain many
     * small "valleys" on surfaces that are supposed to be "smooth"
     *
     * \verbatim
     *              ∙∙∙∙∙∙∙
     *            ∙∙       ∙∙
     *\         ∙           ∙         /
     *  \\      ∙             ∙      //
     *    \\    ∙             ∙    //
     *      \\  ∙             ∙  //
     *        \\ ∙           ∙ //
     *          \\∙∙       ∙∙//
     *            \\∙∙∙∙∙∙∙//
     *              \\   //
     *                \ /
     * \endverbatim
     *
     * Each gap is given a confidence score in the range [0, 1] that shows how likely that gap
     * corresponds to a real diameter. That score is determined from:
     * - the aspect ratios of involved faces
     * - self-intersections of the surface mesh
     * - relative sizes of involved faces
     * - similarity of normals at the touching points of the MIS
     * - topological distance (distance across the surface between the touching points of the MIS)
     * - dihedral angles near the touching points of the MIS
     *
     * Generally gaps with a score above 0.95 have a good change of corresponding to real diameters.
     *
     * \returns A vector containing one gap candiate for each face of the surface mesh
     */
    std::vector<Gap> gaps();

    /**
     * \brief Compute the smallest Gap of the surface mesh
     *
     * See SurfaceMesh::gaps() for details on how gaps area measured and scored.
     *
     * \returns The smallest gap, by diameter, of all gaps with a score above 0.95,
     * or the smallest gap among the top 10 percent of all gaps, if no such gaps exist.
     */
    Gap min_gap();

    /**
     * \brief Create a non-fitted volume mesh for this surface mesh.
     *
     * \param[in] levels Number of regular refinements in intended multi-grid hierarchy
     *
     * The volume mesh is created such that after the adaptive refinement prescribed
     * by the volume mesh and after \c levels further global refinement steps,
     * the final volume mesh roughly hits the local min gap of the surface mesh.
     *
     * Assumes that adaptive refinements are 3-refinements and global refinements are 2-refinements.
     */
    VolumeMesh fbm_mesh(std::uint64_t levels);

    /**
     * \brief Create warnings for this surface mesh
     *
     * Warns about:
     * - self-intersections,
     * - degenerate triangles,
     * - highly anisotropic triangles,
     *
     * \returns Lists of warnings, separated by types
     */
    MeshWarnings warnings() const;

    /**
     * \brief Write the surface mesh to disk
     *
     * \param[in] filename Filename to write to. Must end in .ply
     *
     * Writes the surface mesh to disk as a .ply file. The written file
     * contains mesh properties that have been calculated as intermediate results,
     * such as maximal inscribed spheres, or topological distances.
     * These properties will be reused if the written mesh is read again.
     *
     * \returns A result indicating success or containing an error message
     */
    Result<void, std::string> write_to_file(const std::string& filename);
  };

  /**
   * \brief Load a surface mesh from a file
   *
   * \param[in] filename     Path to mesh file
   * \param[in] triangulate  If set, input mesh will automatically be triangulated
   *
   * Supported file formats are
   * - .off
   * - .obj
   * - .stl
   * - .ply
   * - .ts
   * - .vtp
   *
   * \returns a result containing either the surface mesh or an error message
   */
  Result<SurfaceMesh, std::string>
  load_from_file(const std::string& filename, bool triangulate = false);

} // namespace MeshHexer
