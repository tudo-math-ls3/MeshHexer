#pragma once

#include <types.hpp>

#include <memory>
#include <optional>
#include <variant>

namespace HexMesher
{
  template<typename T, typename E>
  class Result
  {
    struct TagOK
    {
    };

    struct TagErr
    {
    };

    std::variant<T, E> content;

  public:
    // Rule of five
    Result() = delete;

    Result(TagOK, const T& value) : content(std::in_place_index<0>, value)
    {
    }

    Result(TagOK, T&& value) : content(std::in_place_index<0>, std::move(value))
    {
    }

    Result(TagErr, const E& value) : content(std::in_place_index<1>, value)
    {
    }

    Result(TagErr, E&& value) :
      content(std::in_place_index<1>, std::move(value))
    {
    }

    Result(const Result&) = delete;
    Result& operator=(const Result&) = delete;

    Result(Result&&) = default;
    Result& operator=(Result&&) = default;

    ~Result() = default;

    // Static factory methods
    static Result ok(const T& value)
    {
      return Result(TagOK{}, value);
    }

    static Result ok(T&& value)
    {
      return Result(TagOK{}, std::forward<T>(value));
    }

    static Result err(const E& value)
    {
      return Result(TagErr{}, value);
    }

    static Result err(E&& value)
    {
      return Result(TagErr{}, std::forward<E>(value));
    }

    bool is_ok() const
    {
      return content.index() == 0;
    }

    bool is_err() const
    {
      return content.index() == 1;
    }

    T ok_value() const
    {
      return std::get<0>(content);
    }

    E err_value() const
    {
      return std::get<1>(content);
    }

    T& ok_ref()
    {
      return std::get<0>(content);
    }

    E& err_ref()
    {
      return std::get<1>(content);
    }

    const T& ok_ref() const
    {
      return std::get<0>(content);
    }

    const E& err_ref() const
    {
      return std::get<1>(content);
    }

    T&& take_ok() &&
    {
      return std::get<0>(std::move(content));
    };

    E&& take_err() &&
    {
      return std::get<1>(std::move(content));
    };
  };

  template<typename E>
  class Result<void, E>
  {
    std::optional<E> error;

  public:
    // Rule of five
    Result() : error(std::nullopt)
    {
    }

    Result(const E& value) : error(value)
    {
    }

    Result(E&& value) : error(std::move(value))
    {
    }

    Result(const Result&) = delete;
    Result& operator=(const Result&) = delete;

    Result(Result&&) = default;
    Result& operator=(Result&&) = default;

    ~Result() = default;

    // Static factory methods
    static Result ok()
    {
      return Result();
    }

    static Result err(const E& value)
    {
      return Result(value);
    }

    static Result err(E&& value)
    {
      return Result(std::forward<E>(value));
    }

    bool is_ok() const
    {
      return !error.has_value();
    }

    bool is_err() const
    {
      return error.has_value();
    }

    E err_value() const
    {
      return error.value();
    }

    E& err_ref()
    {
      return error.value();
    }

    const E& err_ref() const
    {
      return error.value();
    }

    E&& take_err() &&
    {
      return std::move(error).value();
    };
  };

  class SurfaceMesh
  {
  public:
    // Forward declaration for PIMPL
    class SurfaceMeshImpl;

  private:
    std::unique_ptr<SurfaceMeshImpl> impl;

  public:
    explicit SurfaceMesh(std::unique_ptr<SurfaceMeshImpl> ptr);
    ~SurfaceMesh();

    SurfaceMesh(SurfaceMesh&&);
    SurfaceMesh(const SurfaceMesh&) = delete;

    SurfaceMesh& operator=(SurfaceMesh&&);
    SurfaceMesh& operator=(const SurfaceMesh&) = delete;

    BoundingBox bounding_box() const;

    std::uint32_t num_vertices() const;
    std::uint32_t num_edges() const;
    std::uint32_t num_faces() const;

    bool is_closed() const;
    bool is_wound_consistently() const;
    bool is_outward_oriented() const;
    double minimal_aspect_ratio() const;
    double maximal_aspect_ratio() const;

    MinGap min_gap();
    MinGap min_gap_percentile(double);

    VolumeMesh fbm_mesh(std::uint64_t levels);

    void warnings(MeshWarnings&) const;

    Result<void, std::string> write_to_file(const std::string& filename);
  };

  Result<SurfaceMesh, std::string>
  load_from_file(const std::string& filename, bool triangulate = false);

  // NOTE(mmuegge): This really belongs in the io header. Maybe introduce a public io header?
  template<typename Iter>
  void write_range_as_mtx(std::ostream& stream, Iter begin, Iter end)
  {
    // NOTE(mmuegge): For compatability with FEAT3 we write everything as
    // real. We could inspect the type produced by the iterator and set
    // integer as type if appropriate, but then we would also need to update
    // the parsing logic in FEAT3.
    stream << "%%MatrixMarket matrix array real general\n";

    const auto size = std::distance(begin, end);
    stream << size << " 1\n";

    for(Iter it = begin; it != end; it++)
    {
      stream << *it << "\n";
    }
  }
} // namespace HexMesher
