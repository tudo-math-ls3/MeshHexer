#pragma once

#include <hexmesher_config.hpp>

#ifndef HEXMESHER_PRAGMA_OMP
#ifdef HEXMESHER_HAVE_OMP
#define HEXMESHER_PRAGMA_OMP_HELPER(x) _Pragma(#x)
#define HEXMESHER_PRAGMA_OMP(x) HEXMESHER_PRAGMA_OMP_HELPER(omp x)
#else
#define HEXMESHER_PRAGMA_OMP(x)
#endif
#endif

