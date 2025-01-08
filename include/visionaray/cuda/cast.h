// This file is distributed under the MIT license.
// See the LICENSE file for details.

#pragma once

#ifndef VSNRAY_CUDA_CAST_H
#define VSNRAY_CUDA_CAST_H 1

#include <cuda_runtime.h>
#include <vector_functions.h>

#include "../gpu/cast.h"

namespace visionaray
{
namespace cuda
{

template <typename Dest, typename Source>
VSNRAY_FUNC
inline Dest cast(Source const& value)
{
    return ::visionaray::gpu::cast<Dest>(value);
}

} // cuda
} // visionaray

#endif // VSNRAY_CUDA_CAST_H
