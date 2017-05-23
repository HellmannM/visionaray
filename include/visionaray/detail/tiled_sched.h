// This file is distributed under the MIT license.
// See the LICENSE file for details.

#pragma once

#ifndef VSNRAY_DETAIL_TILED_SCHED_H
#define VSNRAY_DETAIL_TILED_SCHED_H 1

#include <memory>

namespace visionaray
{

template <typename R>
class tiled_sched
{
public:

    explicit tiled_sched(unsigned num_threads);
   ~tiled_sched();

    template <typename K, typename SP>
    void frame(K kernel, SP sched_params, unsigned frame_num = 0);

    void reset(unsigned num_threads);

private:

    struct impl;
    std::unique_ptr<impl> const impl_;

};

} // visionaray

#include "tiled_sched.inl"

#endif // VSNRAY_DETAIL_TILED_SCHED_H
