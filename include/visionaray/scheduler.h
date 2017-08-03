// This file is distributed under the MIT license.
// See the LICENSE file for details.

#pragma once

#ifndef VSNRAY_SCHEDULER_H
#define VSNRAY_SCHEDULER_H 1

#include <cstddef>
#include <type_traits>
#include <utility>

#include "detail/sched_common.h"
#include "pinhole_camera.h"

namespace visionaray
{

//-------------------------------------------------------------------------------------------------
// Base classes for scheduler params
//

template <typename Rect>
struct sched_params_base
{
    sched_params_base(Rect sb)
        : scissor_box(sb)
    {
    }

    Rect scissor_box;
};

template <typename Rect, typename Intersector>
struct sched_params_intersector_base : sched_params_base<Rect>
{
    using has_intersector = void;

    sched_params_intersector_base(Rect sb, Intersector& i)
        : sched_params_base<Rect>(sb)
        , intersector(i)
    {
    }

    Intersector& intersector;
};


//-------------------------------------------------------------------------------------------------
// Param structs for different pixel sampling strategies
//

template <typename... Args>
struct sched_params;

template <typename Base, typename RT, typename PxSamplerT>
struct sched_params<Base, RT, PxSamplerT> : Base
{
    using has_pinhole_camera = void;

    using rt_type               = RT;
    using pixel_sampler_type    = PxSamplerT;

    template <typename ...Args>
    sched_params(pinhole_camera const& c, RT& r, Args&&... args)
        : Base( std::forward<Args>(args)... )
        , cam(c)
        , rt(r)
    {
    }

    pinhole_camera cam;
    RT& rt;
};

template <typename Base, typename MT, typename RT, typename PxSamplerT>
struct sched_params<Base, MT, RT, PxSamplerT> : Base
{
    using has_camera_matrices = void;

    using rt_type               = RT;
    using pixel_sampler_type    = PxSamplerT;

    template <typename ...Args>
    sched_params(MT const& vm, MT const& pm, RT& r, Args&&... args)
        : Base( std::forward<Args>(args)... )
        , view_matrix(vm)
        , proj_matrix(pm)
        , rt(r)
    {
    }

    MT view_matrix;
    MT proj_matrix;
    RT& rt;
};


//-------------------------------------------------------------------------------------------------
// Deduce sched params type from members
//

namespace detail
{

template <typename SP>
class sched_params_has_cam
{
private:

    template <typename U>
    static std::true_type  test(typename U::has_camera*);

    template <typename U>
    static std::false_type test(...);

public:

    using type = decltype( test<typename std::decay<SP>::type>(nullptr) );

};

template <typename SP>
class sched_params_has_intersector
{
private:

    template <typename U>
    static std::true_type  test(typename U::has_intersector*);

    template <typename U>
    static std::false_type test(...);

public:

    using type = decltype( test<typename std::decay<SP>::type>(nullptr) );

};

template <typename SP>
class sched_params_has_view_matrix
{
private:

    template <typename U>
    static std::true_type  test(typename U::has_camera_matrices*);

    template <typename U>
    static std::false_type test(...);

public:

    using type = decltype( test<typename std::decay<SP>::type>(nullptr) );

};

} // detail


//-------------------------------------------------------------------------------------------------
// Sched params factory
//

template <typename PxSamplerT, typename RT>
auto make_sched_params(
        PxSamplerT              /* */,
        pinhole_camera const&   cam,
        RT&                     rt
        )
    -> sched_params<sched_params_base<recti>, RT, PxSamplerT>
{
    return sched_params<sched_params_base<recti>, RT, PxSamplerT>(
            cam,
            rt,
            recti(0, 0, rt.width(), rt.height())
            );
}

template <typename PxSamplerT, typename Intersector, typename RT>
auto make_sched_params(
        PxSamplerT              /* */,
        pinhole_camera const&   cam,
        RT&                     rt,
        Intersector&            isect
        )
    -> sched_params<sched_params_intersector_base<recti, Intersector>, RT, PxSamplerT>
{
    return sched_params<sched_params_intersector_base<recti, Intersector>, RT, PxSamplerT>{
            cam,
            rt,
            recti(0, 0, rt.width(), rt.height()),
            isect
            };
}


template <typename PxSamplerT, typename MT, typename RT>
auto make_sched_params(
        PxSamplerT      /* */,
        MT const&       view_matrix,
        MT const&       proj_matrix,
        RT&             rt
        )
    -> sched_params<sched_params_base<recti>, MT, RT, PxSamplerT>
{
    return sched_params<sched_params_base<recti>, MT, RT, PxSamplerT>(
            view_matrix,
            proj_matrix,
            rt,
            recti(0, 0, rt.width(), rt.height())
            );
}

template <typename PxSamplerT, typename MT, typename RT, typename Intersector>
auto make_sched_params(
        PxSamplerT      /* */,
        MT const&       view_matrix,
        MT const&       proj_matrix,
        RT&             rt,
        Intersector&    isect
        )
    -> sched_params<sched_params_intersector_base<recti, Intersector>, MT, RT, PxSamplerT>
{
    return sched_params<sched_params_intersector_base<recti, Intersector>, MT, RT, PxSamplerT>{
            view_matrix,
            proj_matrix,
            rt,
            recti(0, 0, rt.width(), rt.height()),
            isect
            };
}

template <
    typename First,
    typename = typename std::enable_if<!std::is_base_of<pixel_sampler::base_type, First>::value>::type,
    typename ...Args
    >
auto make_sched_params(First first, Args&&... args)
    -> decltype(make_sched_params(pixel_sampler::uniform_type{}, first, std::forward<Args>(args)...))
{
    return make_sched_params(pixel_sampler::uniform_type{}, first, std::forward<Args>(args)...);
}

} // visionaray

#ifdef __CUDACC__
#include "detail/cuda_sched.h"
#endif
#include "detail/simple_sched.h"
#if !defined(__MINGW32__) && !defined(__MINGW64__)
#include "detail/tiled_sched.h"
#endif

#endif // VSNRAY_SCHEDULER_H
