// This file is distributed under the MIT license.
// See the LICENSE file for details.

#pragma once

#ifndef VSNRAY_GET_SURFACE_H
#define VSNRAY_GET_SURFACE_H 1

#include <iterator>
#include <type_traits>
#include <utility>

#include "texture/texture.h"
#include "array.h"
#include "bvh.h"
#include "generic_primitive.h"
#include "get_color.h"
#include "get_normal.h"
#include "get_shading_normal.h"
#include "get_tex_coord.h"
#include "prim_traits.h"
#include "surface.h"
#include "tags.h"

namespace visionaray
{
namespace detail
{

//-------------------------------------------------------------------------------------------------
// Helper functions
//

// deduce simd surface type from params -------------------

template <typename Params, typename T, typename Enable = void>
struct simd_decl_surface;

template <typename Params, typename T>
struct simd_decl_surface<Params, T, typename std::enable_if<
        !(has_colors<Params>::value || has_textures<Params>::value)
        >::type>
{
private:

    enum { Size_ = simd::num_elements<T>::value };
    using N_     = typename Params::normal_type;
    using M_     = typename Params::material_type;

public:
    using type = surface<
        decltype(simd::pack(std::declval<array<N_, Size_>>())),
        decltype(simd::pack(std::declval<array<M_, Size_>>()))
        >;

    using array_type = array<surface<N_, M_>, Size_>;
};

template <typename Params, typename T>
struct simd_decl_surface<Params, T, typename std::enable_if<
        has_colors<Params>::value || has_textures<Params>::value
        >::type>
{
private:

    enum { Size_ = simd::num_elements<T>::value };
    using N_     = typename Params::normal_type;
    using M_     = typename Params::material_type;
    using C_     = typename Params::color_type;

public:
    using type = surface<
        decltype(simd::pack(std::declval<array<N_, Size_>>())),
        decltype(simd::pack(std::declval<array<M_, Size_>>())),
        decltype(simd::pack(std::declval<array<C_, Size_>>()))
        >;

    using array_type = array<surface<N_, M_, C_>, Size_>;
};


//-------------------------------------------------------------------------------------------------
// Struct containing both the geometric and the shading normal
//

template <typename V>
struct normal_pair
{
    V geometric_normal;
    V shading_normal;
};

template <typename V>
VSNRAY_FUNC
inline normal_pair<V> make_normal_pair(V const& gn, V const& sn)
{
    return normal_pair<V>{ gn, sn };
}


// TODO: consolidate the following with get_normal()?
// TODO: consolidate interface with get_color() and get_tex_coord()?

//-------------------------------------------------------------------------------------------------
// get_normal_pair()
//

template <typename Normals, typename HR, typename Primitive, typename NormalBinding>
VSNRAY_FUNC
inline auto get_normal_pair(
        Normals                     normals,
        HR const&                   hr,
        Primitive                   prim,
        NormalBinding               /* */,
        typename std::enable_if<num_normals<Primitive, NormalBinding>::value >= 2>::type* = 0
        )
    -> decltype( make_normal_pair(
            get_normal(normals, hr, prim, NormalBinding{}),
            get_shading_normal(normals, hr, prim, NormalBinding{})
            ) )
{
    return make_normal_pair(
            get_normal(normals, hr, prim, NormalBinding{}),
            get_shading_normal(normals, hr, prim, NormalBinding{})
        );
}

template <typename Normals, typename HR, typename Primitive, typename NormalBinding>
VSNRAY_FUNC
inline auto get_normal_pair(
        Normals                     normals,
        HR const&                   hr,
        Primitive                   /* */,
        NormalBinding               /* */,
        typename std::enable_if<num_normals<Primitive, NormalBinding>::value == 1>::type* = 0
        )
    -> decltype( make_normal_pair(
            get_normal(normals, hr, Primitive{}, NormalBinding{}),
            get_shading_normal(normals, hr, Primitive{}, NormalBinding{})
            ) )
{
    return make_normal_pair(
            get_normal(normals, hr, Primitive{}, NormalBinding{}),
            get_shading_normal(normals, hr, Primitive{}, NormalBinding{})
            );
}

template <typename Normals, typename HR, typename Primitive, typename NormalBinding>
VSNRAY_FUNC
inline auto get_normal_pair(
        Normals                     normals,
        HR const&                   hr,
        Primitive                   prim,
        NormalBinding               /* */,
        typename std::enable_if<num_normals<Primitive, NormalBinding>::value == 0>::type* = 0
        )
    -> decltype( make_normal_pair(
            get_normal(hr, prim),
            get_shading_normal(hr, prim)
            ) )
{
    VSNRAY_UNUSED(normals);

    return make_normal_pair(
            get_normal(hr, prim),
            get_shading_normal(hr, prim)
            );
}

template <typename HR, typename Primitive>
VSNRAY_FUNC
inline auto get_normal_pair(
        HR const&                   hr,
        Primitive                   prim
        )
    -> decltype( make_normal_pair(
            get_normal(hr, prim),
            get_shading_normal(hr, prim)
            ) )
{
    return make_normal_pair(
            get_normal(hr, prim),
            get_shading_normal(hr, prim)
            );
}


// get_normal_pair as functor for template arguments
struct get_normal_pair_t
{
    template <typename Normals, typename HR, typename Primitive, typename NormalBinding>
    VSNRAY_FUNC
    inline auto operator()(
            Normals                     normals,
            HR const&                   hr,
            Primitive                   prim,
            NormalBinding               /* */
            ) const
        -> decltype( get_normal_pair(normals, hr, prim, NormalBinding{}) )
    {
        return get_normal_pair(normals, hr, prim, NormalBinding{});
    }

    template <typename HR, typename Primitive>
    VSNRAY_FUNC
    inline auto operator()(
            HR const&                   hr,
            Primitive                   prim
            ) const
        -> decltype( get_normal_pair(hr, prim) )
    {
        return get_normal_pair(hr, prim);
    }
};


// overload for generic_primitive
template <
    typename Normals,
    typename HR,
    typename ...Ts,
    typename NormalBinding
    >
VSNRAY_FUNC
inline auto get_normal_pair(
        Normals                     normals,
        HR const&                   hr,
        generic_primitive<Ts...>    prim,
        NormalBinding               /* */
        )
    -> normal_pair<typename std::iterator_traits<Normals>::value_type>
{
    get_normal_from_generic_primitive_visitor<
        get_normal_pair_t,
        normal_pair<typename std::iterator_traits<Normals>::value_type>,
        NormalBinding,
        Normals,
        HR
        >visitor(
            normals,
            hr
            );

    return apply_visitor( visitor, prim );
}


//-------------------------------------------------------------------------------------------------
// dispatch function for get_normal()
//

template <
    typename Params,
    typename Normals,
    typename HR,
    typename Primitive = typename Params::primitive_type,
    typename NormalBinding = typename Params::normal_binding,
    typename = typename std::enable_if<!is_any_bvh<Primitive>::value>::type
    >
VSNRAY_FUNC
inline auto get_normal_dispatch(
        Params const&   params,
        Normals         normals,
        HR const&       hr,
        typename std::enable_if<
            num_normals<Primitive, NormalBinding>::value == 1
            >::type* = 0
        )
    -> decltype( get_normal_pair(normals, hr, Primitive{}, NormalBinding{}) )
{
    VSNRAY_UNUSED(params);

    return get_normal_pair(normals, hr, Primitive{}, NormalBinding{});
}

template <
    typename Params,
    typename Normals,
    typename HR,
    typename Primitive = typename Params::primitive_type,
    typename NormalBinding = typename Params::normal_binding,
    typename = typename std::enable_if<!is_any_bvh<Primitive>::value>::type
    >
VSNRAY_FUNC
inline auto get_normal_dispatch(
        Params const&   params,
        Normals         normals,
        HR const&       hr,
        typename std::enable_if<
            num_normals<Primitive, NormalBinding>::value != 1
            >::type* = 0
        )
    -> decltype( get_normal_pair(normals, hr, params.prims.begin[hr.prim_id], NormalBinding{}) )
{
    return get_normal_pair(normals, hr, params.prims.begin[hr.prim_id], NormalBinding{});
}

// overload for BVHs
template <
    typename Params,
    typename Normals,
    typename R,
    typename Base,
    typename Primitive = typename Params::primitive_type,
    typename NormalBinding = typename Params::normal_binding,
    typename = typename std::enable_if<is_any_bvh<Primitive>::value>::type
    >
VSNRAY_FUNC
inline auto get_normal_dispatch(
        Params const&                  params,
        Normals                        normals,
        hit_record_bvh<R, Base> const& hr,
        typename std::enable_if<
            num_normals<typename Primitive::primitive_type, NormalBinding>::value == 1
            >::type* = 0
        )
    -> decltype( get_normal_pair(
                normals,
                static_cast<Base const&>(hr),
                typename Primitive::primitive_type{},
                NormalBinding{}
            ) )
{
    VSNRAY_UNUSED(params);

    return get_normal_pair(
                normals,
                static_cast<Base const&>(hr),
                typename Primitive::primitive_type{},
                NormalBinding{}
            );
}

template <
    typename Params,
    typename Normals,
    typename R,
    typename Base,
    typename Primitive = typename Params::primitive_type,
    typename NormalBinding = typename Params::normal_binding,
    typename = typename std::enable_if<is_any_bvh<Primitive>::value>::type
    >
VSNRAY_FUNC
inline auto get_normal_dispatch(
        Params const&                  params,
        Normals                        normals,
        hit_record_bvh<R, Base> const& hr,
        typename std::enable_if<
            num_normals<typename Primitive::primitive_type, NormalBinding>::value != 1
            >::type* = 0
        )
    -> decltype( get_normal_pair(
            normals,
            static_cast<Base const&>(hr),
            typename Primitive::primitive_type{},
            NormalBinding{}
            ) )
{
    // Find the BVH that contains prim_id
    size_t num_primitives_total = 0;

    size_t i = 0;
    while (static_cast<size_t>(hr.prim_id) >= num_primitives_total + params.prims.begin[i].num_primitives())
    {
        num_primitives_total += params.prims.begin[i++].num_primitives();
    }

    return get_normal_pair(
            normals,
            static_cast<Base const&>(hr),
            params.prims.begin[i].primitive(hr.primitive_list_index),
            typename Params::normal_binding{}
            );
}


//-------------------------------------------------------------------------------------------------
// Sample textures with range check
//

template <typename HR, typename Params>
VSNRAY_FUNC
inline typename Params::color_type get_tex_color(
        HR const&                      hr,
        Params const&                  params,
        std::integral_constant<int, 1> /* */
        )
{
    using P = typename Params::primitive_type;
    using C = typename Params::color_type;

    auto coord = get_tex_coord(params.tex_coords, hr, P{});

    auto const& tex = params.textures[hr.geom_id];
    return tex.width() > 0 ? C(visionaray::tex1D(tex, coord)) : C(1.0);
}

template <typename HR, typename Params>
VSNRAY_FUNC
inline typename Params::color_type get_tex_color(
        HR const&                      hr,
        Params const&                  params,
        std::integral_constant<int, 2> /* */
        )
{
    using P = typename Params::primitive_type;
    using C = typename Params::color_type;

    auto coord = get_tex_coord(params.tex_coords, hr, P{});

    auto const& tex = params.textures[hr.geom_id];
    return tex.width() > 0 && tex.height() > 0
            ? C(visionaray::tex2D(tex, coord))
            : C(1.0);
}

template <typename HR, typename Params>
VSNRAY_FUNC
inline typename Params::color_type get_tex_color(
        HR const&                      hr,
        Params const&                  params,
        std::integral_constant<int, 3> /* */
        )
{
    using P = typename Params::primitive_type;
    using C = typename Params::color_type;

    auto coord = get_tex_coord(params.tex_coords, hr, P{});

    auto const& tex = params.textures[hr.geom_id];
    return tex.width() > 0 && tex.height() > 0 && tex.depth() > 0
            ? C(visionaray::tex3D(tex, coord))
            : C(1.0);
}


//-------------------------------------------------------------------------------------------------
//
//

template <typename HR, typename Params>
VSNRAY_FUNC
inline auto get_surface_impl(
        has_no_normals_tag  /* */,
        has_no_colors_tag   /* */,
        has_no_textures_tag /* */,
        HR const&           hr,
        Params const&       params
        )
    -> surface<typename Params::normal_type, typename Params::material_type>
{
    auto ns = get_normal_dispatch(params, nullptr, hr);
    return make_surface(
            ns.geometric_normal,
            ns.shading_normal,
            params.materials[hr.geom_id]
            );
}

template <typename HR, typename Params>
VSNRAY_FUNC
inline auto get_surface_impl(
        has_normals_tag     /* */,
        has_no_colors_tag   /* */,
        has_no_textures_tag /* */,
        HR const&           hr,
        Params const&       params
        )
    -> surface<typename Params::normal_type, typename Params::material_type>
{
    auto ns = get_normal_dispatch(params, params.normals, hr);
    return make_surface(
            ns.geometric_normal,
            ns.shading_normal,
            params.materials[hr.geom_id]
            );
}

template <typename HR, typename Params>
VSNRAY_FUNC
inline auto get_surface_impl(
        has_normals_tag     /* */,
        has_no_colors_tag   /* */,
        has_textures_tag    /* */,
        HR const&           hr,
        Params const&       params
        )
    -> surface<
            typename Params::normal_type,
            typename Params::material_type,
            typename Params::color_type
            >
{
    auto ns = get_normal_dispatch(params, params.normals, hr);
    auto tc = get_tex_color(
                    hr,
                    params,
                    std::integral_constant<int, Params::texture_type::dimensions>{}
                    );

    return make_surface(
            ns.geometric_normal,
            ns.shading_normal,
            params.materials[hr.geom_id],
            tc
            );
}

template <typename HR, typename Params>
VSNRAY_FUNC
inline auto get_surface_impl(
        has_no_normals_tag  /* */,
        has_colors_tag      /* */,
        has_textures_tag    /* */,
        HR const&           hr,
        Params const&       params
        )
    -> surface<
            typename Params::normal_type,
            typename Params::material_type,
            typename Params::color_type
            >
{
    using P = typename Params::primitive_type;

    auto ns    = get_normal_dispatch(params, nullptr, hr);
    auto color = get_color(params.colors, hr, P{}, typename Params::color_binding{});
    auto tc    = get_tex_color(
                        hr,
                        params,
                        std::integral_constant<int, Params::texture_type::dimensions>{}
                        );

    return make_surface(
            ns.geometric_normal,
            ns.shading_normal,
            params.materials[hr.geom_id],
            color * tc
            );
}

template <typename HR, typename Params>
VSNRAY_FUNC
inline auto get_surface_impl(
        has_normals_tag     /* */,
        has_colors_tag      /* */,
        has_textures_tag    /* */,
        HR const&           hr,
        Params const&       params
        )
    -> surface<
            typename Params::normal_type,
            typename Params::material_type,
            typename Params::color_type
            >
{
    using P = typename Params::primitive_type;

    auto ns    = get_normal_dispatch(params, params.normals, hr);
    auto color = get_color(params.colors, hr, P{}, typename Params::color_binding{});
    auto tc    = get_tex_color(
                        hr,
                        params,
                        std::integral_constant<int, Params::texture_type::dimensions>{}
                        );

    return make_surface(
            ns.geometric_normal,
            ns.shading_normal,
            params.materials[hr.geom_id],
            color * tc
            );
}


//-------------------------------------------------------------------------------------------------
// SIMD
//

template <
    typename NormalsTag,
    typename ColorsTag,
    typename TexturesTag,
    typename HR,
    typename Params,
    typename = typename std::enable_if<simd::is_simd_vector<typename HR::scalar_type>::value>::type
    >
VSNRAY_FUNC
inline auto get_surface_impl(
        NormalsTag    /* */,
        ColorsTag     /* */,
        TexturesTag   /* */,
        HR const&     hr,
        Params const& params
        )
    -> typename simd_decl_surface<Params, typename HR::scalar_type>::type
{
    using T = typename HR::scalar_type;

    auto hrs = unpack(hr);

    typename simd_decl_surface<Params, T>::array_type surfs;

    for (int i = 0; i < simd::num_elements<T>::value; ++i)
    {
        if (hrs[i].hit)
        {
            surfs[i] = get_surface_impl(
                    NormalsTag{},
                    ColorsTag{},
                    TexturesTag{},
                    hrs[i],
                    params
                    );
        }
    }

    return simd::pack(surfs);
}

} // detail


template <typename HR, typename Params>
VSNRAY_FUNC
inline auto get_surface(HR const& hr, Params const& p)
    -> decltype( detail::get_surface_impl(
            detail::has_normals<Params>{},
            detail::has_colors<Params>{},
            detail::has_textures<Params>{},
            hr,
            p
            ) )
{
    return detail::get_surface_impl(
            detail::has_normals<Params>{},
            detail::has_colors<Params>{},
            detail::has_textures<Params>{},
            hr,
            p
            );
}

} // visionaray

#endif // VSNRAY_SURFACE_H
