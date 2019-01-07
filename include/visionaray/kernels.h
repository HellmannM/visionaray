// This file is distributed under the MIT license.
// See the LICENSE file for details.

#pragma once

#ifndef VSNRAY_KERNELS_H
#define VSNRAY_KERNELS_H 1

#include <iterator>
#include <limits>

#include "math/forward.h"
#include "math/vector.h"
#include "prim_traits.h"
#include "tags.h"

namespace visionaray
{

//-------------------------------------------------------------------------------------------------
// Parameter structs for built-in kernels
//
// Use the make_kernel_params() factory function to create:
//
// Parameters:
//      NormalBindingTag:   [normals_per_face_binding|normals_per_vertex_binding]
//      ColorBindingTag:    [colors_per_face_binding|colors_per_vertex_binding] (TODO: per geometry)
//      primitives [..), normals, materials, lights [..)
//      num_bounces:        depth of the ray tree (applies only to recursive algorithms)
//      scene_epsilon:      used as an offset to bias ray origins to avoid self intersections
//      background_color:   RGBA
//      ambient_color:      RGBA
//
// w/ normals:
//
//  make_kernel_params(
//      NormalBindingTag,
//      primitives_begin,
//      primitives_end,
//      normals,
//      materials,
//      lights_begin,
//      lights_end,
//      num_bounces,
//      scene_epsilon,
//      background_color,
//      ambient_color
//      );
//
//
// w/ normals and textures:
//
//  make_kernel_params(
//      NormalBindingTag,
//      primitives_begin,
//      primitives_end,
//      normals,
//      texture_coordinates,
//      materials,
//      textures,
//      lights_begin,
//      lights_end,
//      num_bounces,
//      scene_epsilon,
//      background_color,
//      ambient_color
//      );
//
// w/ textures, normals and colors:
//
//  make_kernel_params(
//      NormalBindingTag,
//      ColorBindingTag,
//      primitives_begin,
//      primitives_end,
//      normals,
//      texture_coordinates,
//      materials,
//      colors,
//      textures,
//      lights_begin,
//      lights_end,
//      num_bounces,
//      scene_epsilon,
//      background_color,
//      ambient_color
//      );
//
//
//-------------------------------------------------------------------------------------------------

//-------------------------------------------------------------------------------------------------
// Param structs
//

template <typename ...Args>
struct kernel_params;

template <
    typename NormalBinding,
    typename Primitives,
    typename Normals,
    typename Materials,
    typename Lights,
    typename Color,
    typename ...Args
    >
struct kernel_params<
        NormalBinding,
        Primitives,
        Normals,
        Materials,
        Lights,
        Color,
        Args...
        >
{
    using normal_binding    = NormalBinding;
    using primitive_type    = typename std::iterator_traits<Primitives>::value_type;
    using normal_type       = typename std::iterator_traits<Normals>::value_type;
    using material_type     = typename std::iterator_traits<Materials>::value_type;
    using light_type        = typename std::iterator_traits<Lights>::value_type;
    using color_type        = vector<3, typename scalar_type<primitive_type>::type>;

    using color_binding     = unspecified_binding;

    struct
    {
        Primitives begin;
        Primitives end;
    } prims;

    Normals   geometric_normals;
    Normals   shading_normals;
    Materials materials;

    struct
    {
        Lights begin;
        Lights end;
    } lights;

    unsigned num_bounces;
    float epsilon;

    Color bg_color;
    Color ambient_color;
};

template <
    typename NormalBinding,
    typename Primitives,
    typename Normals,
    typename TexCoords,
    typename Materials,
    typename Textures,
    typename Lights,
    typename Color,
    typename ...Args
    >
struct kernel_params<
        NormalBinding,
        Primitives,
        Normals,
        TexCoords,
        Materials,
        Textures,
        Lights,
        Color,
        Args...
        >
{
    using has_textures      = void;

    using normal_binding    = NormalBinding;
    using primitive_type    = typename std::iterator_traits<Primitives>::value_type;
    using normal_type       = typename std::iterator_traits<Normals>::value_type;
    using tex_coords_type   = typename std::iterator_traits<TexCoords>::value_type;
    using material_type     = typename std::iterator_traits<Materials>::value_type;
    using texture_type      = typename std::iterator_traits<Textures>::value_type;
    using light_type        = typename std::iterator_traits<Lights>::value_type;
    using color_type        = vector<3, typename scalar_type<primitive_type>::type>;

    using color_binding     = unspecified_binding;

    struct
    {
        Primitives begin;
        Primitives end;
    } prims;

    Normals   geometric_normals;
    Normals   shading_normals;
    TexCoords tex_coords;
    Materials materials;
    Textures  textures;

    struct
    {
        Lights begin;
        Lights end;
    } lights;

    unsigned num_bounces;
    float epsilon;

    Color bg_color;
    Color ambient_color;
};

template <
    typename NormalBinding,
    typename ColorBinding,
    typename Primitives,
    typename Normals,
    typename TexCoords,
    typename Materials,
    typename Colors,
    typename Textures,
    typename Lights,
    typename Color,
    typename ...Args
    >
struct kernel_params<
        NormalBinding,
        ColorBinding,
        Primitives,
        Normals,
        TexCoords,
        Materials,
        Colors,
        Textures,
        Lights,
        Color,
        Args...
        >
{
    using has_textures      = void;
    using has_colors        = void;

    using normal_binding    = NormalBinding;
    using color_binding     = ColorBinding;
    using primitive_type    = typename std::iterator_traits<Primitives>::value_type;
    using normal_type       = typename std::iterator_traits<Normals>::value_type;
    using tex_coords_type   = typename std::iterator_traits<TexCoords>::value_type;
    using material_type     = typename std::iterator_traits<Materials>::value_type;
    using color_type        = typename std::iterator_traits<Colors>::value_type;
    using texture_type      = typename std::iterator_traits<Textures>::value_type;
    using light_type        = typename std::iterator_traits<Lights>::value_type;

    struct
    {
        Primitives begin;
        Primitives end;
    } prims;

    Normals   geometric_normals;
    Normals   shading_normals;
    TexCoords tex_coords;
    Materials materials;
    Colors    colors;
    Textures  textures;

    struct
    {
        Lights begin;
        Lights end;
    } lights;

    unsigned num_bounces;
    float epsilon;

    Color bg_color;
    Color ambient_color;
};


//-------------------------------------------------------------------------------------------------
// Factory for param structs
//

// default ------------------------------------------------

template <
    typename Primitives,
    typename Materials,
    typename Lights
    >
auto make_kernel_params(
        Primitives const&   begin,
        Primitives const&   end,
        Materials const&    materials,
        Lights const&       lbegin,
        Lights const&       lend,
        unsigned            num_bounces     = 5,
        float               epsilon         = std::numeric_limits<float>::epsilon(),
        vec4 const&         bg_color        = vec4(0.0),
        vec4 const&         ambient_color   = vec4(0.0)
        )
    -> kernel_params<
            normals_per_vertex_binding,
            Primitives,
            vector<3, typename scalar_type<typename std::iterator_traits<Primitives>::value_type>::type>*,
            Materials,
            Lights,
            vec4
            >
{
    return {
        { begin, end },
        nullptr,
        nullptr,
        materials,
        { lbegin, lend },
        num_bounces,
        epsilon,
        bg_color,
        ambient_color
        };
}


// w/ normals ---------------------------------------------

template <
    typename NormalBinding,
    typename Primitives,
    typename Normals,
    typename Materials,
    typename Lights,
    typename = typename std::enable_if<is_normal_binding<NormalBinding>::value>::type
    >
auto make_kernel_params(
        NormalBinding       /* */,
        Primitives const&   begin,
        Primitives const&   end,
        Normals const&      geometric_normals,
        Normals const&      shading_normals,
        Materials const&    materials,
        Lights const&       lbegin,
        Lights const&       lend,
        unsigned            num_bounces     = 5,
        float               epsilon         = std::numeric_limits<float>::epsilon(),
        vec4 const&         bg_color        = vec4(0.0),
        vec4 const&         ambient_color   = vec4(0.0)
        )
    -> kernel_params<NormalBinding, Primitives, Normals, Materials, Lights, vec4>
{
    return {
        { begin, end },
        geometric_normals,
        shading_normals,
        materials,
        { lbegin, lend },
        num_bounces,
        epsilon,
        bg_color,
        ambient_color
        };
}


// w/ normals and textures --------------------------------

template <
    typename NormalBinding,
    typename Primitives,
    typename Normals,
    typename TexCoords,
    typename Materials,
    typename Textures,
    typename Lights,
    typename = typename std::enable_if<is_normal_binding<NormalBinding>::value>::type
    >
auto make_kernel_params(
        NormalBinding       /* */,
        Primitives const&   begin,
        Primitives const&   end,
        Normals const&      geometric_normals,
        Normals const&      shading_normals,
        TexCoords const&    tex_coords,
        Materials const&    materials,
        Textures const&     textures,
        Lights const&       lbegin,
        Lights const&       lend,
        unsigned            num_bounces     = 5,
        float               epsilon         = std::numeric_limits<float>::epsilon(),
        vec4 const&         bg_color        = vec4(0.0),
        vec4 const&         ambient_color   = vec4(1.0)
        )
    -> kernel_params<NormalBinding, Primitives, Normals, TexCoords, Materials, Textures, Lights, vec4>
{
    return {
        { begin, end },
        geometric_normals,
        shading_normals,
        tex_coords,
        materials,
        textures,
        { lbegin, lend },
        num_bounces,
        epsilon,
        bg_color,
        ambient_color
        };
}


// w/ normals, colors and textures ------------------------

template <
    typename NormalBinding,
    typename ColorBinding,
    typename Primitives,
    typename Normals,
    typename TexCoords,
    typename Materials,
    typename Colors,
    typename Textures,
    typename Lights,
    typename = typename std::enable_if<is_normal_binding<NormalBinding>::value>::type,
    typename = typename std::enable_if<is_color_binding<ColorBinding>::value>::type
    >
auto make_kernel_params(
        NormalBinding       /* */,
        ColorBinding        /* */,
        Primitives const&   begin,
        Primitives const&   end,
        Normals const&      geometric_normals,
        Normals const&      shading_normals,
        TexCoords const&    tex_coords,
        Materials const&    materials,
        Colors const&       colors,
        Textures const&     textures,
        Lights const&       lbegin,
        Lights const&       lend,
        unsigned            num_bounces     = 5,
        float               epsilon         = std::numeric_limits<float>::epsilon(),
        vec4 const&         bg_color        = vec4(0.0),
        vec4 const&         ambient_color   = vec4(1.0)
        )
    -> kernel_params<
        NormalBinding,
        ColorBinding,
        Primitives,
        Normals,
        TexCoords,
        Materials,
        Colors,
        Textures,
        Lights,
        vec4
        >
{
    return {
        { begin, end },
        geometric_normals,
        shading_normals,
        tex_coords,
        materials,
        colors,
        textures,
        { lbegin, lend },
        num_bounces,
        epsilon,
        bg_color,
        ambient_color
        };
}

} // visionaray

#include "detail/pathtracing.inl"
#include "detail/simple.inl"
#include "detail/whitted.inl"

#endif // VSNRAY_KERNELS_H
