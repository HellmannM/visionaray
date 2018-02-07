// This file is distributed under the MIT license.
// See the LICENSE file for details.

#pragma once

#ifndef VSNRAY_MATH_INTERSECT_H
#define VSNRAY_MATH_INTERSECT_H 1

#include <type_traits>

#include <visionaray/array.h>

#include "simd/type_traits.h"

#include "aabb.h"
#include "limits.h"
#include "plane.h"
#include "ray.h"
#include "sphere.h"
#include "triangle.h"
#include "vector.h"

//#define HLS_NO_XIL_FPO_LIB
//#include <ap_int.h>
#include "fixed.h"

namespace MATH_NAMESPACE
{

template <typename T1, typename T2>
struct hit_record;


//-------------------------------------------------------------------------------------------------
// ray / aabb
//

template <typename T, typename U>
struct hit_record<basic_ray<T>, basic_aabb<U>>
{
    using scalar_type = T;
    using mask_type = simd::mask_type_t<T>;

    MATH_FUNC hit_record()
        : hit(false)
        , tnear(numeric_limits<float>::max())
        , tfar(-numeric_limits<float>::max())
    {
    }

    mask_type   hit;
    scalar_type tnear;
    scalar_type tfar;

};

//template <typename T, typename U>
//MATH_FUNC
//inline hit_record<basic_ray<T>, basic_aabb<U>> intersect(
//        basic_ray<T> const&     ray,
//        basic_aabb<U> const&    aabb,
//        vector<3, T>            inv_dir
//        )
//{
//#pragma HLS ALLOCATION instances=mul limit=2 operation
//	hit_record<basic_ray<T>, basic_aabb<U>> result;
//
//	typedef fixed<16,16> F; // TODO find optimal values
//	vector<3, F> ori(F(ray.ori.x), F(ray.ori.y), F(ray.ori.z));
//	vector<3, F> inv(F(inv_dir.x), F(inv_dir.y), F(inv_dir.z));
//    vector<3, F> t1 = (vector<3, F>(aabb.min) - ori) * inv;
//    vector<3, F> t2 = (vector<3, F>(aabb.max) - ori) * inv;
////    vector<3, T> t1 = (vector<3, T>(aabb.min) - ray.ori) * inv_dir;
////    vector<3, T> t2 = (vector<3, T>(aabb.max) - ray.ori) * inv_dir;
//
//    result.tnear = min_max( t1.x, t2.x, min_max(t1.y, t2.y, min(t1.z, t2.z)) );
//    result.tfar  = max_min( t1.x, t2.x, max_min(t1.y, t2.y, max(t1.z, t2.z)) );
//    result.hit   = result.tfar >= result.tnear;
//
//    return result;
//}
//
//template <typename T, typename U>
//MATH_FUNC
//inline hit_record<basic_ray<T>, basic_aabb<U>> intersect(
//        basic_ray<T> const& ray,
//        basic_aabb<U> const& aabb
//        )
//{
////    vector<3, T> inv_dir = T(1.0) / ray.dir;
//    vector<3, T> inv_dir = ray.inv_dir;
//    return intersect(ray, aabb, inv_dir);
//}


//-------------------------------------------------------------------------------------------------
// general ray primitive hit record
//

template <typename T>
struct hit_record<basic_ray<T>, primitive<unsigned>>
{
    using scalar_type = T;
    using int_type = simd::int_type_t<T>;
    using mask_type = simd::mask_type_t<T>;

    MATH_FUNC hit_record()
        : hit(false)
        , prim_id(0)
        , geom_id(0)
        , t(numeric_limits<float>::max())
        , u(0.0)
        , v(0.0)
    {
    }


    mask_type hit;
    int_type prim_id;
    int_type geom_id;

    T t;
    vector<3, T> isect_pos;

    T u;
    T v;
};


////-------------------------------------------------------------------------------------------------
//// ray / triangle
////
//template <typename T, typename U>
//MATH_FUNC
//inline hit_record<basic_ray<T>, primitive<unsigned>> intersect(
//        basic_ray<T> const&                     ray,
//        hanika_triangle<3, U, unsigned> const&   tri
//        )
//{
//#pragma HLS ALLOCATION instances=mul limit=2 operation
//
///*
//    typedef vector<3, T> vec_type;
//
//    hit_record<basic_ray<T>, primitive<unsigned>> result;
//    result.t = T(-1.0);
//
//    // case T != U
//    vec_type v1(tri.v1);
//    vec_type e1(tri.e1);
//    vec_type e2(tri.e2);
//
//    vec_type s1 = cross(ray.dir, e2);
//    T div = dot(s1, e1);
//
//    result.hit = ( div != T(0.0) );
//
//    if ( !any(result.hit) )
//    {
//        return result;
//    }
//
////    T inv_div = T(1.0) / div;
//
//    vec_type d = ray.ori - v1;
////    T b1 = dot(d, s1) * inv_div;
//    T b1 = dot(d, s1) / div;
//
//    result.hit &= ( b1 >= T(0.0) && b1 <= T(1.0) );
//
//    if ( !any(result.hit) )
//    {
//        return result;
//    }
//
//    vec_type s2 = cross(d, e1);
////    T b2 = dot(ray.dir, s2) * inv_div;
//    T b2 = dot(ray.dir, s2) / div;
//
//    result.hit &= ( b2 >= T(0.0) && b1 + b2 <= T(1.0) );
//
//    if ( !any(result.hit) )
//    {
//        return result;
//    }
//
//    result.prim_id = tri.prim_id;
//    result.geom_id = tri.geom_id;
////    result.t = dot(e2, s2) * inv_div;
//    result.t = dot(e2, s2) / div;
//    result.u = b1;
//    result.v = b2;
//    return result;
//*/
//
//// hanika
//    hit_record<basic_ray<T>, primitive<unsigned>> result;
//    result.t = T(-1.0);
//    result.hit = false;
//
//    int i0 = (int)tri.i0;
//
//    typedef fixed<16, 16> F; // TODO: find optimal values
//
//// kann nicht mit [] auf vector<ap_fixed> zugreifen. warum? -> TODO
//    F ori0, ori1, ori2, dir0, dir1, dir2;
//    if(i0 == 0){
//        ori0 = (F)ray.ori.x;
//        ori1 = (F)ray.ori.y;
//        ori2 = (F)ray.ori.z;
//        dir0 = (F)ray.dir.x;
//        dir1 = (F)ray.dir.y;
//        dir2 = (F)ray.dir.z;
//    }else if(i0 == 1){
//        ori0 = (F)ray.ori.y;
//        ori1 = (F)ray.ori.z;
//        ori2 = (F)ray.ori.x;
//        dir0 = (F)ray.dir.y;
//        dir1 = (F)ray.dir.z;
//        dir2 = (F)ray.dir.x;
//    }else{
//        ori0 = (F)ray.ori.z;
//        ori1 = (F)ray.ori.x;
//        ori2 = (F)ray.ori.y;
//        dir0 = (F)ray.dir.z;
//        dir1 = (F)ray.dir.x;
//        dir2 = (F)ray.dir.y;
//	}
//
//    F denominator(dir0 + dir1 * (F)tri.np + dir2 * (F)tri.nq);
//    if (denominator == F(0.0)){
////    	std::cout << "div by 0 in prim intersect!" << std::endl; //TODO: why does this happen
//    	return result;
//    }
//
//    F t = ((F)tri.d - ori0 - ori1 * (F)tri.np - ori2 * (F)tri.nq)
//            / denominator;
////    F t = ((F)tri.d - (F)ray.ori[tri.i0] - (F)ray.ori[i1] * (F)tri.np - (F)ray.ori[i2] * (F)tri.nq)
////            / ((F)ray.dir[tri.i0] + (F)ray.dir[i1] * (F)tri.np + (F)ray.dir[i2] * (F)tri.nq);
//
//    if (t < (F)0)
//        return result;
//
//    F kp = ori1 + t * dir1 - (F)tri.pp;
//    F kq = ori2 + t * dir2 - (F)tri.pq;
//    F u = (F)tri.e1p * kq - (F)tri.e1q * kp;
//    F v = (F)tri.e2q * kp - (F)tri.e2p * kq;
////    F kp = (F)ray.ori[i1] + t * (F)ray.dir[i1] - (F)tri.pp;
////    F kq = (F)ray.ori[i2] + t * (F)ray.dir[i2] - (F)tri.pq;
////    F u = (F)tri.e[0] * kq - (F)tri.e[2] * kp;
////    F v = (F)tri.e[3] * kp - (F)tri.e[1] * kq;
//
//    result.hit = (u >= (F)-0.02 && v >= (F)-0.02 && (u + v) <= (F)1.04);
//    if (result.hit)
//        result.t = (T)t;
//
//    return result;
//
//}


//-------------------------------------------------------------------------------------------------
// ray / sphere
//

template <typename T, typename U>
MATH_FUNC
inline hit_record<basic_ray<T>, primitive<unsigned>> intersect(
        basic_ray<T> const&                 ray,
        basic_sphere<U, unsigned> const&    sphere
        )
{

    typedef basic_ray<T> ray_type;
    typedef vector<3, T> vec_type;

    ray_type r = ray;
    r.ori -= vec_type( sphere.center );

    auto A = dot(r.dir, r.dir);
    auto B = dot(r.dir, r.ori) * T(2.0);
    auto C = dot(r.ori, r.ori) - sphere.radius * sphere.radius;

    // solve Ax**2 + Bx + C
    auto disc = B * B - T(4.0) * A * C;
    auto valid = disc >= T(0.0);

    auto root_disc = select(valid, sqrt(disc), disc);

    auto q = select( B < T(0.0), T(-0.5) * (B - root_disc), T(-0.5) * (B + root_disc) );

    auto t1 = q / A;
    auto t2 = C / q;

    hit_record<basic_ray<T>, primitive<unsigned>> result;
    result.hit = valid;
    result.prim_id = sphere.prim_id;
    result.geom_id = sphere.geom_id;
    result.t = select( valid, select( t1 > t2, t2, t1 ), T(-1.0) );
    return result;
}


//-------------------------------------------------------------------------------------------------
// ray / plane
//

template <typename T, typename U>
MATH_FUNC
inline hit_record<basic_ray<T>, primitive<unsigned>> intersect(
        basic_ray<T> const&                 ray,
        basic_plane<3, U, unsigned> const&  p
        )
{
    hit_record<basic_ray<T>, primitive<unsigned>> result;

    vector<3, T> normal(p.normal);
    T offset(p.offset);

    T s = dot(normal, ray.dir);

    result.hit = (s != T(0.0));

    result.prim_id = p.prim_id;
    result.geom_id = p.geom_id;

    result.t = select(
            result.hit,
            ( offset - dot(normal, ray.ori) ) / s,
            T(-1.0)
            );

    result.isect_pos = select(
            result.hit,
            ray.ori + result.t * ray.dir,
            result.isect_pos
            ); // TODO: sphere/triangle don't assign isect_pos in intersect

    return result;
}


//-------------------------------------------------------------------------------------------------
// plane / sphere
//

template <typename T, typename P>
struct hit_record<basic_plane<3, T>, basic_sphere<T, P>>
{
    using scalar_type = T;
    using mask_type = simd::mask_type_t<T>;

    // TODO: make this part of the API?
    struct basic_circle
    {
        vector<3, T> center;
        T radius;
    };

    MATH_FUNC hit_record()
        : hit(false)
    {
    }

    mask_type       hit;
    basic_circle    circle;
};

template <typename T, typename P>
MATH_FUNC
inline hit_record<basic_plane<3, T>, basic_sphere<T, P>> intersect(
        basic_plane<3, T> const&    plane,
        basic_sphere<T, P> const&   sphere
        )
{
    hit_record<basic_plane<3, T>, basic_sphere<T, P>> result;

    T dist = plane.offset - dot(sphere.center, plane.normal);

    auto hit = abs(dist) <= sphere.radius;

    result.hit           = hit;
    result.circle.center = select(hit, sphere.center + plane.normal * dist, vector<3, T>());
    result.circle.radius = select(hit, sqrt(sphere.radius * sphere.radius - dist * dist), T(0.0));

    return result;
}

template <typename T, typename P>
MATH_FUNC
inline hit_record<basic_sphere<T, P>, basic_plane<3, T>> intersect(
        basic_sphere<T, P> const&   sphere,
        basic_plane<3, T> const&    plane
        )
{
    return intersect(plane, sphere);
}


//-------------------------------------------------------------------------------------------------
// pack / unpack functions for hit records
//

namespace simd
{

// general primitive --------------------------------------

template <
    size_t N,
    typename T = float_from_simd_width_t<N>
    >
MATH_FUNC
inline hit_record<basic_ray<T>, primitive<unsigned>> pack(
        array<hit_record<ray, primitive<unsigned>>, N> const& hrs
        )
{
    hit_record<basic_ray<T>, primitive<unsigned>> result;

    int* hit = reinterpret_cast<int*>(&result.hit);
    int* prim_id = reinterpret_cast<int*>(&result.prim_id);
    int* geom_id = reinterpret_cast<int*>(&result.geom_id);
    float* t = reinterpret_cast<float*>(&result.t);
    array<vec3, N> isect_pos;
    float* u = reinterpret_cast<float*>(&result.u);
    float* v = reinterpret_cast<float*>(&result.v);

    for (size_t i = 0; i < N; ++i)
    {
        hit[i]       = hrs[i].hit ? 0xFFFFFFFF : 0x00000000;
        prim_id[i]   = hrs[i].prim_id;
        geom_id[i]   = hrs[i].geom_id;
        t[i]         = hrs[i].t;
        isect_pos[i] = hrs[i].isect_pos;
        u[i]         = hrs[i].u;
        v[i]         = hrs[i].v;
    }

    result.isect_pos = pack(isect_pos);

    return result;
}

template <
    typename FloatT,
    typename = typename std::enable_if<is_simd_vector<FloatT>::value>::type
    >
MATH_FUNC
inline array<hit_record<ray, primitive<unsigned>>, num_elements<FloatT>::value> unpack(
        hit_record<basic_ray<FloatT>, primitive<unsigned>> const& hr
        )
{
    using float_array = aligned_array_t<FloatT>;
    using int_array = aligned_array_t<int_type_t<FloatT>>;

    int_array hit;
    store(hit, convert_to_int(hr.hit));

    int_array prim_id;
    store(prim_id, hr.prim_id);

    int_array geom_id;
    store(geom_id, hr.geom_id);

    float_array t;
    store(t, hr.t);

    auto isect_pos = unpack(hr.isect_pos);

    float_array u;
    store(u, hr.u);

    float_array v;
    store(v, hr.v);

    array<hit_record<ray, primitive<unsigned>>, num_elements<FloatT>::value> result;
    for (size_t i = 0; i < num_elements<FloatT>::value; ++i)
    {
        result[i].hit       = hit[i] != 0;
        result[i].prim_id   = prim_id[i];
        result[i].geom_id   = geom_id[i];
        result[i].t         = t[i];
        result[i].isect_pos = isect_pos[i];
        result[i].u         = u[i];
        result[i].v         = v[i];
    }
    return result;
}

// TODO: consolidate!
#if VSNRAY_SIMD_ISA_GE(VSNRAY_SIMD_ISA_AVX512F)
template <>
inline array<hit_record<ray, primitive<unsigned>>, num_elements<simd::float16>::value> unpack(
        hit_record<basic_ray<simd::float16>, primitive<unsigned>> const& hr
        )
{
    using FloatT = simd::float16;
    using float_array = aligned_array_t<FloatT>;
    using int_array = aligned_array_t<int_type_t<FloatT>>;
    using mask_array = aligned_array_t<mask_type_t<FloatT>>;

    mask_array hit;
    store(hit, hr.hit);

    int_array prim_id;
    store(prim_id, hr.prim_id);

    int_array geom_id;
    store(geom_id, hr.geom_id);

    float_array t;
    store(t, hr.t);

    auto isect_pos = unpack(hr.isect_pos);

    float_array u;
    store(u, hr.u);

    float_array v;
    store(v, hr.v);

    array<hit_record<ray, primitive<unsigned>>, num_elements<FloatT>::value> result;
    for (size_t i = 0; i < num_elements<FloatT>::value; ++i)
    {
        result[i].hit       = hit[i] != 0;
        result[i].prim_id   = prim_id[i];
        result[i].geom_id   = geom_id[i];
        result[i].t         = t[i];
        result[i].isect_pos = isect_pos[i];
        result[i].u         = u[i];
        result[i].v         = v[i];
    }
    return result;
}
#endif

} // simd

} // MATH_NAMESPACE

#endif // VSNRAY_MATH_INTERSECT_H
