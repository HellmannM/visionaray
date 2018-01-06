// This file is distributed under the MIT license.
// See the LICENSE file for details.

#pragma once

#ifndef VSNRAY_MATH_TRIANGLE_H
#define VSNRAY_MATH_TRIANGLE_H 1

#include <cstddef>

#include "primitive.h"
#include "vector.h"

#define HLS_NO_XIL_FPO_LIB
#include <ap_int.h>

namespace MATH_NAMESPACE
{
/*
template <size_t Dim, typename T, typename P>
class basic_triangle : public primitive<P>
{
public:

    using scalar_type =  T;
    using vec_type    =  vector<Dim, T>;

public:

    MATH_FUNC basic_triangle() = default;
    MATH_FUNC basic_triangle(
            vector<Dim, T> const& v1,
            vector<Dim, T> const& e1,
            vector<Dim, T> const& e2
            );

    vec_type v1;
    vec_type e1;
    vec_type e2;
};
*/

template <size_t Dim, typename T, typename P>
class basic_triangle : public primitive<P>
{
public:

    using F = T;
    using D = ap_fixed<32, 1, AP_RND>;
    using S = ap_fixed<32, 16, AP_RND>;
    using scalar_type =  T;
    using vec_type    =  vector<Dim, T>;

public:

    MATH_FUNC basic_triangle() = default;
    MATH_FUNC basic_triangle(
            vector<Dim, T> const& v1,
            vector<Dim, T> const& e1,
            vector<Dim, T> const& e2
            );

    vec_type v1;
    vec_type e1;
    vec_type e2;

    uint32_t i0;
    vector<3, F> n;
    F pp, pq;
    D np, nq;
    S e1p, e2p, e1q, e2q;
    S d;

};

} // MATH_NAMESPACE

#include "detail/triangle.inl"

#endif // VSNRAY_MATH_TRIANGLE_H
