// This file is distributed under the MIT license.
// See the LICENSE file for details.

#pragma once

#ifndef VSNRAY_MATH_UNORM_H
#define VSNRAY_MATH_UNORM_H 1

#include <cstdint>

#include "config.h"

namespace MATH_NAMESPACE
{
namespace detail
{

//-------------------------------------------------------------------------------------------------
// Best unsigned integer type from bit depth. TODO: support arbitrary bit depth
//

template <unsigned Bits> struct best_uint { /* No best_uint::type causes a compiler error */ };
template <> struct best_uint< 8> { typedef uint8_t  type; };
template <> struct best_uint<16> { typedef uint16_t type; };
template <> struct best_uint<32> { typedef uint32_t type; };

} // detail


//-------------------------------------------------------------------------------------------------
// unorm
//

template <unsigned Bits>
class unorm
{
public:

    using value_type = typename detail::best_uint<Bits>::type;

public:

    value_type value;

    unorm() = default;

    MATH_FUNC /* implicit */ unorm(float f);

    MATH_FUNC operator value_type() const;
    MATH_FUNC operator float() const;
};

} // MATH_NAMESPACE

#include "detail/unorm.inl"

#endif // VSNRAY_MATH_UNORM_H
