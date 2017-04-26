// This file is distributed under the MIT license.
// See the LICENSE file for details.

namespace MATH_NAMESPACE
{
namespace simd
{

//-------------------------------------------------------------------------------------------------
// Basic arithmetic
//

template <typename T, typename U>
MATH_FUNC
VSNRAY_FORCE_INLINE basic_int<T>& operator+=(basic_int<T>& a, U const& b)
{
    a = a + b;
    return a;
}

template <typename T, typename U>
MATH_FUNC
VSNRAY_FORCE_INLINE basic_int<T>& operator-=(basic_int<T>& a, U const& b)
{
    a = a - b;
    return a;
}

template <typename T, typename U>
MATH_FUNC
VSNRAY_FORCE_INLINE basic_int<T>& operator*=(basic_int<T>& a, U const& b)
{
    a = a * b;
    return a;
}

template <typename T, typename U>
MATH_FUNC
VSNRAY_FORCE_INLINE basic_int<T>& operator/=(basic_int<T>& a, U const& b)
{
    a = a / b;
    return a;
}


//-------------------------------------------------------------------------------------------------
// Bitwise operators
//

template <typename T, typename U>
MATH_FUNC
VSNRAY_FORCE_INLINE basic_int<T>& operator&=(basic_int<T>& a, U const& b)
{
    a = a & b;
    return a;
}

template <typename T, typename U>
MATH_FUNC
VSNRAY_FORCE_INLINE basic_int<T>& operator|=(basic_int<T>& a, U const& b)
{
    a = a | b;
    return a;
}

template <typename T, typename U>
MATH_FUNC
VSNRAY_FORCE_INLINE basic_int<T>& operator^=(basic_int<T>& a, U const& b)
{
    a = a ^ b;
    return a;
}

template <typename T>
MATH_FUNC
VSNRAY_FORCE_INLINE basic_int<T>& operator<<=(basic_int<T>& a, int count)
{
    a << count;
    return a;
}

template <typename T>
MATH_FUNC
VSNRAY_FORCE_INLINE basic_int<T>& operator>>=(basic_int<T>& a, int count)
{
    a >> count;
    return a;
}

} // simd
} // MATH_NAMESPACE
