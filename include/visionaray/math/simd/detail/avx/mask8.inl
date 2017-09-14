// This file is distributed under the MIT license.
// See the LICENSE file for details.

namespace MATH_NAMESPACE
{
namespace simd
{

//-------------------------------------------------------------------------------------------------
// mask8 members
//

VSNRAY_FORCE_INLINE mask8::basic_mask(__m256 const& m)
    : f(m)
{
}

VSNRAY_FORCE_INLINE mask8::basic_mask(__m256i const& m)
    : i(m)
{
}

VSNRAY_FORCE_INLINE mask8::basic_mask(bool b)
    : i( basic_int<__m256i>(b ? 0xFFFFFFFF : 0x00000000) )
{
}

VSNRAY_FORCE_INLINE mask8::basic_mask(
        bool x1,
        bool x2,
        bool x3,
        bool x4,
        bool x5,
        bool x6,
        bool x7,
        bool x8
        )
    : i( basic_int<__m256i>(
        x1 ? 0xFFFFFFFF : 0x00000000,
        x2 ? 0xFFFFFFFF : 0x00000000,
        x3 ? 0xFFFFFFFF : 0x00000000,
        x4 ? 0xFFFFFFFF : 0x00000000,
        x5 ? 0xFFFFFFFF : 0x00000000,
        x6 ? 0xFFFFFFFF : 0x00000000,
        x7 ? 0xFFFFFFFF : 0x00000000,
        x8 ? 0xFFFFFFFF : 0x00000000
        ) )
{
}

VSNRAY_FORCE_INLINE mask8::basic_mask(bool const v[8])
    : i( basic_int<__m256i>(
        v[0] ? 0xFFFFFFFF : 0x00000000,
        v[1] ? 0xFFFFFFFF : 0x00000000,
        v[2] ? 0xFFFFFFFF : 0x00000000,
        v[3] ? 0xFFFFFFFF : 0x00000000,
        v[4] ? 0xFFFFFFFF : 0x00000000,
        v[5] ? 0xFFFFFFFF : 0x00000000,
        v[6] ? 0xFFFFFFFF : 0x00000000,
        v[7] ? 0xFFFFFFFF : 0x00000000
        ) )
{
}

//-------------------------------------------------------------------------------------------------
// Static cast
//

VSNRAY_FORCE_INLINE int8 convert_to_int(mask8 const& a)
{
    return a.i;
}


//-------------------------------------------------------------------------------------------------
// any / all intrinsics
//

VSNRAY_FORCE_INLINE bool any(mask8 const& m)
{
    return _mm256_movemask_ps(m.f) != 0;
}

VSNRAY_FORCE_INLINE bool all(mask8 const& m)
{
    return _mm256_movemask_ps(m.f) == 0xFF;
}


//-------------------------------------------------------------------------------------------------
// select intrinsic
//

VSNRAY_FORCE_INLINE mask8 select(mask8 const& m, mask8 const& a, mask8 const& b)
{
    return _mm256_blendv_ps(b.f, a.f, m.f);
}


//-------------------------------------------------------------------------------------------------
// Load / store
//

template <typename S, typename T>
VSNRAY_FORCE_INLINE void store(S dst[8], T const& v, mask8 const& mask)
{
    T old(dst);
    store( dst, select(mask, v, old) );
}

template <typename S, typename T>
VSNRAY_FORCE_INLINE void store(S dst[8], T const& v, mask8 const& mask, T const& old)
{
    store( dst, select(mask, v, old) );
}


//-------------------------------------------------------------------------------------------------
// Bitwise operations
//

VSNRAY_FORCE_INLINE mask8 operator!(mask8 const& a)
{
    return _mm256_xor_ps(a.f, mask8(true).f);
}

VSNRAY_FORCE_INLINE mask8 operator&(mask8 const& a, mask8 const& b)
{
    return _mm256_and_ps(a.f, b.f);
}

VSNRAY_FORCE_INLINE mask8 operator|(mask8 const& a, mask8 const& b)
{
    return _mm256_or_ps(a.f, b.f);
}

VSNRAY_FORCE_INLINE mask8 operator^(mask8 const& a, mask8 const& b)
{
    return _mm256_xor_ps(a.f, b.f);
}


//-------------------------------------------------------------------------------------------------
// Logical operations
//

VSNRAY_FORCE_INLINE mask8 operator&&(mask8 const& a, mask8 const& b)
{
    // Ok because masks only store booleans
    return a & b;
}

VSNRAY_FORCE_INLINE mask8 operator||(mask8 const& a, mask8 const& b)
{
    // Ok because masks only store booleans
    return a | b;
}


//-------------------------------------------------------------------------------------------------
// Comparisons
//

VSNRAY_FORCE_INLINE mask8 operator==(mask8 const& u, mask8 const& v)
{
    return _mm256_cmp_ps(float8(u.i), float8(v.i), _CMP_EQ_OQ);
}

VSNRAY_FORCE_INLINE mask8 operator!=(mask8 const& u, mask8 const& v)
{
    return !(u == v);
}

} // simd


//-------------------------------------------------------------------------------------------------
// Import SIMD intrinsics into namespace visionaray.
// Enable ADL!
//

using simd::select;
using simd::store;
using simd::any;
using simd::all;

} // MATH_NAMESPACE
