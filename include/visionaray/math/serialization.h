// This file is distributed under the MIT license.
// See the LICENSE file for details.

#ifndef VISIONARAY_MATH_SERIALIZATION_H
#define VISIONARAY_MATH_SERIALIZATION_H


#include "math.h"


namespace boost
{


namespace serialization
{


template <typename A, size_t D, typename T>
inline void serialize(A& a, MATH_NAMESPACE::vector<D, T>& v, unsigned /* version */ )
{
    for (size_t d = 0; d < D; ++d)
    {
        a & v[d];
    }
}


template <typename A, typename T>
inline void serialize(A& a, MATH_NAMESPACE::vector<2, T>& v, unsigned /* version */ )
{
    a & v.x;
    a & v.y;
}


template <typename A, typename T>
inline void serialize(A& a, MATH_NAMESPACE::vector<3, T>& v, unsigned /* version */ )
{
    a & v.x;
    a & v.y;
    a & v.z;
}


template <typename A, typename T>
inline void serialize(A& a, MATH_NAMESPACE::vector<4, T>& v, unsigned /* version */ )
{
    a & v.x;
    a & v.y;
    a & v.z;
    a & v.w;
}


template <typename A, typename T>
inline void serialize(A& a, MATH_NAMESPACE::matrix<4, 4, T>& m, unsigned /* version */ )
{
    a & m.col0;
    a & m.col1;
    a & m.col2;
    a & m.col3;
}


template <typename A, typename T>
inline void serialize(A& a, MATH_NAMESPACE::base_aabb<T>& b, unsigned /* version */ )
{
    a & b.min;
    a & b.max;
}


template <typename A, typename T>
inline void serialize(A& a, MATH_NAMESPACE::rectangle<MATH_NAMESPACE::xywh_layout, T>& r, unsigned /* version */ )
{
    a & r.x;
    a & r.y;
    a & r.w;
    a & r.h;
}


} // serialization


} // boost


#endif // VISIONARAY_MATH_SERIALIZATION_H


