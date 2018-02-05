// This file is distributed under the MIT license.
// See the LICENSE file for details.

#include <cmath>
#include <type_traits>

namespace visionaray
{
namespace detail
{

template <typename T>
struct widen;

template <>
struct widen<int8_t>
{
    using type = int16_t;
};

template <>
struct widen<int16_t>
{
    using type = int32_t;
};

template <>
struct widen<int32_t>
{
    using type = int64_t;
};

template <>
struct widen<int64_t>
{
    using type = __int128;
};

} // detail


//-------------------------------------------------------------------------------------------------
// fixed members
//

template <unsigned I, unsigned F>
MATH_FUNC
inline typename fixed<I, F>::rep_type fixed<I, F>::get_rep_() const
{
    return rep_;
}

template <unsigned I, unsigned F>
MATH_FUNC
inline void fixed<I,F>::set_rep_(typename fixed<I, F>::rep_type v)
{
    rep_ = v;
}


template <unsigned I, unsigned F>
MATH_FUNC
inline fixed<I, F>::fixed(char c)
    : rep_(c << F)
{
}

template <unsigned I, unsigned F>
MATH_FUNC
inline fixed<I, F>::fixed(short s)
    : rep_(s << F)
{
}

//template <unsigned I, unsigned F>
//MATH_FUNC
//inline fixed<I, F>::fixed(int i)
//    : rep_(i << F)
//{
//}

template <unsigned I, unsigned F>
MATH_FUNC
inline fixed<I, F>::fixed(int i)
{
	rep_ = i;
	rep_ <<= F;
}

template <unsigned I, unsigned F>
MATH_FUNC
inline fixed<I, F>::fixed(long l)
    : rep_(l << F)
{
}

template <unsigned I, unsigned F>
MATH_FUNC
inline fixed<I, F>::fixed(long long ll)
    : rep_(ll << F)
{
}

template <unsigned I, unsigned F>
MATH_FUNC
inline fixed<I, F>::fixed(unsigned char uc)
    : rep_(uc << F)
{
}

template <unsigned I, unsigned F>
MATH_FUNC
inline fixed<I, F>::fixed(unsigned short us)
    : rep_(us << F)
{
}

template <unsigned I, unsigned F>
MATH_FUNC
inline fixed<I, F>::fixed(unsigned int ui)
    : rep_(ui << F)
{
}

template <unsigned I, unsigned F>
MATH_FUNC
inline fixed<I, F>::fixed(unsigned long ul)
    : rep_(ul << F)
{
}

template <unsigned I, unsigned F>
MATH_FUNC
inline fixed<I, F>::fixed(unsigned long long ull)
    : rep_(ull << F)
{
}

template <unsigned I, unsigned F>
MATH_FUNC
inline fixed<I, F>::fixed(float f)
    : rep_(f * (rep_type(1) << F))
{
}

template <unsigned I, unsigned F>
MATH_FUNC
inline fixed<I, F>::fixed(double d)
    : rep_(d * (rep_type(1) << F))
{
}

template <unsigned I, unsigned F>
MATH_FUNC
inline fixed<I, F>::fixed(long double ld)
    : rep_(ld * (rep_type(1) << F))
{
}

template <unsigned I, unsigned F>
template <unsigned A, unsigned B>
MATH_FUNC
inline fixed<I, F>::fixed(const fixed<A, B>& fix)
{
	using dst = fixed<I, F>;
	using src = fixed<A, B>;
	
	int diff = F - B;	// number of added decimal bits
	
	if (sizeof(typename src::rep_type) >= sizeof(typename dst::rep_type))
//	if ((I+F) > (A+B)) // which is faster?
	{	// shrink: shift first then store
		if (diff < 0)
			rep_ = static_cast<typename dst::rep_type>(fix.get_rep_() >> -diff);
		else
			rep_ = static_cast<typename dst::rep_type>(fix.get_rep_() << diff);
	}
	else
	{	// widen: store first then shift
		if (diff < 0)
			rep_ = static_cast<typename dst::rep_type>(fix.get_rep_()) >> -diff;
		else
			rep_ = static_cast<typename dst::rep_type>(fix.get_rep_()) << diff;
	}
}

template <unsigned I, unsigned F>
MATH_FUNC
inline fixed<I, F>::operator char() const
{
    return static_cast<char>(rep_) >> F;
}

template <unsigned I, unsigned F>
MATH_FUNC
inline fixed<I, F>::operator short() const
{
    return static_cast<short>(rep_) >> F;
}

template <unsigned I, unsigned F>
MATH_FUNC
inline fixed<I, F>::operator int() const
{
    return static_cast<int>(rep_) >> F;
}

template <unsigned I, unsigned F>
MATH_FUNC
inline fixed<I, F>::operator long() const
{
    return static_cast<long>(rep_) >> F;
}

template <unsigned I, unsigned F>
MATH_FUNC
inline fixed<I, F>::operator long long() const
{
    return static_cast<long long>(rep_) >> F;
}

template <unsigned I, unsigned F>
MATH_FUNC
inline fixed<I, F>::operator unsigned char() const
{
    return static_cast<unsigned char>(rep_) >> F;
}

template <unsigned I, unsigned F>
MATH_FUNC
inline fixed<I, F>::operator unsigned short() const
{
    return static_cast<unsigned short>(rep_) >> F;
}

template <unsigned I, unsigned F>
MATH_FUNC
inline fixed<I, F>::operator unsigned int() const
{
    return static_cast<unsigned int>(rep_) >> F;
}

template <unsigned I, unsigned F>
MATH_FUNC
inline fixed<I, F>::operator unsigned long() const
{
    return static_cast<unsigned long>(rep_) >> F;
}

template <unsigned I, unsigned F>
MATH_FUNC
inline fixed<I, F>::operator unsigned long long() const
{
    return static_cast<unsigned long long>(rep_) >> F;
}

template <unsigned I, unsigned F>
MATH_FUNC
inline fixed<I, F>::operator float() const
{
    return static_cast<float>(rep_) / (rep_type(1) << F);
}

template <unsigned I, unsigned F>
MATH_FUNC
inline fixed<I, F>::operator double() const
{
    return static_cast<double>(rep_) / (rep_type(1) << F);
}

template <unsigned I, unsigned F>
MATH_FUNC
inline fixed<I, F>::operator long double() const
{
    return static_cast<long double>(rep_) / (rep_type(1) << F);
}
/*
template <unsigned I, unsigned F>
template <unsigned A, unsigned B>
MATH_FUNC
inline fixed<I, F>::operator fixed<A, B>() const
{
	using src = fixed<I, F>;
	using dst = fixed<A, B>;
	
	int diff = B - F;	// number of added decimal bits
	fixed<A, B> res;
	
	if (sizeof(typename src::rep_type) >= sizeof(typename dst::rep_type))
//	if ((I+F) > (A+B)) // which is faster?
	{	// shrink: shift first then store
		if (diff < 0)
			res.rep_ = static_cast<typename dst::rep_type>(rep_ >> -diff);
		else
			res.rep_ = static_cast<typename dst::rep_type>(rep_ << diff);
	}
	else
	{	// widen: store first then shift
		if (diff < 0)
			res.rep_ = static_cast<typename dst::rep_type>(rep_) >> -diff;
		else
			res.rep_ = static_cast<typename dst::rep_type>(rep_) << diff;
	}
	return res;
}
*/




//-------------------------------------------------------------------------------------------------
// Basic arithmetic
//

template <unsigned I, unsigned F>
MATH_FUNC
inline fixed<I, F> operator-(fixed<I, F> const& a)
{
    using T = typename fixed<I, F>::rep_type;

    T i = *reinterpret_cast<T const*>(&a);
    T k = ~i;

    return *reinterpret_cast<fixed<I, F>*>(&k);
}

template <unsigned I, unsigned F>
MATH_FUNC
inline fixed<I, F> operator+(fixed<I, F> const& a, fixed<I, F> const& b)
{
    using T = typename fixed<I, F>::rep_type;

    T i = *reinterpret_cast<T const*>(&a);
    T j = *reinterpret_cast<T const*>(&b);

    T k = i + j;

    return *reinterpret_cast<fixed<I, F>*>(&k);
}

template <unsigned I, unsigned F>
MATH_FUNC
inline fixed<I, F> operator-(fixed<I, F> const& a, fixed<I, F> const& b)
{
    using T = typename fixed<I, F>::rep_type;

    T i = *reinterpret_cast<T const*>(&a);
    T j = *reinterpret_cast<T const*>(&b);

    T k = i - j;

    return *reinterpret_cast<fixed<I, F>*>(&k);
}

template <unsigned I, unsigned F>
MATH_FUNC
inline fixed<I, F> operator*(fixed<I, F> const& a, fixed<I, F> const& b)
{
    using T = typename fixed<I, F>::rep_type;
    using WT = typename detail::widen<T>::type;

    T i = *reinterpret_cast<T const*>(&a);
    T j = *reinterpret_cast<T const*>(&b);

    WT k = static_cast<WT>(i) * j >> F;

    return *reinterpret_cast<fixed<I, F>*>(&k);
}

template <unsigned I, unsigned F>
MATH_FUNC
inline fixed<I, F> operator/(fixed<I, F> const& a, fixed<I, F> const& b)
{
    using T = typename fixed<I, F>::rep_type;
    using WT = typename detail::widen<T>::type;

    T i = *reinterpret_cast<T const*>(&a);
    T j = *reinterpret_cast<T const*>(&b);

    WT k = (static_cast<WT>(i) << F) / j;

    return *reinterpret_cast<fixed<I, F>*>(&k);
}


//--------------------------------------------------------------------------------------------------
// Comparisons
//

template <unsigned I, unsigned F>
MATH_FUNC
inline bool operator==(fixed<I, F> const& a, fixed<I, F> const& b)
{
    using T = typename fixed<I, F>::rep_type;

    T i = *reinterpret_cast<T const*>(&a);
    T j = *reinterpret_cast<T const*>(&b);

    return i == j;
}

template <unsigned I, unsigned F>
MATH_FUNC
inline bool operator!=(fixed<I, F> const& a, fixed<I, F> const& b)
{
    using T = typename fixed<I, F>::rep_type;

    T i = *reinterpret_cast<T const*>(&a);
    T j = *reinterpret_cast<T const*>(&b);

    return i != j;
}

template <unsigned I, unsigned F>
MATH_FUNC
inline bool operator>(fixed<I, F> const& a, fixed<I, F> const& b)
{
    using T = typename fixed<I, F>::rep_type;

    T i = *reinterpret_cast<T const*>(&a);
    T j = *reinterpret_cast<T const*>(&b);

    return i > j;
}

template <unsigned I, unsigned F>
MATH_FUNC
inline bool operator<(fixed<I, F> const& a, fixed<I, F> const& b)
{
    using T = typename fixed<I, F>::rep_type;

    T i = *reinterpret_cast<T const*>(&a);
    T j = *reinterpret_cast<T const*>(&b);

    return i < j;
}

template <unsigned I, unsigned F>
MATH_FUNC
inline bool operator>=(fixed<I, F> const& a, fixed<I, F> const& b)
{
    using T = typename fixed<I, F>::rep_type;

    T i = *reinterpret_cast<T const*>(&a);
    T j = *reinterpret_cast<T const*>(&b);

    return i >= j;
}

template <unsigned I, unsigned F>
MATH_FUNC
inline bool operator<=(fixed<I, F> const& a, fixed<I, F> const& b)
{
    using T = typename fixed<I, F>::rep_type;

    T i = *reinterpret_cast<T const*>(&a);
    T j = *reinterpret_cast<T const*>(&b);

    return i <= j;
}


//-------------------------------------------------------------------------------------------------
// Math functions
//

template <unsigned I, unsigned F>
MATH_FUNC
inline fixed<I, F> abs(fixed<I, F> const& a)
{
    using T = typename fixed<I, F>::rep_type;

    T i = *reinterpret_cast<T const*>(&a);

    T bits = i < 0 ? T(-1) : T(0);
    T k = (i + bits) ^ bits;

    return *reinterpret_cast<fixed<I, F>*>(&k);
}

template <unsigned I, unsigned F>
MATH_FUNC
inline fixed<I, F> floor(fixed<I, F> const& a)
{
    using T = typename fixed<I, F>::rep_type;

    T i = *reinterpret_cast<T const*>(&a);

    T k = ((i >> F) << F);

    return *reinterpret_cast<fixed<I, F>*>(&k);
}

template <unsigned I, unsigned F>
MATH_FUNC
inline fixed<I, F> rsqrt(fixed<I, F> const& a)
{
    static const int K = 5;
    static const fixed<I, F> threehalf = 1.5f;
    fixed<I, F> ahalf = a * fixed<I, F>(0.5f);
    fixed<I, F> t = a;

    for (int k = 0; k < K; ++k)
    {
        t = t * (threehalf - ahalf * t * t);
    }

    return t;
}

template <unsigned I, unsigned F>
MATH_FUNC
inline fixed<I, F> sqrt(fixed<I, F> const& a)
{
    static const int K = 8;
    static const fixed<I, F> half = 0.5f;

    fixed<I, F> x = 1.0f;

    for (int k = 0; k < K; ++k)
    {
        x = half * (x + a / x);
    }

    return x;
}

template <unsigned I, unsigned F>
MATH_FUNC
inline fixed<I, F> trunc(fixed<I, F> const& a)
{
    using T = typename fixed<I, F>::rep_type;

    T i = *reinterpret_cast<T const*>(&a);

    T sgn = (i >> (I + F - 1)) & 1;
    T k = ((i >> F) << F) + sgn * (T(1) << F);

    return *reinterpret_cast<fixed<I, F>*>(&k);
}

} // visionaray
