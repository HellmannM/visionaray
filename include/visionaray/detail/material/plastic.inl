// This file is distributed under the MIT license.
// See the LICENSE file for details.

#include <visionaray/math/constants.h>

namespace visionaray
{

//-------------------------------------------------------------------------------------------------
// Public interface
//

template <typename T>
VSNRAY_FUNC
inline spectrum<T> plastic<T>::ambient() const
{
    return ca_ * ka_;
}

template <typename T>
template <typename SR>
VSNRAY_FUNC
inline spectrum<typename SR::scalar_type> plastic<T>::shade(SR const& sr) const
{
    using U = typename SR::scalar_type;

    auto l = sr.light;
    auto wi = sr.light_dir;
    auto wo = sr.view_dir;
    auto n = sr.normal;
#if 1 // two-sided
    n = faceforward( n, sr.view_dir, sr.geometric_normal );
#endif
    auto ndotl = max( U(0.0), dot(n, wi) );

    return spectrum<U>(
            constants::pi<U>()
          * ( cd_impl(sr, n, wo, wi) + specular_brdf_.f(n, wo, wi) )
          * spectrum<U>(from_rgb(l.intensity(sr.isect_pos)))
          * ndotl
            );
}

namespace detail
{

template <typename T, typename SR>
VSNRAY_FUNC
spectrum<T> tex_color_from_shade_record(SR)
{
    return spectrum<T>(1.0f);
}

template <typename T, typename L, typename C, typename U>
VSNRAY_FUNC
spectrum<T> tex_color_from_shade_record(shade_record<L, C, U> const& sr)
{
    return spectrum<T>(from_rgb(sr.tex_color));
}

} // detail

template <typename T>
template <typename SR, typename U, typename Sampler>
VSNRAY_FUNC
inline spectrum<U> plastic<T>::sample(
        SR const&       shade_rec,
        vector<3, U>&   refl_dir,
        U&              pdf,
        Sampler&        sampler
        ) const
{
    U pdf1(0.0);
    U pdf2(0.0);

    vector<3, U> refl1(0.0);
    vector<3, U> refl2(0.0);

    spectrum<U>  diff(0.0);
    spectrum<U>  spec(0.0);

    auto prob_diff = mean_value( diffuse_brdf_.cd ) * diffuse_brdf_.kd;
    auto prob_spec = mean_value( specular_brdf_.cs ) * specular_brdf_.ks;

    auto all_zero  = prob_diff == U(0.0) && prob_spec == U(0.0);

    prob_diff      = select( all_zero, U(0.5), prob_diff );
    prob_spec      = select( all_zero, U(0.5), prob_spec );

    prob_diff      = prob_diff / (prob_diff + prob_spec);


    auto u         = sampler.next();

    auto n = shade_rec.normal;
#if 1 // two-sided
    n = faceforward( n, shade_rec.view_dir, shade_rec.geometric_normal );
#endif

    if (any(u < U(prob_diff)))
    {
        diff       = detail::tex_color_from_shade_record<U>(shade_rec) * diffuse_brdf_.sample_f(n, shade_rec.view_dir, refl1, pdf1, sampler);
    }

    if (any(u >= U(prob_diff)))
    {
        spec       = specular_brdf_.sample_f(n, shade_rec.view_dir, refl2, pdf2, sampler);
    }

    pdf            = select( u < U(prob_diff), pdf1,  pdf2  );
    refl_dir       = select( u < U(prob_diff), refl1, refl2 );

    return           select( u < U(prob_diff), diff,  spec  );
}

//--- deprecated begin ------------------------------------

template <typename T>
VSNRAY_FUNC
inline void plastic<T>::set_ca(spectrum<T> const& ca)
{
    ca_ = ca;
}

template <typename T>
VSNRAY_FUNC
inline spectrum<T> plastic<T>::get_ca() const
{
    return ca_;
}

template <typename T>
VSNRAY_FUNC
inline void plastic<T>::set_ka(T ka)
{
    ka_ = ka;
}

template <typename T>
VSNRAY_FUNC
inline T plastic<T>::get_ka() const
{
    return ka_;
}

template <typename T>
VSNRAY_FUNC
inline void plastic<T>::set_cd(spectrum<T> const& cd)
{
    diffuse_brdf_.cd = cd;
}

template <typename T>
VSNRAY_FUNC
inline spectrum<T> plastic<T>::get_cd() const
{
    return diffuse_brdf_.cd;
}

template <typename T>
VSNRAY_FUNC
inline void plastic<T>::set_kd(T kd)
{
    diffuse_brdf_.kd = kd;
}

template <typename T>
VSNRAY_FUNC
inline T plastic<T>::get_kd() const
{
    return diffuse_brdf_.kd;
}

template <typename T>
VSNRAY_FUNC
inline void plastic<T>::set_cs(spectrum<T> const& cs)
{
    specular_brdf_.cs = cs;
}

template <typename T>
VSNRAY_FUNC
inline spectrum<T> plastic<T>::get_cs() const
{
    return specular_brdf_.cs;
}

template <typename T>
VSNRAY_FUNC
inline void plastic<T>::set_ks(T ks)
{
    specular_brdf_.ks = ks;
}

template <typename T>
VSNRAY_FUNC
inline T plastic<T>::get_ks() const
{
    return specular_brdf_.ks;
}

template <typename T>
VSNRAY_FUNC
inline void plastic<T>::set_specular_exp(T exp)
{
    specular_brdf_.exp = exp;
}

template <typename T>
VSNRAY_FUNC
inline T plastic<T>::get_specular_exp() const
{
    return specular_brdf_.exp;
}

//--- deprecated end --------------------------------------

template <typename T>
VSNRAY_FUNC
inline spectrum<T>& plastic<T>::ca()
{
    return ca_;
}

template <typename T>
VSNRAY_FUNC
inline spectrum<T> const& plastic<T>::ca() const
{
    return ca_;
}

template <typename T>
VSNRAY_FUNC
inline T& plastic<T>::ka()
{
    return ka_;
}

template <typename T>
VSNRAY_FUNC
inline T const& plastic<T>::ka() const
{
    return ka_;
}

template <typename T>
VSNRAY_FUNC
inline spectrum<T>& plastic<T>::cd()
{
    return diffuse_brdf_.cd;
}

template <typename T>
VSNRAY_FUNC
inline spectrum<T> const& plastic<T>::cd() const
{
    return diffuse_brdf_.cd;
}

template <typename T>
VSNRAY_FUNC
inline T& plastic<T>::kd()
{
    return diffuse_brdf_.kd;
}

template <typename T>
VSNRAY_FUNC
inline T const& plastic<T>::kd() const
{
    return diffuse_brdf_.kd;
}

template <typename T>
VSNRAY_FUNC
inline spectrum<T>& plastic<T>::cs()
{
    return specular_brdf_.cs;
}

template <typename T>
VSNRAY_FUNC
inline spectrum<T> const& plastic<T>::cs() const
{
    return specular_brdf_.cs;
}

template <typename T>
VSNRAY_FUNC
inline T& plastic<T>::ks()
{
    return specular_brdf_.ks;
}

template <typename T>
VSNRAY_FUNC
inline T const& plastic<T>::ks() const
{
    return specular_brdf_.ks;
}

template <typename T>
VSNRAY_FUNC
inline T& plastic<T>::specular_exp()
{
    return specular_brdf_.exp;
}

template <typename T>
VSNRAY_FUNC
inline T const& plastic<T>::specular_exp() const
{
    return specular_brdf_.exp;
}


//-------------------------------------------------------------------------------------------------
// Private functions
//

template <typename T>
template <typename SR, typename V>
VSNRAY_FUNC
inline spectrum<T> plastic<T>::cd_impl(SR const& sr, V const& n, V const& wo, V const& wi) const
{
    VSNRAY_UNUSED(sr);
    return diffuse_brdf_.f(n, wo, wi);
}

template <typename T>
template <typename L, typename C, typename S, typename V>
VSNRAY_FUNC
inline spectrum<T> plastic<T>::cd_impl(shade_record<L, C, S> const& sr, V const& n, V const& wo, V const& wi) const
{
    return spectrum<T>(from_rgb(sr.tex_color)) * diffuse_brdf_.f(n, wo, wi);
}

} // visionaray
