// Copyright (C) 2018-2024 - DevSH Graphics Programming Sp. z O.O.
// This file is part of the "Nabla Engine".
// For conditions of distribution and use, see copyright notice in nabla.h

#ifndef _NBL_BUILTIN_HLSL_BXDF_BSDF_SPECULAR_GGX_INCLUDED_
#define _NBL_BUILTIN_HLSL_BXDF_BSDF_SPECULAR_GGX_INCLUDED_

#include <nbl/builtin/hlsl/bxdf/transmission.hlsl>
#include <nbl/builtin/hlsl/bxdf/ndf/ggx.hlsl>
#include <nbl/builtin/hlsl/bxdf/geom/smith/ggx.hlsl>
#include <nbl/builtin/hlsl/bxdf/brdf/specular/ggx.hlsl>


namespace nbl
{
namespace hlsl
{
namespace bxdf
{
namespace bsdf
{
namespace specular
{

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar ggx_height_correlated_aniso_cos_eval_wo_clamps(
    Scalar NdotH2, Scalar TdotH2, Scalar BdotH2,
    Scalar absNdotL, Scalar NdotL2, Scalar TdotL2, Scalar BdotL2,
    Scalar absNdotV, Scalar NdotV2, Scalar TdotV2, Scalar BdotV2,
    bool transmitted, Scalar VdotH, Scalar LdotH, Scalar VdotHLdotH,
    Scalar orientedEta, Scalar orientedEta2,
    Scalar ax, Scalar ax2, Scalar ay, Scalar ay2)
{
    Scalar NG_already_in_reflective_dL_measure = brdf::specular::ggx_height_correlated_aniso_cos_eval_DG_wo_clamps(NdotH2, TdotH2, BdotH2, absNdotL, NdotL2, TdotL2, BdotL2, absNdotV, NdotV2, TdotV2, BdotV2, ax, ax2, ay, ay2);

    const Scalar reflectance = fresnel::dielectric_common(orientedEta2, abs(VdotH));

    return reflectance * ndf::ggx::microfacet_to_light_measure_transform(NG_already_in_reflective_dL_measure, absNdotL, transmitted, VdotH, LdotH, VdotHLdotH, orientedEta);
}

// before calling you must ensure that `AnisotropicMicrofacetCache` is valid (if a given V vector can "see" the L vector)
template<class IncomingRayDirInfo>
    NBL_REQUIRES(ray_dir_info::basic<IncomingRayDirInfo>)
typename IncomingRayDirInfo::scalar_t ggx_height_correlated_aniso_cos_eval(
    NBL_CONST_REF_ARG(LightSample<IncomingRayDirInfo>) _sample,
    NBL_CONST_REF_ARG(surface_interactions::Anisotropic<IncomingRayDirInfo>) interaction,
    NBL_CONST_REF_ARG(AnisotropicMicrofacetCache) _cache,
    typename IncomingRayDirInfo::scalar_t  eta,
    typename IncomingRayDirInfo::scalar_t  ax,
    typename IncomingRayDirInfo::scalar_t  ay)
{
    using scalar_t = typename IncomingRayDirInfo::scalar_t ;

    const scalar_t TdotH2 = _cache.TdotH * _cache.TdotH;
    const scalar_t BdotH2 = _cache.BdotH * _cache.BdotH;

    const scalar_t TdotL2 = _sample.TdotL * _sample.TdotL;
    const scalar_t BdotL2 = _sample.BdotL * _sample.BdotL;

    const scalar_t TdotV2 = interaction.TdotV * interaction.TdotV;
    const scalar_t BdotV2 = interaction.BdotV * interaction.BdotV;

    const scalar_t VdotH = _cache.VdotH;

    scalar_t orientedEta, dummy;
    const bool backside = math::getOrientedEtas(orientedEta, dummy, VdotH, eta);
    const scalar_t orientedEta2 = orientedEta * orientedEta;

    const scalar_t VdotHLdotH = VdotH * _cache.LdotH;
    const bool transmitted = VdotHLdotH < 0.0;

    return ggx_height_correlated_aniso_cos_eval_wo_clamps(
        _cache.NdotH2, TdotH2, BdotH2,
        abs(_sample.NdotL), _sample.NdotL2, TdotL2, BdotL2,
        abs(interaction.NdotV), interaction.NdotV_squared, TdotV2, BdotV2,
        transmitted, VdotH, _cache.LdotH, VdotHLdotH, orientedEta, orientedEta2,
        ax, ax * ax, ay, ay * ay
    );
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar ggx_height_correlated_cos_eval_wo_clamps(
    Scalar NdotH2, Scalar absNdotL, Scalar NdotL2,
    Scalar absNdotV, Scalar NdotV2,
    bool transmitted, Scalar VdotH, Scalar LdotH, Scalar VdotHLdotH,
    Scalar orientedEta, Scalar orientedEta2, Scalar a2)
{
    const Scalar NG_already_in_reflective_dL_measure = brdf::specular::ggx_height_correlated_cos_eval_DG_wo_clamps(NdotH2, absNdotL, NdotL2, absNdotV, NdotV2, a2);

    const Scalar reflectance = fresnel::dielectric_common(orientedEta2, abs(VdotH));

    return reflectance * ndf::ggx::microfacet_to_light_measure_transform(NG_already_in_reflective_dL_measure, absNdotL, transmitted, VdotH, LdotH, VdotHLdotH, orientedEta);
}

// before calling you must ensure that `AnisotropicMicrofacetCache` is valid (if a given V vector can "see" the L vector)
template<class IncomingRayDirInfo>
    NBL_REQUIRES(ray_dir_info::basic<IncomingRayDirInfo>)
typename IncomingRayDirInfo::scalar_t ggx_height_correlated_cos_eval(
    NBL_CONST_REF_ARG(LightSample<IncomingRayDirInfo>) _sample,
    NBL_CONST_REF_ARG(surface_interactions::Isotropic<IncomingRayDirInfo>) interaction,
    NBL_CONST_REF_ARG(IsotropicMicrofacetCache) _cache,
    typename IncomingRayDirInfo::scalar_t eta,
    typename IncomingRayDirInfo::scalar_t a2)
{
    using scalar_t = typename IncomingRayDirInfo::scalar_t;

    scalar_t orientedEta, dummy;
    const bool backside = math::getOrientedEtas(orientedEta, dummy, _cache.VdotH, eta);
    const scalar_t orientedEta2 = orientedEta * orientedEta;

    const scalar_t VdotHLdotH = _cache.VdotH * _cache.LdotH;
    const bool transmitted = VdotHLdotH < 0.0;

    return ggx_height_correlated_cos_eval_wo_clamps(
        _cache.NdotH2, abs(_sample.NdotL), _sample.NdotL2,
        abs(interaction.NdotV), interaction.NdotV_squared,
        transmitted, _cache.VdotH, _cache.LdotH, VdotHLdotH, orientedEta, orientedEta2, a2
    );
}

// TODO: unifty the two following functions into `microfacet_BSDF_cos_generate_wo_clamps(float3 H,...)` and `microfacet_BSDF_cos_generate` or at least a auto declaration macro in lieu of a template
template<class IncomingRayDirInfo>
    NBL_REQUIRES(ray_dir_info::basic<IncomingRayDirInfo>)
LightSample<IncomingRayDirInfo> ggx_cos_generate_wo_clamps(
    NBL_CONST_REF_ARG(vector<typename IncomingRayDirInfo::scalar_t, 3>) localV,
    bool backside,
    NBL_CONST_REF_ARG(vector<typename IncomingRayDirInfo::scalar_t, 3>) upperHemisphereLocalV,
    NBL_CONST_REF_ARG(matrix<typename IncomingRayDirInfo::scalar_t, 3, 3>) m,
    NBL_REF_ARG(vector<typename IncomingRayDirInfo::scalar_t, 3>) u,
    typename IncomingRayDirInfo::scalar_t _ax,
    typename IncomingRayDirInfo::scalar_t _ay,
    typename IncomingRayDirInfo::scalar_t rcpOrientedEta,
    typename IncomingRayDirInfo::scalar_t orientedEta2,
    typename IncomingRayDirInfo::scalar_t rcpOrientedEta2,
    NBL_REF_ARG(AnisotropicMicrofacetCache) _cache)
{
    using scalar_t = typename IncomingRayDirInfo::scalar_t;

    // thanks to this manoeuvre the H will always be in the upper hemisphere (NdotH>0.0)
    const vector<scalar_t, 3> H = brdf::specular::ggx_cos_generate(upperHemisphereLocalV, u.xy, _ax, _ay);

    const scalar_t VdotH = dot(localV, H);
    const scalar_t reflectance = fresnel::dielectric_common(orientedEta2, abs(VdotH));

    scalar_t rcpChoiceProb;
    bool transmitted = math::partitionRandVariable(reflectance, u.z, rcpChoiceProb);

    vector<scalar_t, 3> localL;
    _cache = AnisotropicMicrofacetCache::create(localV, H, transmitted, rcpOrientedEta, rcpOrientedEta2);
    localL = math::reflect_refract(transmitted, localV, H, VdotH, _cache.LdotH, rcpOrientedEta);

    return LightSample<IncomingRayDirInfo>::createTangentSpace(localV, IncomingRayDirInfo::create(localL), m);
}

template<class IncomingRayDirInfo>
    NBL_REQUIRES(ray_dir_info::basic<IncomingRayDirInfo>)
LightSample<IncomingRayDirInfo> ggx_cos_generate(
    NBL_CONST_REF_ARG(surface_interactions::Anisotropic<IncomingRayDirInfo>) interaction,
    NBL_REF_ARG(vector<typename IncomingRayDirInfo::scalar_t, 3>) u,
    typename IncomingRayDirInfo::scalar_t ax,
    typename IncomingRayDirInfo::scalar_t ay,
    typename IncomingRayDirInfo::scalar_t eta,
    NBL_REF_ARG(AnisotropicMicrofacetCache) _cache)
{
    using scalar_t = typename IncomingRayDirInfo::scalar_t;

    const vector<scalar_t, 3> localV = interaction.getTangentSpaceV();

    scalar_t orientedEta, rcpOrientedEta;
    const bool backside = math::getOrientedEtas(orientedEta, rcpOrientedEta, interaction.NdotV, eta);

    const vector<scalar_t, 3> upperHemisphereV = backside ? (-localV) : localV;

    const matrix<scalar_t, 3, 3> m = interaction.getTangentFrame();
    return ggx_cos_generate_wo_clamps<IncomingRayDirInfo>(localV, backside, upperHemisphereV, m, u, ax, ay, rcpOrientedEta, orientedEta * orientedEta, rcpOrientedEta * rcpOrientedEta, _cache);
}



template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar ggx_pdf_wo_clamps(
    bool transmitted, Scalar reflectance, Scalar ndf, Scalar devsh_v,
    Scalar absNdotV, Scalar VdotH, Scalar LdotH, Scalar VdotHLdotH, Scalar orientedEta)
{
    return geom_smith::VNDF_pdf_wo_clamps(ndf, geom_smith::ggx::G1_wo_numerator(absNdotV, devsh_v), absNdotV, transmitted, VdotH, LdotH, VdotHLdotH, orientedEta, reflectance);
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar ggx_pdf_wo_clamps(
    bool transmitted, Scalar reflectance, Scalar NdotH2, Scalar absNdotV,
    Scalar NdotV2, Scalar VdotH, Scalar LdotH, Scalar VdotHLdotH, Scalar a2, Scalar orientedEta)
{
    const Scalar ndf = ndf::ggx::trowbridge_reitz(a2, NdotH2);
    const Scalar devsh_v = geom_smith::ggx::devsh_part(NdotV2, a2, 1.0 - a2);

    return ggx_pdf_wo_clamps(transmitted, reflectance, ndf, devsh_v, absNdotV, VdotH, LdotH, VdotHLdotH, orientedEta);
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar ggx_pdf_wo_clamps(
    bool transmitted, Scalar reflectance, Scalar NdotH2, Scalar TdotH2,
    Scalar BdotH2, Scalar absNdotV, Scalar NdotV2, Scalar TdotV2, Scalar BdotV2,
    Scalar VdotH, Scalar LdotH, Scalar VdotHLdotH, Scalar ax, Scalar ay, Scalar ax2, Scalar ay2, Scalar orientedEta)
{
    const Scalar ndf = ndf::ggx::aniso(TdotH2, BdotH2, NdotH2, ax, ay, ax2, ay2);
    const Scalar devsh_v = geom_smith::ggx::devsh_part(TdotV2, BdotV2, NdotV2, ax2, ay2);

    return ggx_pdf_wo_clamps(transmitted, reflectance, ndf, devsh_v, absNdotV, VdotH, LdotH, VdotHLdotH, orientedEta);
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar ggx_cos_quotient_and_pdf_wo_clamps(
    NBL_REF_ARG(Scalar) pdf, Scalar ndf, bool transmitted, Scalar absNdotL,
    Scalar NdotL2, Scalar absNdotV, Scalar NdotV2, Scalar VdotH,
    Scalar LdotH, Scalar VdotHLdotH, Scalar reflectance, Scalar orientedEta, Scalar a2)
{
    const Scalar one_minus_a2 = 1.0 - a2;
    const Scalar devsh_v = geom_smith::ggx::devsh_part(NdotV2, a2, one_minus_a2);
    pdf = ggx_pdf_wo_clamps(transmitted, reflectance, ndf, devsh_v, absNdotV, VdotH, LdotH, VdotHLdotH, orientedEta);

    return geom_smith::ggx::G2_over_G1_devsh(absNdotL, NdotL2, absNdotV, devsh_v, a2, one_minus_a2);
}

template<class IncomingRayDirInfo>
    NBL_REQUIRES(ray_dir_info::basic<IncomingRayDirInfo>)
typename IncomingRayDirInfo::scalar_t ggx_cos_quotient_and_pdf(
    NBL_REF_ARG(typename IncomingRayDirInfo::scalar_t) pdf,
    NBL_CONST_REF_ARG(LightSample<IncomingRayDirInfo>) _sample,
    NBL_CONST_REF_ARG(surface_interactions::Isotropic<IncomingRayDirInfo>) interaction,
    NBL_CONST_REF_ARG(IsotropicMicrofacetCache) _cache,
    typename IncomingRayDirInfo::scalar_t eta,
    typename IncomingRayDirInfo::scalar_t a2)
{
    using scalar_t = typename IncomingRayDirInfo::scalar_t;
    
    const scalar_t ndf = ndf::ggx::trowbridge_reitz(a2, _cache.NdotH2);

    scalar_t orientedEta, dummy;
    const bool backside = math::getOrientedEtas(orientedEta, dummy, _cache.VdotH, eta);
    const scalar_t orientedEta2 = orientedEta * orientedEta;

    const scalar_t VdotHLdotH = _cache.VdotH * _cache.LdotH;
    const bool transmitted = VdotHLdotH < 0.0;

    const scalar_t reflectance = fresnel::dielectric_common(orientedEta2, abs(_cache.VdotH));

    const scalar_t absNdotV = abs(interaction.NdotV);
    return ggx_cos_quotient_and_pdf_wo_clamps(pdf, ndf, transmitted, abs(_sample.NdotL), _sample.NdotL2, absNdotV, interaction.NdotV_squared, _cache.VdotH, _cache.LdotH, VdotHLdotH, reflectance, orientedEta, a2);
}

template<class IncomingRayDirInfo>
    NBL_REQUIRES(ray_dir_info::basic<IncomingRayDirInfo>)
typename IncomingRayDirInfo::scalar_t ggx_aniso_cos_quotient_and_pdf_wo_clamps(
    NBL_REF_ARG(typename IncomingRayDirInfo::scalar_t) pdf,
    typename IncomingRayDirInfo::scalar_t ndf,
    bool transmitted,
    typename IncomingRayDirInfo::scalar_t absNdotL,
    typename IncomingRayDirInfo::scalar_t NdotL2,
    typename IncomingRayDirInfo::scalar_t TdotL2,
    typename IncomingRayDirInfo::scalar_t BdotL2,
    typename IncomingRayDirInfo::scalar_t absNdotV,
    typename IncomingRayDirInfo::scalar_t TdotV2,
    typename IncomingRayDirInfo::scalar_t BdotV2,
    typename IncomingRayDirInfo::scalar_t NdotV2,
    typename IncomingRayDirInfo::scalar_t VdotH,
    typename IncomingRayDirInfo::scalar_t LdotH,
    typename IncomingRayDirInfo::scalar_t VdotHLdotH,
    typename IncomingRayDirInfo::scalar_t reflectance,
    typename IncomingRayDirInfo::scalar_t orientedEta,
    typename IncomingRayDirInfo::scalar_t ax2,
    typename IncomingRayDirInfo::scalar_t ay2)
{
    const typename IncomingRayDirInfo::scalar_t devsh_v = geom_smith::ggx::devsh_part(TdotV2, BdotV2, NdotV2, ax2, ay2);
    pdf = ggx_pdf_wo_clamps(transmitted, reflectance, ndf, devsh_v, absNdotV, VdotH, LdotH, VdotHLdotH, orientedEta);

    return geom_smith::ggx::G2_over_G1_devsh(
        absNdotL, TdotL2, BdotL2, NdotL2,
        absNdotV, devsh_v,
        ax2, ay2
    );
}

template<class IncomingRayDirInfo>
    NBL_REQUIRES(ray_dir_info::basic<IncomingRayDirInfo>)
typename IncomingRayDirInfo::scalar_t ggx_aniso_cos_quotient_and_pdf(
    NBL_REF_ARG(typename IncomingRayDirInfo::scalar_t) pdf,
    NBL_CONST_REF_ARG(LightSample<IncomingRayDirInfo>) _sample,
    NBL_CONST_REF_ARG(surface_interactions::Anisotropic<IncomingRayDirInfo>) interaction,
    NBL_CONST_REF_ARG(AnisotropicMicrofacetCache) _cache,
    typename IncomingRayDirInfo::scalar_t eta,
    typename IncomingRayDirInfo::scalar_t ax,
    typename IncomingRayDirInfo::scalar_t ay)
{
    using scalar_t = typename IncomingRayDirInfo::scalar_t;

    const scalar_t ax2 = ax * ax;
    const scalar_t ay2 = ay * ay;
    const scalar_t TdotH2 = _cache.TdotH * _cache.TdotH;
    const scalar_t BdotH2 = _cache.BdotH * _cache.BdotH;
    const scalar_t ndf = ndf::ggx::aniso(TdotH2, BdotH2, _cache.NdotH2, ax, ay, ax2, ay2);

    const scalar_t TdotL2 = _sample.TdotL * _sample.TdotL;
    const scalar_t BdotL2 = _sample.BdotL * _sample.BdotL;

    const scalar_t TdotV2 = interaction.TdotV * interaction.TdotV;
    const scalar_t BdotV2 = interaction.BdotV * interaction.BdotV;

    const scalar_t VdotH = _cache.VdotH;

    scalar_t orientedEta, dummy;
    const bool backside = math::getOrientedEtas(orientedEta, dummy, VdotH, eta);
    const scalar_t orientedEta2 = orientedEta * orientedEta;

    const scalar_t VdotHLdotH = VdotH * _cache.LdotH;
    const bool transmitted = VdotHLdotH < 0.0;

    const scalar_t reflectance = fresnel::dielectric_common(orientedEta2, abs(VdotH));

    const scalar_t absNdotV = abs(interaction.NdotV);
    return ggx_aniso_cos_quotient_and_pdf_wo_clamps(pdf, ndf, transmitted, abs(_sample.NdotL), _sample.NdotL2, TdotL2, BdotL2, absNdotV, TdotV2, BdotV2, interaction.NdotV_squared, VdotH, _cache.LdotH, VdotHLdotH, reflectance, orientedEta, ax2, ay2);
}

}
}
}
}
}

#endif