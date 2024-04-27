// Copyright (C) 2018-2024 - DevSH Graphics Programming Sp. z O.O.
// This file is part of the "Nabla Engine".
// For conditions of distribution and use, see copyright notice in nabla.h

#ifndef _NBL_BUILTIN_HLSL_BXDF_BSDF_SPECULAR_BECKMANN_INCLUDED_
#define _NBL_BUILTIN_HLSL_BXDF_BSDF_SPECULAR_BECKMANN_INCLUDED_

#include <nbl/builtin/hlsl/bxdf/ndf/beckmann.hlsl>
#include <nbl/builtin/hlsl/bxdf/geom/smith/beckmann.hlsl>
#include <nbl/builtin/hlsl/bxdf/brdf/specular/beckmann.hlsl>
#include <nbl/builtin/hlsl/bxdf/bsdf/specular/common.hlsl>

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

template<class IncomingRayDirInfo>
    NBL_REQUIRES(ray_dir_info::basic<IncomingRayDirInfo>)
LightSample<IncomingRayDirInfo> beckmann_cos_generate_wo_clamps(
    NBL_CONST_REF_ARG(vector<typename IncomingRayDirInfo::scalar_t, 3>) localV,
    bool backside,
    NBL_CONST_REF_ARG(vector<typename IncomingRayDirInfo::scalar_t, 3>) upperHemisphereLocalV,
    NBL_CONST_REF_ARG(matrix<typename IncomingRayDirInfo::scalar_t, 3, 3>) m,
    NBL_CONST_REF_ARG(vector<typename IncomingRayDirInfo::scalar_t, 3>) u,
    typename IncomingRayDirInfo::scalar_t ax,
    typename IncomingRayDirInfo::scalar_t ay,
    typename IncomingRayDirInfo::scalar_t rcpOrientedEta,
    typename IncomingRayDirInfo::scalar_t orientedEta2,
    typename IncomingRayDirInfo::scalar_t rcpOrientedEta2,
    NBL_REF_ARG(AnisotropicMicrofacetCache<typename IncomingRayDirInfo::scalar_t>) _cache)
{
    using scalar_t = typename IncomingRayDirInfo::scalar_t;
    // thanks to this manoeuvre the H will always be in the upper hemisphere (NdotH>0.0)
    const vector<scalar_t, 3> H = brdf::specular::beckmann_cos_generate_wo_clamps(upperHemisphereLocalV,u.xy,ax,ay);

    const scalar_t VdotH = dot(localV,H);
    const scalar_t reflectance = fresnel::dielectric_common(orientedEta2,abs(VdotH));
    
    scalar_t rcpChoiceProb;
    const bool transmitted = math::partitionRandVariable(reflectance, u.z, rcpChoiceProb);
    
    vector<scalar_t, 3> localL;
    _cache = AnisotropicMicrofacetCache::create(localV, H, transmitted, rcpOrientedEta, rcpOrientedEta2);
    localL = math::reflect_refract(transmitted, localV, H, VdotH, _cache.LdotH, rcpOrientedEta);
    
    return LightSample<IncomingRayDirInfo>::createTangentSpace(localV, IncomingRayDirInfo::create(localL), m);
}

template<class IncomingRayDirInfo>
    NBL_REQUIRES(ray_dir_info::basic<IncomingRayDirInfo>)
LightSample<IncomingRayDirInfo> beckmann_cos_generate(
    NBL_CONST_REF_ARG(surface_interactions::Anisotropic<IncomingRayDirInfo>) interaction,
    NBL_REF_ARG(vector<typename IncomingRayDirInfo::scalar_t, 3>) u,
    typename IncomingRayDirInfo::scalar_t ax,
    typename IncomingRayDirInfo::scalar_t ay,
    typename IncomingRayDirInfo::scalar_t eta,
    NBL_REF_ARG(AnisotropicMicrofacetCache<typename IncomingRayDirInfo::scalar_t>) _cache)
{
    using scalar_t = typename IncomingRayDirInfo::scalar_t;

    const vector<scalar_t, 3> localV = interaction.getTangentSpaceV();
    
    float orientedEta, rcpOrientedEta;
    const bool backside = getOrientedEtas(orientedEta, rcpOrientedEta, interaction.NdotV, eta);
    
    const vector<scalar_t, 3> upperHemisphereV = backside ? (-localV):localV;

    const matrix<scalar_t, 3, 3> m = interaction.getTangentFrame();
    return beckmann_cos_generate_wo_clamps<IncomingRayDirInfo>(localV,backside,upperHemisphereV,m, u,ax,ay, rcpOrientedEta,orientedEta*orientedEta,rcpOrientedEta*rcpOrientedEta,_cache);
}


// isotropic PDF
template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar beckmann_pdf_wo_clamps(
    bool transmitted, Scalar reflectance, Scalar ndf, Scalar absNdotV, Scalar NdotV2, Scalar VdotH,
    Scalar LdotH, Scalar VdotHLdotH, Scalar a2, Scalar orientedEta, NBL_REF_ARG(Scalar) onePlusLambda_V)
{
    const Scalar lambda = geom_smith::beckmann::Lambda(NdotV2, a2);
    return geom_smith::VNDF_pdf_wo_clamps(ndf,lambda,absNdotV,transmitted,VdotH,LdotH,VdotHLdotH,orientedEta,reflectance,onePlusLambda_V);
}

// anisotropic PDF
template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar beckmann_pdf_wo_clamps(
    bool transmitted, Scalar reflectance, Scalar ndf, Scalar absNdotV, Scalar TdotV2, Scalar BdotV2,
    Scalar NdotV2, Scalar VdotH, Scalar LdotH, Scalar VdotHLdotH, Scalar ax2, Scalar ay2,
    Scalar orientedEta, NBL_REF_ARG(Scalar) onePlusLambda_V)
{
    Scalar c2 = geom_smith::beckmann::C2(TdotV2, BdotV2, NdotV2, ax2, ay2);
    Scalar lambda = geom_smith::beckmann::Lambda(c2);
    return geom_smith::VNDF_pdf_wo_clamps(ndf,lambda,absNdotV,transmitted,VdotH,LdotH,VdotHLdotH,orientedEta,reflectance,onePlusLambda_V);
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
quotient_and_pdf<Scalar, Scalar> beckmann_cos_quotient_and_pdf_wo_clamps(
    Scalar ndf, bool transmitted, Scalar NdotL2, Scalar absNdotV, Scalar NdotV2, Scalar VdotH,
    Scalar LdotH, Scalar VdotHLdotH, Scalar reflectance, Scalar orientedEta, Scalar a2)
{
    Scalar onePlusLambda_V;
    const Scalar pdf = beckmann_pdf_wo_clamps(transmitted, reflectance, ndf, absNdotV, NdotV2, VdotH, LdotH, VdotHLdotH, a2, orientedEta, onePlusLambda_V);

    return quotient_and_pdf<Scalar, Scalar>::create( geom_smith::beckmann::G2_over_G1(onePlusLambda_V, NdotL2, a2), pdf );
}

template<class IncomingRayDirInfo>
    NBL_REQUIRES(ray_dir_info::basic<IncomingRayDirInfo>)
quotient_and_pdf<typename IncomingRayDirInfo::scalar_t, typename IncomingRayDirInfo::scalar_t> beckmann_cos_quotient_and_pdf(
    NBL_CONST_REF_ARG(LightSample<IncomingRayDirInfo) _sample,
    NBL_CONST_REF_ARG(surface_interactions::Isotropic<IncomingRayDirInfo>) interaction,
    NBL_CONST_REF_ARG(IsotropicMicrofacetCache<typename IncomingRayDirInfo::scalar_t>) _cache,
    typename IncomingRayDirInfo::scalar_t eta,
    typename IncomingRayDirInfo::scalar_tt a2)
{
    using scalar_t = typename IncomingRayDirInfo::scalar_t;
    
    const scalar_t ndf = ndf::beckmann(a2,_cache.NdotH2);
    
    scalar_t orientedEta, dummy;
    const bool backside = math::getOrientedEtas(orientedEta, dummy, _cache.VdotH, eta);
    const scalar_t orientedEta2 = orientedEta*orientedEta;

    const scalar_t VdotHLdotH = _cache.VdotH*_cache.LdotH;
    const bool transmitted = VdotHLdotH<0.0;
    
    const scalar_t reflectance = fresnel::dielectric_common(orientedEta2,abs(_cache.VdotH));

    const scalar_t absNdotV = abs(interaction.NdotV);

    return beckmann_cos_quotient_and_pdf_wo_clamps(ndf, transmitted, _sample.NdotL2, absNdotV, interaction.NdotV_squared, _cache.VdotH, _cache.LdotH, VdotHLdotH, reflectance, orientedEta, a2);
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
quotient_and_pdf<Scalar, Scalar> beckmann_aniso_dielectric_cos_quotient_and_pdf_wo_clamps(
    Scalar ndf, in bool transmitted, Scalar NdotL2, Scalar TdotL2, Scalar BdotL2, Scalar absNdotV,
    Scalar TdotV2, Scalar BdotV2, Scalar NdotV2, Scalar VdotH, Scalar LdotH, Scalar VdotHLdotH,
    Scalar reflectance, Scalar orientedEta, Scalar ax2, Scalar ay2)
{
    Scalar onePlusLambda_V;
    const Scalar pdf = beckmann_pdf_wo_clamps(transmitted,reflectance, ndf,absNdotV,TdotV2,BdotV2,NdotV2, VdotH,LdotH,VdotHLdotH, ax2,ay2,orientedEta,onePlusLambda_V);

    return quotient_and_pdf<Scalar, Scalar>::create( geom_smith::beckmann::G2_over_G1(onePlusLambda_V, TdotL2, BdotL2, NdotL2, ax2, ay2), pdf );
}

template<class IncomingRayDirInfo>
    NBL_REQUIRES(ray_dir_info::basic<IncomingRayDirInfo>)
quotient_and_pdf<typename IncomingRayDirInfo::scalar_t, typename IncomingRayDirInfo::scalar_t> beckmann_aniso_dielectric_cos_quotient_and_pdf(
    NBL_CONST_REF_ARG(LightSample<IncomingRayDirInfo>) _sample,
    NBL_CONST_REF_ARG(surface_interactions::Anisotropic<IncomingRayDirInfo>) interaction,
    NBL_CONST_REF_ARG(AnisotropicMicrofacetCache<typename IncomingRayDirInfo::scalar_t>) _cache,
    typename IncomingRayDirInfo::scalar_t eta,
    typename IncomingRayDirInfo::scalar_t ax,
    typename IncomingRayDirInfo::scalar_t ay)
{
    using scalar_t = typename IncomingRayDirInfo::scalar_t;

    const scalar_t ax2 = ax*ax;
    const scalar_t ay2 = ay*ay;
    const scalar_t TdotH2 = _cache.TdotH*_cache.TdotH;
    const scalar_t BdotH2 = _cache.BdotH*_cache.BdotH;
    const scalar_t ndf = ndf::beckmann(ax,ay,ax2,ay2, TdotH2,BdotH2,_cache.NdotH2);

    const scalar_t TdotL2 = _sample.TdotL*_sample.TdotL;
    const scalar_t BdotL2 = _sample.BdotL*_sample.BdotL;
    
    const scalar_t TdotV2 = interaction.TdotV*interaction.TdotV;
    const scalar_t BdotV2 = interaction.BdotV*interaction.BdotV;
    
    const scalar_t VdotH = _cache.VdotH;

    scalar_t orientedEta, dummy;
    const bool backside = math::getOrientedEtas(orientedEta, dummy, VdotH, eta);
    const scalar_t orientedEta2 = orientedEta*orientedEta;
    
    const scalar_t VdotHLdotH = VdotH*_cache.LdotH;
    const bool transmitted = VdotHLdotH<0.0;
    
    const scalar_t reflectance = fresnel::dielectric_common(orientedEta2,abs(VdotH));
	return beckmann_aniso_dielectric_cos_quotient_and_pdf_wo_clamps(ndf, transmitted, _sample.NdotL2,TdotL2,BdotL2, abs(interaction.NdotV),TdotV2,BdotV2, interaction.NdotV_squared, VdotH,_cache.LdotH,VdotHLdotH, reflectance,orientedEta, ax2,ay2);
}


template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar beckmann_smith_height_correlated_dielectric_cos_eval_wo_clamps(
    Scalar NdotH2, Scalar NdotL2, Scalar absNdotV, Scalar NdotV2,
    bool transmitted, Scalar VdotH, Scalar LdotH, Scalar VdotHLdotH,
    Scalar orientedEta, Scalar orientedEta2, Scalar a2)
{
    const Scalar scalar_part = brdf::specular::beckmann_height_correlated_cos_eval_DG_wo_clamps(NdotH2, NdotL2, NdotV2, a2);
    
    const Scalar reflectance = fresnel::dielectric_common(orientedEta2,abs(VdotH));

    return reflectance*ndf::microfacet_to_light_measure_transform(scalar_part,absNdotV,transmitted,VdotH,LdotH,VdotHLdotH,orientedEta);
}

// before calling you must ensure that `AnisotropicMicrofacetCache` is valid (if a given V vector can "see" the L vector)
template<class IncomingRayDirInfo>
    NBL_REQUIRES(ray_dir_info::basic<IncomingRayDirInfo>)
typename IncomingRayDirInfo::scalar_t beckmann_smith_height_correlated_dielectric_cos_eval_wo_cache_validation(
    NBL_CONST_REF_ARG(LightSample<IncomingRayDirInfo>) _sample,
    NBL_CONST_REF_ARG(surface_interactions::Isotropic<IncomingRayDirInfo>) interaction,
    NBL_CONST_REF_ARG(IsotropicMicrofacetCache<typename IncomingRayDirInfo::scalar_t>) _cache,
    typename IncomingRayDirInfo::scalar_t eta,
    typename IncomingRayDirInfo::scalar_t a2)
{
    using scalar_t = typename IncomingRayDirInfo::scalar_t;

    scalar_t orientedEta, dummy;
    const bool backside = math::getOrientedEtas(orientedEta, dummy, _cache.VdotH, eta);
    const scalar_t orientedEta2 = orientedEta*orientedEta;
    
    const scalar_t VdotHLdotH = _cache.VdotH*_cache.LdotH;
    const bool transmitted = VdotHLdotH<0.0;

    return beckmann_smith_height_correlated_dielectric_cos_eval_wo_clamps(
        _cache.NdotH2,_sample.NdotL2,abs(interaction.NdotV),interaction.NdotV_squared,
        transmitted,_cache.VdotH,_cache.LdotH,VdotHLdotH,
        orientedEta,orientedEta2,a2);
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar beckmann_aniso_smith_height_correlated_dielectric_cos_eval_wo_clamps(
    Scalar NdotH2, Scalar TdotH2, Scalar BdotH2,
    Scalar NdotL2, Scalar TdotL2, Scalar BdotL2,
    Scalar absNdotV, Scalar NdotV2, Scalar TdotV2, Scalar BdotV2,
    bool transmitted, Scalar VdotH, Scalar LdotH, Scalar VdotHLdotH,
    Scalar orientedEta, Scalar orientedEta2,
    Scalar ax, Scalar ax2, Scalar ay, Scalar ay2)
{
    const Scalar scalar_part = brdf::specular::beckmann_aniso_height_correlated_cos_eval_DG_wo_clamps(NdotH2,TdotH2,BdotH2, NdotL2,TdotL2,BdotL2, NdotV2,TdotV2,BdotV2, ax, ax2, ay, ay2);
    
    const Scalar reflectance = fresnel::dielectric_common(orientedEta2,abs(VdotH));
    
    return reflectance*ndf::microfacet_to_light_measure_transform(scalar_part,absNdotV,transmitted,VdotH,LdotH,VdotHLdotH,orientedEta);
}

// before calling you must ensure that `AnisotropicMicrofacetCache` is valid (if a given V vector can "see" the L vector)
template<class IncomingRayDirInfo>
    NBL_REQUIRES(ray_dir_info::basic<IncomingRayDirInfo>)
typename IncomingRayDirInfo::scalar_t beckmann_aniso_smith_height_correlated_cos_eval_wo_cache_validation(
    NBL_CONST_REF_ARG(LightSample<IncomingRayDirInfo>) _sample,
    NBL_CONST_REF_ARG(surface_interactions::Anisotropic<IncomingRayDirInfo>) interaction,
    NBL_CONST_REF_ARG(AnisotropicMicrofacetCache<typename IncomingRayDirInfo::scalar_t>) _cache,
    typename IncomingRayDirInfo::scalar_t eta
    typename IncomingRayDirInfo::scalar_t ax,
    typename IncomingRayDirInfo::scalar_t ax2,
    typename IncomingRayDirInfo::scalar_t ay,
    typename IncomingRayDirInfo::scalar_t ay2)
{
    using scalar_t = typename IncomingRayDirInfo::scalar_t;

    const scalar_t TdotH2 = _cache.TdotH*_cache.TdotH;
    const scalar_t BdotH2 = _cache.BdotH*_cache.BdotH;

    const scalar_t TdotL2 = _sample.TdotL*_sample.TdotL;
    const scalar_t BdotL2 = _sample.BdotL*_sample.BdotL;

    const scalar_t TdotV2 = interaction.TdotV*interaction.TdotV;
    const scalar_t BdotV2 = interaction.BdotV*interaction.BdotV;

    const scalar_t VdotH = _cache.VdotH;

    scalar_t orientedEta, dummy;
    const bool backside = math::getOrientedEtas(orientedEta, dummy, VdotH, eta);
    const scalar_t orientedEta2 = orientedEta*orientedEta;
    
    const scalar_t VdotHLdotH = VdotH*_cache.LdotH;
    const bool transmitted = VdotHLdotH<0.0;

    return beckmann_aniso_smith_height_correlated_dielectric_cos_eval_wo_clamps(
        _cache.NdotH2,TdotH2,BdotH2,
        _sample.NdotL2,TdotL2,BdotL2,
        abs(interaction.NdotV),interaction.NdotV_squared,TdotV2,BdotV2,
        transmitted,VdotH,_cache.LdotH,VdotHLdotH,
        orientedEta,orientedEta2,ax,ax*ax,ay,ay*ay);
}

}
}
}
}
}

#endif