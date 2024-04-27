// Copyright (C) 2018-2024 - DevSH Graphics Programming Sp. z O.O.
// This file is part of the "Nabla Engine".
// For conditions of distribution and use, see copyright notice in nabla.h

#ifndef _NBL_BUILTIN_HLSL_BXDF_BRDF_SPECULAR_BLINN_PHONG_INCLUDED_
#define _NBL_BUILTIN_HLSL_BXDF_BRDF_SPECULAR_BLINN_PHONG_INCLUDED_

#include <nbl/builtin/hlsl/bxdf/common.hlsl>
#include <nbl/builtin/hlsl/bxdf/reflection.hlsl>
#include <nbl/builtin/hlsl/bxdf/fresnel.hlsl>
#include <nbl/builtin/hlsl/bxdf/ndf/blinn_phong.hlsl>
#include <nbl/builtin/hlsl/bxdf/geom/smith/beckmann.hlsl>

namespace nbl
{
namespace hlsl
{
namespace bxdf
{
namespace brdf
{
namespace specular
{

//conversion between alpha and Phong exponent, Walter et.al.
template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
float phong_exp_to_alpha2(Scalar n)
{
    return 2.0/(n+2.0);
}
//+INF for a2==0.0
template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
float alpha2_to_phong_exp(Scalar a2)
{
    return 2.0/a2 - 2.0;
}

//https://zhuanlan.zhihu.com/p/58205525
//only NDF sampling
//however we dont really care about phong sampling
template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
vector<Scalar, 3> blinn_phong_cos_generate(vector<Scalar, 2> u, Scalar n)
{
    Scalar phi = 2.0*math::PI*u.y;
    Scalar cosTheta = pow(u.x, 1.0/(n+1.0));
    Scalar sinTheta = sqrt(1.0 - cosTheta*cosTheta);
    Scalar cosPhi = cos(phi);
    Scalar sinPhi = sin(phi);
    return vector<Scalar, 3>(cosPhi*sinTheta, sinPhi*sinTheta, cosTheta);
}

template<class RayDirInfo>
    NBL_REQUIRES(ray_dir_info::basic<RayDirInfo>)
LightSample<RayDirInfo> blinn_phong_cos_generate(NBL_CONST_REF_ARG(
    NBL_CONST_REF_ARG(surface_interactions::Anisotropic<RayDirInfo>) interaction,
    NBL_CONST_REF_ARG(vector<typename RayDirInfo::scalar_t, 2>) u,
    typename RayDirInfo::scalar_t n,
    NBL_REF_ARG(AnisotropicMicrofacetCache<typename RayDirInfo::scalar_t>) _cache)
{
    using scalar_t = typename RayDirInfo::scalar_t;
    const vector<scalar_t, 3> H = blinn_phong_cos_generate(u,n);
    const vector<scalar_t, 3> localV = interaction.getTangentSpaceV();

    _cache = AnisotropicMicrofacetCache::create(localV, H);
    vector<scalar_t, 3> localL;
    localL = math::reflect(localV, H, _cache.VdotH);
    
    const matrix<scalar_t, 3, 3> m = interaction.getTangentFrame();

    return LightSample<RayDirInfo>::createTangentSpace(localV, RayDirInfo::create(localL), m);
}

/*
float3 blinn_phong_dielectric_cos_quotient_and_pdf(out float& pdf, in BxDFSample s, in surface_interactions::Isotropic<RayDirInfo> interaction, in float n, in float3 ior)
{
	pdf = (n+1.0)*0.5*RECIPROCAL_PI * 0.25*pow(s.NdotH,n)/s.VdotH;

    float3 fr = fresnel_dielectric(ior, s.VdotH);
    return fr * s.NdotL * (n*(n + 6.0) + 8.0) * s.VdotH / ((pow(0.5,0.5*n) + n) * (n + 1.0));
}

float3 blinn_phong_conductor_cos_quotient_and_pdf(out float& pdf, in BxDFSample s, in surface_interactions::Isotropic<RayDirInfo> interaction, in float n, in float2x3 ior)
{
	pdf = (n+1.0)*0.5*RECIPROCAL_PI * 0.25*pow(s.NdotH,n)/s.VdotH;

    float3 fr = fresnel_conductor(ior[0], ior[1], s.VdotH);
    return fr * s.NdotL * (n*(n + 6.0) + 8.0) * s.VdotH / ((pow(0.5,0.5*n) + n) * (n + 1.0));
}
*/

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar blinn_phong_cos_eval_DG_wo_clamps(Scalar NdotH, Scalar NdotV_squared, Scalar NdotL2, Scalar n, Scalar a2)
{
    Scalar NG = blinn_phong(NdotH, n);
    if (a2>FLT_MIN)
        NG *= geom_smith::beckmann::correlated(NdotV_squared, NdotL2, a2);
    return NG;
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar blinn_phong_cos_eval_DG_wo_clamps(Scalar NdotH, Scalar NdotV_squared, Scalar NdotL2, Scalar n)
{
    Scalar a2 = phong_exp_to_alpha2(n);
    return blinn_phong_cos_eval_DG_wo_clamps(NdotH, NdotV_squared, NdotL2, n, a2);
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
vector<Scalar, 3> blinn_phong_cos_eval_wo_clamps(Scalar NdotH, Scalar maxNdotV, Scalar NdotV_squared, Scalar NdotL2,
    Scalar VdotH, Scalar n, NBL_CONST_REF_ARG(matrix<Scalar, 2,3>) ior, Scalar a2)
{
    Scalar scalar_part = blinn_phong_cos_eval_DG_wo_clamps(NdotH, NdotV_squared, NdotL2, n, a2);
    return fresnel::conductor(ior[0], ior[1], VdotH)*ndf::microfacet_to_light_measure_transform(scalar_part,maxNdotV);
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
vector<Scalar, 3> blinn_phong_cos_eval_wo_clamps(Scalar NdotH, Scalar maxNdotV, Scalar NdotV_squared, Scalar NdotL2, Scalar VdotH, Scalar n, NBL_CONST_REF_ARG(matrix<Scalar, 2,3>) ior)
{
    Scalar a2 = phong_exp_to_alpha2(n);
    return blinn_phong_cos_eval_wo_clamps(NdotH, maxNdotV, NdotV_squared, NdotL2, VdotH, n, ior, a2);
}

template <class RayDirInfo>
vector<Scalar, 3> blinn_phong_cos_eval(
    NBL_CONST_REF_ARG(LightSample<RayDirInfo>) _sample,
    NBL_CONST_REF_ARG(surface_interactions::Isotropic<RayDirInfo>) interaction,
    NBL_CONST_REF_ARG(IsotropicMicrofacetCache<Scalar>) _cache,
    Scalar n,
    NBL_CONST_REF_ARG(matrix<Scalar, 2,3>) ior)
{
    if (interaction.NdotV>FLT_MIN)
        return blinn_phong_cos_eval_wo_clamps(_cache.NdotH, interaction.NdotV, interaction.NdotV_squared, _sample.NdotL2, _cache.VdotH, n, ior);
    else
        return vector<Scalar, 3>(0.0,0.0,0.0);
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar blinn_phong_cos_eval_DG_wo_clamps(Scalar NdotH, Scalar NdotH2, Scalar TdotH2, Scalar BdotH2,
    Scalar TdotL2, Scalar BdotL2, Scalar TdotV2, Scalar BdotV2, Scalar NdotV_squared, Scalar NdotL2, Scalar nx, Scalar ny, Scalar ax2, Scalar ay2)
{
    Scalar DG = blinn_phong(NdotH, 1.0/(1.0-NdotH2), TdotH2, BdotH2, nx, ny);
    if (ax2>FLT_MIN || ay2>FLT_MIN)
        DG *= geom_smith::beckmann::correlated(TdotV2, BdotV2, NdotV_squared, TdotL2, BdotL2, NdotL2, ax2, ay2);
    return DG;
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar blinn_phong_cos_eval_DG_wo_clamps(Scalar NdotH, Scalar NdotH2, Scalar TdotH2, Scalar BdotH2,
    Scalar TdotL2, Scalar BdotL2, Scalar TdotV2, Scalar BdotV2, Scalar NdotV_squared, Scalar NdotL2, Scalar nx, Scalar ny)
{
    Scalar ax2 = phong_exp_to_alpha2(nx);
    Scalar ay2 = phong_exp_to_alpha2(ny);

    return blinn_phong_cos_eval_DG_wo_clamps(NdotH, NdotH2, TdotH2, BdotH2, TdotL2, BdotL2, TdotV2, BdotV2, NdotV_squared, NdotL2, nx, ny, ax2, ay2);
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
vector<Scalar, 3> blinn_phong_cos_eval_wo_clamps(Scalar NdotH, Scalar NdotH2, Scalar TdotH2, Scalar BdotH2, Scalar TdotL2,
    Scalar BdotL2, Scalar maxNdotV, Scalar TdotV2, Scalar BdotV2, Scalar NdotV_squared, Scalar NdotL2, Scalar VdotH, Scalar nx,
    Scalar ny, NBL_CONST_REF_ARG(matrix<Scalar, 2,3>) ior, Scalar ax2, Scalar ay2)
{
    Scalar scalar_part = blinn_phong_cos_eval_DG_wo_clamps(NdotH, NdotH2, TdotH2, BdotH2, TdotL2, BdotL2, TdotV2, BdotV2, NdotV_squared, NdotL2, nx, ny, ax2, ay2);

    return fresnel::conductor(ior[0], ior[1], VdotH)*ndf::microfacet_to_light_measure_transform(scalar_part,maxNdotV);
}
template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
vector<Scalar, 3> blinn_phong_cos_eval_wo_clamps(Scalar NdotH, Scalar NdotH2, Scalar TdotH2, Scalar BdotH2, Scalar TdotL2,
    Scalar BdotL2, Scalar maxNdotV, Scalar TdotV2, Scalar BdotV2, Scalar NdotV_squared, Scalar NdotL2, Scalar VdotH,
    Scalar nx, Scalar ny, NBL_CONST_REF_ARG(matrix<Scalar, 2,3>) ior)
{
    Scalar ax2 = phong_exp_to_alpha2(nx);
    Scalar ay2 = phong_exp_to_alpha2(ny);

    return blinn_phong_cos_eval_wo_clamps(NdotH, NdotH2, TdotH2, BdotH2, TdotL2, BdotL2, maxNdotV, TdotV2, BdotV2, NdotV_squared, NdotL2, VdotH, nx, ny, ior, ax2, ay2);
}

template<class RayDirInfo>
    NBL_REQUIRES(ray_dir_info::basic<RayDirInfo>)
vector<typename RayDirInfo::scalar_t, 3> blinn_phong_cos_eval(
    NBL_CONST_REF_ARG(LightSample<RayDirInfo>) _sample,
    NBL_CONST_REF_ARG(surface_interactions::Anisotropic<RayDirInfo>)interaction,
    NBL_CONST_REF_ARG(AnisotropicMicrofacetCache<typename RayDirInfo::scalar_t>) _cache,
    NBL_CONST_REF_ARG(typename RayDirInfo::scalar_t) nx,
    NBL_CONST_REF_ARG(typename RayDirInfo::scalar_t) ny,
    NBL_CONST_REF_ARG(matrix<typename RayDirInfo::scalar_t, 2, 3>) ior)
{   
    using scalar_t = typename RayDirInfo::scalar_t;
    if (interaction.NdotV>FLT_MIN)
    {
        const scalar_t TdotH2 = _cache.TdotH*_cache.TdotH;
        const scalar_t BdotH2 = _cache.BdotH*_cache.BdotH;

        const scalar_t TdotL2 = _sample.TdotL*_sample.TdotL;
        const scalar_t BdotL2 = _sample.BdotL*_sample.BdotL;

        const scalar_t TdotV2 = interaction.TdotV*interaction.TdotV;
        const scalar_t BdotV2 = interaction.BdotV*interaction.BdotV;
        return blinn_phong_cos_eval_wo_clamps(_cache.NdotH, _cache.NdotH2, TdotH2, BdotH2, TdotL2, BdotL2, interaction.NdotV, TdotV2, BdotV2, interaction.NdotV_squared, _sample.NdotL2, _cache.VdotH, nx, ny, ior);
    }
    else
        return vector<scalar_t, 3>(0.0,0.0,0.0);
}

}
}
}
}
}

#endif