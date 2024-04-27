// Copyright (C) 2018-2024 - DevSH Graphics Programming Sp. z O.O.
// This file is part of the "Nabla Engine".
// For conditions of distribution and use, see copyright notice in nabla.h
#ifndef _NBL_BUILTIN_HLSL_BXDF_BRDF_SPECULAR_BECKMANN_INCLUDED_
#define _NBL_BUILTIN_HLSL_BXDF_BRDF_SPECULAR_BECKMANN_INCLUDED_

#include <nbl/builtin/hlsl/bxdf/reflection.hlsl>
#include <nbl/builtin/hlsl/bxdf/fresnel.hlsl>
#include <nbl/builtin/hlsl/bxdf/ndf/beckmann.hlsl>
#include <nbl/builtin/hlsl/bxdf/geom/smith/common.hlsl>
#include <nbl/builtin/hlsl/bxdf/geom/smith/beckmann.hlsl>
#include <nbl/builtin/hlsl/math/functions.hlsl>
#include <nbl/builtin/hlsl/bxdf/brdf/specular/common.hlsl>

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

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
vector<Scalar, 3> beckmann_cos_generate_wo_clamps(NBL_CONST_REF_ARG(vector<Scalar, 3>) localV, NBL_CONST_REF_ARG(vector<Scalar, 2>) u, Scalar ax, Scalar ay)
{
    //stretch
    vector<Scalar, 3> V = normalize(vector<Scalar, 3> (ax*localV.x, ay*localV.y, localV.z));

    vector<Scalar, 2> slope;
    if (V.z > 0.9999)//V.z=NdotV=cosTheta in tangent space
    {
        Scalar r = sqrt(-log(1.0-u.x));
        Scalar sinPhi = sin(2.0*math::PI*u.y);
        Scalar cosPhi = cos(2.0*math::PI*u.y);
        slope = vector<Scalar, 2> (r,r)*vector<Scalar, 2> (cosPhi,sinPhi);
    }
    else
    {
        Scalar cosTheta = V.z;
        Scalar sinTheta = sqrt(1.0 - cosTheta*cosTheta);
        Scalar tanTheta = sinTheta/cosTheta;
        Scalar cotTheta = 1.0/tanTheta;
        
        Scalar a = -1.0;
        Scalar c = math::erf(cosTheta);
        Scalar sample_x = max(u.x, 1.0e-6f);
        Scalar theta = acos(cosTheta);
        Scalar fit = 1.0 + theta * (-0.876 + theta * (0.4265 - 0.0594*theta));
        Scalar b = c - (1.0 + c) * pow(1.0-sample_x, fit);
        
        Scalar normalization = 1.0 / (1.0 + c + math::SQRT_RECIPROCAL_PI * tanTheta * exp(-cosTheta*cosTheta));

        const int ITER_THRESHOLD = 10;
		const Scalar MAX_ACCEPTABLE_ERR = 1.0e-5;
        int it = 0;
        Scalar value=1000.0;
        while (++it<ITER_THRESHOLD && abs(value)>MAX_ACCEPTABLE_ERR)
        {
            if (!(b>=a && b<=c))
                b = 0.5 * (a+c);

            Scalar invErf = math::erfInv(b);
            value = normalization * (1.0 + b + math::SQRT_RECIPROCAL_PI * tanTheta * exp(-invErf*invErf)) - sample_x;
            Scalar derivative = normalization * (1.0 - invErf*cosTheta);

            if (value > 0.0)
                c = b;
            else
                a = b;

            b -= value/derivative;
        }
        // TODO: investigate if we can replace these two erf^-1 calls with a box muller transform
        slope.x = math::erfInv(b);
        slope.y = math::erfInv(2.0f * max(u.y,1.0e-6f) - 1.0f);
    }
    
    Scalar sinTheta = sqrt(1.0f - V.z*V.z);
    Scalar cosPhi = sinTheta==0.0f ? 1.0f : clamp(V.x/sinTheta, -1.0f, 1.0f);
    Scalar sinPhi = sinTheta==0.0f ? 0.0f : clamp(V.y/sinTheta, -1.0f, 1.0f);
    //rotate
    Scalar tmp = cosPhi*slope.x - sinPhi*slope.y;
    slope.y = sinPhi*slope.x + cosPhi*slope.y;
    slope.x = tmp;

    //unstretch
    slope = vector<Scalar, 2> (ax,ay)*slope;

    return normalize(vector<Scalar, 3> (-slope, 1.0));
}

// TODO: unifty the two following functions into `microfacet_BRDF_cos_generate_wo_clamps(float3 H,...)` and `microfacet_BRDF_cos_generate` or at least a auto declaration macro in lieu of a template
template<class RayDirInfo>
    NBL_REQUIRES(ray_dir_info::basic<RayDirInfo>)
LightSample<RayDirInfo> beckmann_cos_generate_wo_clamps(
    NBL_CONST_REF_ARG(vector<typename RayDirInfo::scalar_t, 3>) localV,
    NBL_CONST_REF_ARG(matrix<typename RayDirInfo::scalar_t, 3, 3>) m,
    NBL_CONST_REF_ARG(vector<typename RayDirInfo::scalar_t, 2>) u,
    typename RayDirInfo::scalar_t ax,
    typename RayDirInfo::scalar_t ay,
    NBL_REF_ARG(AnisotropicMicrofacetCache<typename RayDirInfo::scalar_t>) _cache)
{
    const vector<typename RayDirInfo::scalar_t, 3> H = beckmann_cos_generate_wo_clamps(localV,u,ax,ay);
    
    _cache = AnisotropicMicrofacetCache::create(localV,H);
    vector<typename RayDirInfo::scalar_t, 3> localL = math::reflect(localV, H, _cache.VdotH);
    
    return LightSample<RayDirInfo>::createTangentSpace(localV, RayDirInfo::create(localL), m);
}

template<class RayDirInfo>
    NBL_REQUIRES(ray_dir_info::basic<RayDirInfo>)
LightSample<RayDirInfo> beckmann_cos_generate(
    NBL_CONST_REF_ARG(surface_interactions::Anisotropic<RayDirInfo>) interaction,
    NBL_CONST_REF_ARG(vector<typename RayDirInfo::scalar_t, 2>) u,
    typename RayDirInfo::scalar_t ax,
    typename RayDirInfo::scalar_t ay,
    NBL_REF_ARG(AnisotropicMicrofacetCache<typename RayDirInfo::scalar_t>) _cache)
{
    const vector<typename RayDirInfo::scalar_t, 3> localV = interaction.getTangentSpaceV();
    const matrix<typename RayDirInfo::scalar_t, 3, 3> m = interaction.getTangentFrame();
    return beckmann_cos_generate_wo_clamps<RayDirInfo>(localV,m,u,ax,ay,_cache);
}


// isotropic PDF
template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar beckmann_pdf_wo_clamps(Scalar ndf, Scalar maxNdotV, Scalar NdotV2, Scalar a2, NBL_REF_ARG(Scalar) onePlusLambda_V)
{
    const Scalar lambda = geom_smith::beckmann::Lambda(NdotV2, a2);
    return geom_smith::VNDF_pdf_wo_clamps(ndf,lambda,maxNdotV,onePlusLambda_V);
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar beckmann_pdf_wo_clamps(Scalar NdotH2, Scalar maxNdotV, Scalar NdotV2, Scalar a2)
{
    Scalar ndf = ndf::beckmann(a2, NdotH2);

    Scalar dummy;
    return beckmann_pdf_wo_clamps(ndf, maxNdotV,NdotV2, a2, dummy);
}

// anisotropic PDF
template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar beckmann_pdf_wo_clamps(Scalar ndf, Scalar maxNdotV, Scalar TdotV2, Scalar BdotV2, Scalar NdotV2, Scalar ax2, Scalar ay2, NBL_REF_ARG(Scalar) onePlusLambda_V)
{
    Scalar c2 = geom_smith::beckmann::C2(TdotV2, BdotV2, NdotV2, ax2, ay2);
    Scalar lambda = geom_smith::beckmann::Lambda(c2);

    return geom_smith::VNDF_pdf_wo_clamps(ndf, lambda, maxNdotV, onePlusLambda_V);
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar beckmann_pdf_wo_clamps(Scalar NdotH2, Scalar TdotH2, Scalar BdotH2, Scalar maxNdotV, Scalar TdotV2, Scalar BdotV2, Scalar NdotV2, Scalar ax, Scalar ax2, Scalar ay, Scalar ay2)
{
    Scalar ndf = ndf::beckmann(ax, ay, ax2, ay2, TdotH2, BdotH2, NdotH2);

    Scalar dummy;
    return beckmann_pdf_wo_clamps(ndf, maxNdotV, TdotV2, BdotV2, NdotV2, ax2, ay2, dummy);
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
quotient_and_pdf<vector<Scalar, 3>, Scalar> beckmann_cos_quotient_and_pdf_wo_clamps(Scalar ndf, Scalar NdotL2, Scalar maxNdotV, Scalar NdotV2,
    NBL_CONST_REF_ARG(vector<Scalar, 3>) reflectance, Scalar a2)
{
    Scalar onePlusLambda_V;
    Scalar pdf = beckmann_pdf_wo_clamps(ndf,maxNdotV,NdotV2,a2,onePlusLambda_V);

    Scalar G2_over_G1 = geom_smith::beckmann::G2_over_G1(onePlusLambda_V, NdotL2, a2);
    return quotient_and_pdf<vector<Scalar, 3>, Scalar>::create(reflectance*G2_over_G1, pdf);
}

template<class RayDirInfo>
    NBL_REQUIRES(ray_dir_info::basic<RayDirInfo>)
quotient_and_pdf<vector<typename RayDirInfo::scalar_t, 3>, typename RayDirInfo::scalar_t> beckmann_cos_quotient_and_pdf(
    NBL_CONST_REF_ARG(LightSample<RayDirInfo>) _sample,
    NBL_CONST_REF_ARG(surface_interactions::Isotropic<RayDirInfo>) interaction,
    NBL_CONST_REF_ARG(IsotropicMicrofacetCache<typename RayDirInfo::scalar_t>) _cache,
    NBL_CONST_REF_ARG(matrix<typename RayDirInfo::scalar_t, 2, 3>) ior,
    NBL_CONST_REF_ARG(typename RayDirInfo::scalar_t)  a2)
{
    using scalar_t = typename RayDirInfo::scalar_t;
    const scalar_t ndf = ndf::beckmann(a2, _cache.NdotH2);
    scalar_t onePlusLambda_V;
    scalar_t pdf = beckmann_pdf_wo_clamps(ndf, interaction.NdotV, interaction.NdotV2, a2, onePlusLambda_V);
    vector<scalar_t, 3> rem = vector<scalar_t, 3>(0.0,0.0,0.0);
    if (_sample.NdotL>FLT_MIN && interaction.NdotV>FLT_MIN)
    {
        const vector<scalar_t, 3> reflectance = fresnel::conductor(ior[0], ior[1], _cache.VdotH);
    
        scalar_t G2_over_G1 = beckmann_smith_G2_over_G1(onePlusLambda_V, _sample.NdotL2, a2);
        rem = reflectance * G2_over_G1;
    }
    
    return quotient_and_pdf<vector<typename RayDirInfo::scalar_t, 3>, typename RayDirInfo::scalar_t>::create(rem, pdf);
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
quotient_and_pdf<vector<Scalar, 3>, Scalar> beckmann_aniso_cos_quotient_and_pdf_wo_clamps(Scalar ndf, Scalar NdotL2, Scalar TdotL2, Scalar BdotL2, Scalar maxNdotV,
    Scalar TdotV2, Scalar BdotV2, Scalar NdotV2, NBL_CONST_REF_ARG(vector<Scalar, 3>) reflectance, Scalar ax2, Scalar ay2)
{
    Scalar onePlusLambda_V;
    Scalar pdf = beckmann_pdf_wo_clamps(ndf,maxNdotV,TdotV2,BdotV2,NdotV2,ax2,ay2,onePlusLambda_V);

    Scalar G2_over_G1 = geom_smith::beckmann::G2_over_G1(onePlusLambda_V, TdotL2, BdotL2, NdotL2, ax2, ay2);
    return quotient_and_pdf<vector<Scalar, 3>, Scalar>::create(reflectance * G2_over_G1, pdf);
}

template<class RayDirInfo>
    NBL_REQUIRES(ray_dir_info::basic<RayDirInfo>)
quotient_and_pdf<vector<typename RayDirInfo::scalar_t, 3>, typename RayDirInfo::scalar_t> beckmann_aniso_cos_quotient_and_pdf(
    NBL_CONST_REF_ARG(LightSample<RayDirInfo>) _sample,
    NBL_CONST_REF_ARG(surface_interactions::Anisotropic<RayDirInfo>) interaction,
    NBL_CONST_REF_ARG(AnisotropicMicrofacetCache<typename RayDirInfo::scalar_t>) _cache,
    NBL_CONST_REF_ARG(matrix<typename RayDirInfo::scalar_t, 2, 3>) ior,
    NBL_CONST_REF_ARG(typename RayDirInfo::scalar_t) ax,
    NBL_CONST_REF_ARG(typename RayDirInfo::scalar_t) ay)
{
    using scalar_t = typename RayDirInfo::scalar_t;
    const scalar_t ax2 = ax * ax;
    const scalar_t ay2 = ay * ay;

    const scalar_t TdotH2 = _cache.TdotH * _cache.TdotH;
    const scalar_t BdotH2 = _cache.BdotH * _cache.BdotH;
    const scalar_t TdotV2 = interaction.TdotV * interaction.TdotV;
    const scalar_t BdotV2 = interaction.BdotV * interaction.BdotV;

    const scalar_t NdotV2 = interaction.NdotV2;

    const scalar_t ndf = ndf::beckmann(ax, ay, ax2, ay2, TdotH2, BdotH2, _cache.NdotH2);
    scalar_t onePlusLambda_V;
    scalar_t pdf = beckmann_pdf_wo_clamps(ndf, interaction.NdotV, TdotV2, BdotV2, NdotV2, ax2, ay2, onePlusLambda_V);
    quotient_and_pdf<vector<typename RayDirInfo::scalar_t, 3>, typename RayDirInfo::scalar_t> qpdf =
        quotient_and_pdf<vector<typename RayDirInfo::scalar_t, 3>, typename RayDirInfo::scalar_t>::create(float3(0.0, 0.0, 0.0), pdf);
    if (_sample.NdotL>FLT_MIN && interaction.NdotV>FLT_MIN)
    {
        const scalar_t TdotL2 = _sample.TdotL*_sample.TdotL;
        const scalar_t BdotL2 = _sample.BdotL*_sample.BdotL;
    
        const vector<scalar_t, 3> reflectance = fresnel::conductor(ior[0], ior[1], _cache.VdotH);

        qpdf = beckmann_aniso_cos_quotient_and_pdf_wo_clamps(ndf, _sample.NdotL2, TdotL2, BdotL2, interaction.NdotV, TdotV2, BdotV2, NdotV2, reflectance, ax2, ay2);
    }
    
    return qpdf;
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar beckmann_height_correlated_cos_eval_DG_wo_clamps(Scalar NdotH2, Scalar NdotL2, Scalar NdotV2, Scalar a2)
{
    Scalar NG = ndf::beckmann(a2, NdotH2);
    if  (a2>FLT_MIN)
        NG *= geom_smith::beckmann::correlated(NdotV2, NdotL2, a2);
    
    return NG;
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
vector<Scalar, 3> beckmann_height_correlated_cos_eval_wo_clamps(Scalar NdotH2, Scalar NdotL2, Scalar maxNdotV, Scalar NdotV2, Scalar VdotH,
    NBL_CONST_REF_ARG(matrix<Scalar, 2, 3>) ior, Scalar a2)
{
    const Scalar NG = beckmann_height_correlated_cos_eval_DG_wo_clamps(NdotH2, NdotL2, NdotV2, a2);

    const vector<Scalar, 3> fr = fresnel::conductor(ior[0], ior[1], VdotH);

    return fr*ndf::microfacet_to_light_measure_transform(NG,maxNdotV);
}

template<class RayDirInfo>
    NBL_REQUIRES(ray_dir_info::basic<RayDirInfo>)
vector<typename RayDirInfo::scalar_t, 3> beckmann_height_correlated_cos_eval(
    NBL_CONST_REF_ARG(LightSample<RayDirInfo>) _sample,
    NBL_CONST_REF_ARG(surface_interactions::Isotropic<RayDirInfo>) interaction,
    NBL_CONST_REF_ARG(IsotropicMicrofacetCache<typename RayDirInfo::scalar_t>) _cache,
    NBL_CONST_REF_ARG(matrix<typename RayDirInfo::scalar_t, 2, 3>) ior,
    typename RayDirInfo::scalar_t a2)
{
    if (interaction.NdotV>FLT_MIN)
        return beckmann_height_correlated_cos_eval_wo_clamps(_cache.NdotH2,_sample.NdotL2,interaction.NdotV,interaction.NdotV2,_cache.VdotH,ior,a2);
    else
        return vector<typename RayDirInfo::scalar_t, 3>(0.0,0.0,0.0);
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar beckmann_aniso_height_correlated_cos_eval_DG_wo_clamps(Scalar NdotH2, Scalar TdotH2, Scalar BdotH2, Scalar NdotL2, Scalar TdotL2,
    Scalar BdotL2, Scalar NdotV2, Scalar TdotV2, Scalar BdotV2, Scalar ax, Scalar ax2, Scalar ay, Scalar ay2)
{
    Scalar NG = ndf::beckmann(ax, ay, ax2, ay2, TdotH2, BdotH2, NdotH2);
    if (ax>FLT_MIN || ay>FLT_MIN)
        NG *= geom_smith::beckmann::correlated(TdotV2, BdotV2, NdotV2, TdotL2, BdotL2, NdotL2, ax2, ay2);
    
    return NG;
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
vector<Scalar, 3> beckmann_aniso_height_correlated_cos_eval_wo_clamps(Scalar NdotH2, Scalar TdotH2, Scalar BdotH2, Scalar NdotL2, Scalar TdotL2,
    Scalar BdotL2, Scalar maxNdotV, Scalar NdotV2, Scalar TdotV2, Scalar BdotV2, Scalar VdotH, NBL_CONST_REF_ARG(matrix<Scalar, 2, 3>) ior,
    Scalar ax, Scalar ax2, Scalar ay, Scalar ay2)
{
    const Scalar NG = beckmann_aniso_height_correlated_cos_eval_DG_wo_clamps(NdotH2,TdotH2,BdotH2, NdotL2,TdotL2,BdotL2, NdotV2,TdotV2,BdotV2, ax, ax2, ay, ay2);

    const vector<Scalar, 3> fr = fresnel::conductor(ior[0], ior[1], VdotH);
    
    return fr*ndf::microfacet_to_light_measure_transform(NG,maxNdotV);
}

template<class RayDirInfo>
    NBL_REQUIRES(ray_dir_info::basic<RayDirInfo>)
vector<typename RayDirInfo::scalar_t, 3> beckmann_aniso_height_correlated_cos_eval(
    NBL_CONST_REF_ARG(LightSample<RayDirInfo>) _sample,
    NBL_CONST_REF_ARG(surface_interactions::Isotropic<RayDirInfo>) interaction,
    NBL_CONST_REF_ARG(AnisotropicMicrofacetCache<typename RayDirInfo::scalar_t>) _cache,
    NBL_CONST_REF_ARG(matrix<typename RayDirInfo::scalar_t, 2, 3>) ior,
    typename RayDirInfo::scalar_t ax,
    typename RayDirInfo::scalar_t ay)
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
        return beckmann_aniso_height_correlated_cos_eval_wo_clamps(_cache.NdotH2,TdotH2,BdotH2, _sample.NdotL2,TdotL2,BdotL2, interaction.NdotV,interaction.NdotV2,TdotV2,BdotV2, _cache.VdotH, ior,ax,ax*ax,ay,ay*ay);
    }
    else
    {
        return vector<scalar_t, 3>(0.0, 0.0, 0.0);
    }
}


template <class IncomingRayDirInfo, class Fresnel>
    NBL_REQUIRES(frensel::frensel<frensel_t> && ray_dir_info::basic<IncomingRayDirInfo>)
using IsotropicBeckmann = IsotropicCookTorrance<IncomingRayDirInfo, Fresnel, ndf::IsotropicBeckmann>;

template <class IncomingRayDirInfo, class Fresnel>
    NBL_REQUIRES(frensel::frensel<frensel_t> && ray_dir_info::basic<IncomingRayDirInfo>)
using AnisotropicBeckmann = AnisotropicCookTorrance<IncomingRayDirInfo, Fresnel, ndf::Beckmann>;

}
}
}
}
}

#endif