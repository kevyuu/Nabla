// Copyright (C) 2018-2024 - DevSH Graphics Programming Sp. z O.O.
// This file is part of the "Nabla Engine".
// For conditions of distribution and use, see copyright notice in nabla.h

#ifndef _NBL_BUILTIN_HLSL_BXDF_BRDF_SPECULAR_GGX_INCLUDED_
#define _NBL_BUILTIN_HLSL_BXDF_BRDF_SPECULAR_GGX_INCLUDED_

#include <nbl/builtin/hlsl/bxdf/reflection.hlsl>
#include <nbl/builtin/hlsl/bxdf/fresnel.hlsl>
#include <nbl/builtin/hlsl/bxdf/ndf/ggx.hlsl>
#include <nbl/builtin/hlsl/bxdf/geom/smith/common.hlsl>
#include <nbl/builtin/hlsl/bxdf/geom/smith/ggx.hlsl>

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

//depr
/*
 float3 ggx_height_correlated_aniso_cos_eval(in BSDFAnisotropicParams params, in surface_interactions::Anisotropic<RayDirInfo> inter, in float2x3 ior, in float a2, in float2 atb, in float aniso)
{
    float g = geom_smith::ggx::height_correlated_aniso_wo_numerator(atb.x, atb.y, params.TdotL, interaction.TdotV, params.BdotL, interaction.BdotV, params.NdotL, interaction.NdotV);
    float ndf = ggx_burley_aniso(aniso, a2, params.TdotH, params.BdotH, params.NdotH);
    float3 fr = fresnel_conductor(ior[0], ior[1], params.VdotH);

    return params.NdotL * g*ndf*fr;
}
*/
//defined using NDF function with better API (compared to burley used above) and new impl of correlated smith
template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
 Scalar ggx_height_correlated_aniso_cos_eval_DG_wo_clamps(Scalar NdotH2, Scalar TdotH2, Scalar BdotH2, Scalar maxNdotL, Scalar NdotL2,
    Scalar TdotL2, Scalar BdotL2, Scalar maxNdotV, Scalar NdotV2, Scalar TdotV2, Scalar BdotV2, Scalar ax, Scalar ax2, Scalar ay, Scalar ay2)
{
    Scalar NG = ndf::ggx::aniso(TdotH2,BdotH2,NdotH2, ax, ay, ax2, ay2);
    if (ax>FLT_MIN || ay>FLT_MIN)
    {
        NG *= geom_smith::ggx::correlated_wo_numerator(
            maxNdotV, TdotV2, BdotV2, NdotV2,
            maxNdotL, TdotL2, BdotL2, NdotL2,
            ax2, ay2
        );
    }

    return NG;
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
 vector<Scalar, 3> ggx_height_correlated_aniso_cos_eval_wo_clamps(Scalar NdotH2, Scalar TdotH2, Scalar BdotH2, Scalar maxNdotL, Scalar NdotL2, Scalar TdotL2,
 Scalar BdotL2, Scalar maxNdotV, Scalar NdotV2, Scalar TdotV2, Scalar BdotV2, Scalar VdotH, in float2x3 ior, Scalar ax, Scalar ax2, Scalar ay, Scalar ay2)
{
    Scalar NG_already_in_reflective_dL_measure =
        ggx_height_correlated_aniso_cos_eval_DG_wo_clamps(NdotH2,TdotH2,BdotH2,maxNdotL,NdotL2,TdotL2,BdotL2,maxNdotV,NdotV2,TdotV2,BdotV2,ax,ax2,ay,ay2);

    vector<Scalar, 3> fr = fresnel::conductor(ior[0], ior[1], VdotH);
    return fr*ndf::ggx::microfacet_to_light_measure_transform(NG_already_in_reflective_dL_measure,maxNdotL);
}

template<class RayDirInfo>
    NBL_REQUIRES(ray_dir_info::basic<RayDirInfo>)
quotient_and_pdf<vector<typename RayDirInfo::scalar_t, 3>, typename RayDirInfo::scalar_t> ggx_height_correlated_aniso_cos_eval(
    NBL_CONST_REF_ARG(LightSample<RayDirInfo>) _sample,
    NBL_CONST_REF_ARG(surface_interactions::Anisotropic<RayDirInfo>) interaction,
    NBL_CONST_REF_ARG(AnisotropicMicrofacetCache<typename RayDirInfo::scalar_t>) _cache,
    NBL_CONST_REF_ARG(matrix<typename RayDirInfo::scalar_t, 2, 3>) ior,
    NBL_CONST_REF_ARG(typename RayDirInfo::scalar_t) ax,
    NBL_CONST_REF_ARG(typename RayDirInfo::scalar_t) ay)
{
    using scalar_t = typename RayDirInfo::scalar_t;
    if (_sample.NdotL>FLT_MIN && interaction.NdotV>FLT_MIN)
    {
        const scalar_t TdotH2 = _cache.TdotH*_cache.TdotH;
        const scalar_t BdotH2 = _cache.BdotH*_cache.BdotH;

        const scalar_t TdotL2 = _sample.TdotL*_sample.TdotL;
        const scalar_t BdotL2 = _sample.BdotL*_sample.BdotL;

        const scalar_t TdotV2 = interaction.TdotV*interaction.TdotV;
        const scalar_t BdotV2 = interaction.BdotV*interaction.BdotV;
        return ggx_height_correlated_aniso_cos_eval_wo_clamps(_cache.NdotH2, TdotH2, BdotH2, _sample.NdotL,_sample.NdotL2,TdotL2,BdotL2, interaction.NdotV,interaction.NdotV_squared,TdotV2,BdotV2, _cache.VdotH, ior, ax,ax*ax,ay,ay*ay);
    }
    else
        return vector<scalar_t, 3>(0.0,0.0,0.0);
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
 Scalar ggx_height_correlated_cos_eval_DG_wo_clamps(Scalar NdotH2, Scalar maxNdotL, Scalar NdotL2, Scalar maxNdotV, Scalar NdotV2, Scalar a2)
{
    Scalar NG = ndf::ggx::trowbridge_reitz(a2, NdotH2);
    if (a2>FLT_MIN)
        NG *= geom_smith::ggx::correlated_wo_numerator(maxNdotV, NdotV2, maxNdotL, NdotL2, a2);

    return NG;
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
vector<Scalar, 3> ggx_height_correlated_cos_eval_wo_clamps(Scalar NdotH2, Scalar maxNdotL, Scalar NdotL2, Scalar maxNdotV, Scalar NdotV2,
    Scalar VdotH, NBL_CONST_REF_ARG(matrix<Scalar, 2, 3>) ior, Scalar a2)
{
    Scalar NG_already_in_reflective_dL_measure = ggx_height_correlated_cos_eval_DG_wo_clamps(NdotH2, maxNdotL, NdotL2, maxNdotV, NdotV2, a2);

    vector<Scalar, 3> fr = fresnel::conductor(ior[0], ior[1], VdotH);

    return fr*ndf::ggx::microfacet_to_light_measure_transform(NG_already_in_reflective_dL_measure, maxNdotL);
}

template<class RayDirInfo>
    NBL_REQUIRES(ray_dir_info::basic<RayDirInfo>)
vector<typename RayDirInfo::scalar_t, 3> ggx_height_correlated_cos_eval(
    NBL_CONST_REF_ARG(LightSample<RayDirInfo>) _sample,
    NBL_CONST_REF_ARG(surface_interactions::Isotropic<RayDirInfo>) interaction,
    NBL_CONST_REF_ARG(IsotropicMicrofacetCache<typename RayDirInfo::scalar_t>) _cache,
    NBL_CONST_REF_ARG(matrix<typename RayDirInfo::scalar_t, 2, 3>) ior,
    typename RayDirInfo::scalar_t a2)
{
    if (_sample.NdotL>FLT_MIN && interaction.NdotV>FLT_MIN)
        return ggx_height_correlated_cos_eval_wo_clamps(_cache.NdotH2,max(_sample.NdotL,0.0),_sample.NdotL2, max(interaction.NdotV,0.0), interaction.NdotV_squared, _cache.VdotH,ior,a2);
    else
        return vector<typename RayDirInfo::scalar_t, 3>(0.0,0.0,0.0);
}



//Heitz's 2018 paper "Sampling the GGX Distribution of Visible Normals"
template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
vector<Scalar, 3> ggx_cos_generate(NBL_CONST_REF_ARG(vector<Scalar, 3>) localV, NBL_CONST_REF_ARG(vector<Scalar, 2>) u, Scalar _ax, Scalar _ay)
{
    vector<Scalar, 3> V = normalize(vector<Scalar, 3>(_ax*localV.x, _ay*localV.y, localV.z));//stretch view vector so that we're sampling as if roughness=1.0

    Scalar lensq = V.x*V.x + V.y*V.y;
    vector<Scalar, 3> T1 = lensq > 0.0 ? vector<Scalar, 3>(-V.y, V.x, 0.0)*rsqrt(lensq) : vector<Scalar, 3>(1.0,0.0,0.0);
    vector<Scalar, 3> T2 = cross(V,T1);

    Scalar r = sqrt(u.x);
    Scalar phi = 2.0 * math::PI * u.y;
    Scalar t1 = r * cos(phi);
    Scalar t2 = r * sin(phi);
    Scalar s = 0.5 * (1.0 + V.z);
    t2 = (1.0 - s)*sqrt(1.0 - t1*t1) + s*t2;
    
    //reprojection onto hemisphere
	//TODO try it wothout the& max(), not sure if -t1*t1-t2*t2>-1.0
    vector<Scalar, 3> H = t1*T1 + t2*T2 + sqrt(max(0.0, 1.0-t1*t1-t2*t2))*V;
    //unstretch
    return normalize(vector<Scalar, 3>(_ax*H.x, _ay*H.y, H.z));
}

// TODO: unifty the two following functions into `microfacet_BRDF_cos_generate_wo_clamps(float3 H,...)` and `microfacet_BRDF_cos_generate` or at least a auto declaration macro in lieu of a template
template<class RayDirInfo>
    NBL_REQUIRES(ray_dir_info::basic<RayDirInfo>)
LightSample<RayDirInfo> ggx_cos_generate_wo_clamps(
    NBL_CONST_REF_ARG(vector<typeneme RayDirInfo::scalar_t, 3>) localV,
    NBL_CONST_REF_ARG(matrix<typeneme RayDirInfo::scalar_t, 3, 3>) m,
    NBL_CONST_REF_ARG(vector<typeneme RayDirInfo::scalar_t, 2>) u,
    typeneme RayDirInfo::scalar_t _ax,
    typeneme RayDirInfo::scalar_t _ay,
    NBL_REF_ARG(AnisotropicMicrofacetCache<typeneme RayDirInfo::scalar_t>) _cache)
{
    const vector<typeneme RayDirInfo::scalar_t, 3> H = ggx_cos_generate(localV,u,_ax,_ay);
    
    _cache = AnisotropicMicrofacetCache::create(localV, H);
    vector<typeneme RayDirInfo::scalar_t, 3> localL;
    localL = math::reflect(localV, H, _cache.VdotH);
    
    return LightSample<RayDirInfo>::createTangentSpace(localV, RayDirInfo::create(localL), m);
}

template<class RayDirInfo>
    NBL_REQUIRES(ray_dir_info::basic<RayDirInfo>)
LightSample<RayDirInfo> ggx_cos_generate(
    NBL_CONST_REF_ARG(surface_interactions::Anisotropic<RayDirInfo>) interaction,
    NBL_CONST_REF_ARG(vector<typeneme RayDirInfo::scalar_t, 2>) u,
    typeneme RayDirInfo::scalar_t _ax,
    typeneme RayDirInfo::scalar_t _ay,
    NBL_REF_ARG(AnisotropicMicrofacetCache<typeneme RayDirInfo::scalar_t>) _cache)
{
    const vector<typeneme RayDirInfo::scalar_t, 3> localV = interaction.getTangentSpaceV();
    const matrix<typeneme RayDirInfo::scalar_t, 3, 3> m = interaction.getTangentFrame();
    return ggx_cos_generate_wo_clamps<RayDirInfo>(localV,m,u,_ax,_ay,_cache);
}


template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar ggx_pdf_wo_clamps(Scalar ndf, Scalar devsh_v, Scalar maxNdotV)
{
    return geom_smith::VNDF_pdf_wo_clamps(ndf, geom_smith::ggx::G1_wo_numerator(maxNdotV,devsh_v));
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar ggx_pdf_wo_clamps(Scalar NdotH2, Scalar maxNdotV, Scalar NdotV2, Scalar a2)
{
    const Scalar ndf = ndf::ggx::trowbridge_reitz(a2, NdotH2);
    const Scalar devsh_v = geom_smith::ggx::devsh_part(NdotV2, a2, 1.0-a2);

    return ggx_pdf_wo_clamps(ndf, devsh_v, maxNdotV);
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar ggx_pdf_wo_clamps(Scalar NdotH2, Scalar TdotH2, Scalar BdotH2, Scalar maxNdotV, Scalar NdotV2, Scalar TdotV2, Scalar BdotV2, Scalar ax, Scalar ay, Scalar ax2, Scalar ay2)
{
    const Scalar ndf = ndf::ggx::aniso(TdotH2,BdotH2,NdotH2, ax, ay, ax2, ay2);
    const Scalar devsh_v = geom_smith::ggx::devsh_part(TdotV2, BdotV2, NdotV2, ax2, ay2);

    return ggx_pdf_wo_clamps(ndf, devsh_v, maxNdotV);
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
quotient_and_pdf<vector<Scalar, 3>, Scalar> ggx_cos_quotient_and_pdf_wo_clamps(Scalar ndf, Scalar maxNdotL, Scalar NdotL2, 
    Scalar maxNdotV, Scalar NdotV2, NBL_CONST_REF_ARG(vector<Scalar, 3>) reflectance, Scalar a2)
{
    const Scalar one_minus_a2 = 1.0 - a2;
    const Scalar devsh_v = geom_smith::ggx::devsh_part(NdotV2, a2, one_minus_a2);
    const Scalar pdf = ggx_pdf_wo_clamps(ndf, devsh_v, maxNdotV);

    const Scalar G2_over_G1 = geom_smith::ggx::G2_over_G1_devsh(maxNdotL, NdotL2, maxNdotV, devsh_v, a2, one_minus_a2);

    return quotient_and_pdf<vector<Scalar, 3>, Scalar>::create(reflectance * G2_over_G1, pdf);
}

template<class RayDirInfo>
    NBL_REQUIRES(ray_dir_info::basic<RayDirInfo>)
quotient_and_pdf<vector<typename RayDirInfo::scalar_t, 3>, typename RayDirInfo::scalar_t> ggx_cos_quotient_and_pdf(
    NBL_CONST_REF_ARG(LightSample<RayDirInfo>) _sample,
    NBL_CONST_REF_ARG(surface_interactions::Isotropic<RayDirInfo>) interaction,
    NBL_CONST_REF_ARG(IsotropicMicrofacetCache<typename RayDirInfo::scalar_t>) _cache,
    NBL_CONST_REF_ARG(matrix<typename RayDirInfo::scalar_t, 2, 3>) ior,
    typename RayDirInfo::scalar_t a2)
{
    using scalar_t = typename RayDirInfo::scalar_t;
    const scalar_t one_minus_a2 = 1.0 - a2;
    const scalar_t ndf = ndf::ggx::trowbridge_reitz(a2, _cache.NdotH2);
    const scalar_t devsh_v = geom_smith::ggx::devsh_part(interaction.NdotV_squared, a2, one_minus_a2);
    const scalar_t pdf = ggx_pdf_wo_clamps(ndf, devsh_v, interaction.NdotV);

    quotient_and_pdf<vector<scalar_t, 3>, scalar_t> qpdf =
        quotient_and_pdf<vector<scalar_t, 3>, scalar_t>::create(vector<scalar_t, 3>(0.0, 0.0, 0.0), pdf);

    if (_sample.NdotL>FLT_MIN && interaction.NdotV>FLT_MIN)
    {
        const vector<scalar_t, 3> reflectance = fresnel::conductor(ior[0], ior[1], _cache.VdotH);
        const float G2_over_G1 = geom_smith::ggx::G2_over_G1_devsh(_sample.NdotL, _sample.NdotL2, interaction.NdotV, devsh_v, a2, one_minus_a2);

        qpdf.quotient = reflectance * G2_over_G1;
    }

    return qpdf;
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
quotient_and_pdf<vector<Scalar, 3>, Scalar> ggx_aniso_cos_quotient_and_pdf_wo_clamps(
    Scalar ndf, Scalar maxNdotL, Scalar NdotL2, Scalar TdotL2, Scalar BdotL2, Scalar maxNdotV,
    Scalar TdotV2, Scalar BdotV2, Scalar NdotV2, NBL_CONST_REF_ARG(vector<Scalar, 3>) reflectance, Scalar ax2,Scalar ay2)
{
    const Scalar devsh_v = geom_smith::ggx::devsh_part(TdotV2, BdotV2, NdotV2, ax2, ay2);
    const Scalar pdf = ggx_pdf_wo_clamps(ndf, devsh_v, maxNdotV);

    const Scalar G2_over_G1 = geom_smith::ggx::G2_over_G1_devsh(
        maxNdotL, TdotL2,BdotL2,NdotL2,
        maxNdotV, devsh_v,
        ax2, ay2
    );

    return quotient_and_pdf<vector<Scalar, 3>, Scalar> ::create(reflectance * G2_over_G1, pdf);
}

template<class RayDirInfo>
    NBL_REQUIRES(ray_dir_info::basic<RayDirInfo>)
quotient_and_pdf<vector<typename RayDirInfo::scalar_t, 3>, typename RayDirInfo::scalar_t>  ggx_aniso_cos_quotient_and_pdf(
    NBL_CONST_REF_ARG(LightSample<RayDirInfo>) _sample,
    NBL_CONST_REF_ARG(surface_interactions::Anisotropic<RayDirInfo>) interaction,
    NBL_CONST_REF_ARG(AnisotropicMicrofacetCache<typename RayDirInfo::scalar_t> _cache,
    NBL_CONST_REF_ARG(matrix<typename RayDirInfo::scalar_t, 2, 3>) ior,
    typename RayDirInfo::scalar_ ax,
    typename RayDirInfo::scalar_ ay)
{
    using scalar_t = typename RayDirInfo::scalar_t;
    const scalar_t ax2 = ax * ax;
    const scalar_t ay2 = ay * ay;

    const scalar_t TdotV2 = interaction.TdotV * interaction.TdotV;
    const scalar_t BdotV2 = interaction.BdotV * interaction.BdotV;
    const scalar_t NdotV2 = interaction.NdotV_squared;

    const scalar_t TdotH2 = _cache.TdotH * _cache.TdotH;
    const scalar_t BdotH2 = _cache.BdotH * _cache.BdotH;

    const scalar_t devsh_v = geom_smith::ggx::devsh_part(TdotV2, BdotV2, NdotV2, ax2, ay2);
    const scalar_t ndf = ndf::ggx::aniso(TdotH2, BdotH2, _cache.NdotH2, ax, ay, ax2, ay2);
    const scalar_t pdf = ggx_pdf_wo_clamps(ndf, devsh_v, interaction.NdotV);
    
    quotient_and_pdf<vector<scalar_t, 3>, scalar_t> qpdf =
        quotient_and_pdf<vector<scalar_t, 3>, scalar_t>::create(vector<scalar_t, 3>(0.0, 0.0, 0.0), pdf);
        
    if (_sample.NdotL>FLT_MIN && interaction.NdotV>FLT_MIN)
    {
        const scalar_t TdotL2 = _sample.TdotL*_sample.TdotL;
        const scalar_t BdotL2 = _sample.BdotL*_sample.BdotL;

        const vector<scalar_t, 3> reflectance = fresnel::conductor(ior[0], ior[1], _cache.VdotH);
        const scalar_t G2_over_G1 = geom_smith::ggx::G2_over_G1_devsh(
            _sample.NdotL, TdotL2, BdotL2, _sample.NdotL2,
            interaction.NdotV, devsh_v,
            ax2, ay2
        );

        qpdf.quotient = reflectance * G2_over_G1;
    }

    return qpdf;
}

}
}
}
}
}

#endif