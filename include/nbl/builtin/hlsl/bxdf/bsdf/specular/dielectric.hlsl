// Copyright (C) 2018-2024 - DevSH Graphics Programming Sp. z O.O.
// This file is part of the "Nabla Engine".
// For conditions of distribution and use, see copyright notice in nabla.h

#ifndef _NBL_BUILTIN_HLSL_BXDF_BSDF_SPECULAR_DIELECTRIC_INCLUDED_
#define _NBL_BUILTIN_HLSL_BXDF_BSDF_SPECULAR_DIELECTRIC_INCLUDED_

#include <nbl/builtin/hlsl/math/functions.hlsl>
#include <nbl/builtin/hlsl/bxdf/common.hlsl>
#include <nbl/builtin/hlsl/bxdf/transmission.hlsl>
#include <nbl/builtin/hlsl/bxdf/fresnel.hlsl>

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

// usually `luminosityContributionHint` would be the Rec.709 luma coefficients (the Y row of the RGB to CIE XYZ matrix)
// its basically a set of weights that determine 
// assert(1.0==luminosityContributionHint.r+luminosityContributionHint.g+luminosityContributionHint.b);
// `remainderMetadata` is a variable in which the generator function returns byproducts of sample generation that would otherwise have to be redundantly calculated in `remainder_and_pdf`
template<class IncomingRayDirInfo>
    NBL_REQUIRES(ray_dir_info::basic<IncomingRayDirInfo>)
LightSample<IncomingRayDirInfo> thin_smooth_dielectric_cos_generate_wo_clamps(
    NBL_CONST_REF_ARG(vector<typename IncomingRayDirInfo::scalar_t, 3>) V,
    NBL_CONST_REF_ARG(vector<typename IncomingRayDirInfo::scalar_t, 3>) T,
    NBL_CONST_REF_ARG(vector<typename IncomingRayDirInfo::scalar_t, 3>) B,
    NBL_CONST_REF_ARG(vector<typename IncomingRayDirInfo::scalar_t, 3>) N,
    typename IncomingRayDirInfo::scalar_t NdotV,
    typename IncomingRayDirInfo::scalar_t absNdotV,
    NBL_REF_ARG(vector<typename IncomingRayDirInfo::scalar_t, 3>) u,
    NBL_CONST_REF_ARG(vector<typename IncomingRayDirInfo::scalar_t, 3>) eta2,
    NBL_CONST_REF_ARG(vector<typename IncomingRayDirInfo::scalar_t, 3>) luminosityContributionHint,
    NBL_REF_ARG(vector<typename IncomingRayDirInfo::scalar_t, 3>) remainderMetadata)
{
    using scalar_t = typename IncomingRayDirInfo::scalar_t;
    // we will only ever intersect from the outside
    const vector<scalar_t, 3> reflectance = fresnel::thindielectric_infinite_scatter(fresnel::dielectric_common(eta2,absNdotV));

    // we are only allowed one choice for the entire ray, so make the probability a weighted sum
    const scalar_t reflectionProb = dot(reflectance, luminosityContributionHint);

    scalar_t rcpChoiceProb;
    const bool transmitted = math::partitionRandVariable(reflectionProb, u.z, rcpChoiceProb);
    remainderMetadata = (transmitted ? (vector<scalar_t, 3>(1.0,1.0,1.0)-reflectance):reflectance)*rcpChoiceProb;
    
    const vector<scalar_t, 3> L = (transmitted ? vector<scalar_t, 3>(0.0,0.0,0.0):(N*2.0*NdotV))-V;
    return LightSample<IncomingRayDirInfo>::create(IncomingRayDirInfo::create(L), dot(V,L), T, B, N);
}


template<class IncomingRayDirInfo>
    NBL_REQUIRES(ray_dir_info::basic<IncomingRayDirInfo>)
LightSample<IncomingRayDirInfo> thin_smooth_dielectric_cos_generate_wo_clamps(
    NBL_CONST_REF_ARG(vector<typename IncomingRayDirInfo::scalar_t, 3>) V,
    NBL_CONST_REF_ARG(vector<typename IncomingRayDirInfo::scalar_t, 3>) T,
    NBL_CONST_REF_ARG(vector<typename IncomingRayDirInfo::scalar_t, 3>) B,
    NBL_CONST_REF_ARG(vector<typename IncomingRayDirInfo::scalar_t, 3>) N,
    typename IncomingRayDirInfo::scalar_t NdotV,
    typename IncomingRayDirInfo::scalar_t absNdotV,
    NBL_REF_ARG(vector<typename IncomingRayDirInfo::scalar_t, 3>) u,
    NBL_CONST_REF_ARG(vector<typename IncomingRayDirInfo::scalar_t, 3>) eta2,
    NBL_CONST_REF_ARG(vector<typename IncomingRayDirInfo::scalar_t, 3>) luminosityContributionHint)
{
    vector<typename IncomingRayDirInfo::scalar_t, 3> dummy;
    return thin_smooth_dielectric_cos_generate_wo_clamps<IncomingRayDirInfo>(V,T,B,N,NdotV,absNdotV,u,eta2,luminosityContributionHint,dummy);
}

template<class IncomingRayDirInfo>
    NBL_REQUIRES(ray_dir_info::basic<IncomingRayDirInfo>)
LightSample<IncomingRayDirInfo> thin_smooth_dielectric_cos_generate(
    NBL_CONST_REF_ARG(surface_interactions::Anisotropic<IncomingRayDirInfo>) interaction,
    NBL_REF_ARG(vector<typename IncomingRayDirInfo::scalar_t, 3>) u,
    NBL_CONST_REF_ARG(vector<typename IncomingRayDirInfo::scalar_t, 3> eta2,
    NBL_CONST_REF_ARG(vector<typename IncomingRayDirInfo::scalar_t, 3> luminosityContributionHint)
{
    return thin_smooth_dielectric_cos_generate_wo_clamps(interaction.V.dir,interaction.T,interaction.B,interaction.N,interaction.NdotV,abs(interaction.NdotV),u,eta2,luminosityContributionHint);
}


template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
quotient_and_pdf<vector<Scalar, 3>, Scalar> thin_smooth_dielectric_cos_quotient_and_pdf_wo_clamps(NBL_CONST_REF_ARG(vector<Scalar, 3>) remainderMetadata)
{
    Scalar pdf = 1.0 / 0.0; // should be reciprocal probability of the fresnel choice divided by 0.0, but would still be an INF.
    return quotient_and_pdf<vector<Scalar, 3>, Scalar>::create(remainderMetadata, pdf);
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
quotient_and_pdf<vector<Scalar, 3>, Scalar> thin_smooth_dielectric_cos_quotient_and_pdf_wo_clamps(
    bool transmitted, Scalar absNdotV, NBL_CONST_REF_ARG(vector<Scalar, 3>) eta2, NBL_CONST_REF_ARG(vector<Scalar, 3>) luminosityContributionHint)
{
    const vector<Scalar, 3> reflectance = fresnel::thindielectric_infinite_scatter(fresnel::dielectric_common(eta2,absNdotV));
    const vector<Scalar, 3> sampleValue = transmitted ? (vector<Scalar, 3>(1.0,1.0,1.0)-reflectance):reflectance;

    const Scalar sampleProb = dot(sampleValue,luminosityContributionHint);

    return thin_smooth_dielectric_cos_quotient_and_pdf_wo_clamps(sampleValue / sampleProb);
}

// for information why we don't check the relation between `V` and `L` or `N` and `H`, see comments for `transmission_cos_quotient_and_pdf` in `irr/builtin/glsl/bxdf/common_samples.hlsl`
template<class IncomingRayDirInfo>
    NBL_REQUIRES(ray_dir_info::basic<IncomingRayDirInfo>)
quotient_and_pdf<typename vector<IncomingRayDirInfo::scalar_t, 3>, typename IncomingRayDirInfo::scalar_t> thin_smooth_dielectric_cos_quotient_and_pdf(
    NBL_CONST_REF_ARG(LightSample<IncomingRayDirInfo>) _sample,
    NBL_CONST_REF_ARG(surface_interactions::Isotropic<IncomingRayDirInfo>) interaction,
    NBL_CONST_REF_ARG(vector<typename IncomingRayDirInfo::scalar_t, 3>) eta2,
    NBL_CONST_REF_ARG(vector<typename IncomingRayDirInfo::scalar_t, 3>) luminosityContributionHint)
{
    const bool transmitted = math::isTransmissionPath(interaction.NdotV,_sample.NdotL);
    return thin_smooth_dielectric_cos_quotient_and_pdf_wo_clamps(transmitted,abs(interaction.NdotV),eta2,luminosityContributionHint);
}

template<class IncomingRayDirInfo>
    NBL_REQUIRES(ray_dir_info::basic<IncomingRayDirInfo>)
LightSample<IncomingRayDirInfo> smooth_dielectric_cos_generate_wo_clamps(
    NBL_CONST_REF_ARG(vector<typename IncomingRayDirInfo::scalar_t, 3>) V,
    NBL_CONST_REF_ARG(vector<typename IncomingRayDirInfo::scalar_t, 3>) T,
    NBL_CONST_REF_ARG(vector<typename IncomingRayDirInfo::scalar_t, 3>) B,
    NBL_CONST_REF_ARG(vector<typename IncomingRayDirInfo::scalar_t, 3>) N,
    bool backside,
    typename IncomingRayDirInfo::scalar_t NdotV,
    typename IncomingRayDirInfo::scalar_t absNdotV,
    typename IncomingRayDirInfo::scalar_t NdotV2,
    NBL_REF_ARG(vector<typename IncomingRayDirInfo::scalar_t, 3>) u,
    typename IncomingRayDirInfo::scalar_t rcpOrientedEta,
    typename IncomingRayDirInfo::scalar_t orientedEta2,
    typename IncomingRayDirInfo::scalar_t rcpOrientedEta2,
    NBL_REF_ARG(bool) transmitted)
{
    using scalar_t = typename IncomingRayDirInfo::scalar_t;
    
    const scalar_t reflectance = fresnel::dielectric_common(orientedEta2,absNdotV);

    scalar_t rcpChoiceProb;
    transmitted = math::partitionRandVariable(reflectance, u.z, rcpChoiceProb);

    const vector<scalar_t, 3> L = math::reflect_refract(transmitted, V, N, backside, NdotV, NdotV2, rcpOrientedEta, rcpOrientedEta2);
    return LightSample<IncomingRayDirInfo>::create(IncomingRayDirInfo::create(L), dot(V,L), T, B, N);
}

template<class IncomingRayDirInfo>
    NBL_REQUIRES(ray_dir_info::basic<IncomingRayDirInfo>)
LightSample<IncomingRayDirInfo> smooth_dielectric_cos_generate(
    NBL_CONST_REF_ARG(surface_interactions::Anisotropic<IncomingRayDirInfo>) interaction,
    NBL_REF_ARG(vector<typename IncomingRayDirInfo::scalar_t, 3>) u,
    typename IncomingRayDirInfo::scalar_t eta)
{
    using scalar_t = typename IncomingRayDirInfo::scalar_t;

    scalar_t orientedEta, rcpOrientedEta;
    const bool backside = getOrientedEtas(orientedEta, rcpOrientedEta, interaction.NdotV, eta);
    
    bool dummy;
    return smooth_dielectric_cos_generate_wo_clamps(
        interaction.V.dir,
        interaction.T,interaction.B,interaction.N,
        backside,
        interaction.NdotV,
        abs(interaction.NdotV),
        interaction.NdotV*interaction.NdotV,
        u,
        rcpOrientedEta, orientedEta*orientedEta, rcpOrientedEta*rcpOrientedEta,
        dummy
    );
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
quotient_and_pdf<Scalar, Scalar> smooth_dielectric_cos_quotient_and_pdf(bool transmitted, Scalar rcpOrientedEta2)
{
    const Scalar pdf = 1.0 / 0.0; // should be reciprocal probability of the fresnel choice divided by 0.0, but would still be an INF.
    return quotient_and_pdf_scalar::create(transmitted ? rcpOrientedEta2 : 1.0, pdf);
}

// for information why we don't check the relation between `V` and `L` or `N` and `H`, see comments for `transmission_cos_quotient_and_pdf` in `irr/builtin/glsl/bxdf/common_samples.hlsl`
template<class IncomingRayDirInfo>
    NBL_REQUIRES(ray_dir_info::basic<IncomingRayDirInfo>)
quotient_and_pdf<typename IncomingRayDirInfo::scalar_t, typename IncomingRayDirInfo::scalar_t> smooth_dielectric_cos_quotient_and_pdf(
    NBL_CONST_REF_ARG(LightSample<IncomingRayDirInfo>) _sample,
    NBL_CONST_REF_ARG(surface_interactions::Isotropic<IncomingRayDirInfo>) interaction,
    typename IncomingRayDirInfo::scalar_t eta)
{
    const bool transmitted = math::isTransmissionPath(interaction.NdotV,_sample.NdotL);
    
    typename IncomingRayDirInfo::scalar_t dummy, rcpOrientedEta;
    const bool backside = math::getOrientedEtas(dummy, rcpOrientedEta, interaction.NdotV, eta);

    return smooth_dielectric_cos_quotient_and_pdf(transmitted,rcpOrientedEta);
}

}
}
}
}
}

#endif