// Copyright (C) 2018-2024 - DevSH Graphics Programming Sp. z O.O.
// This file is part of the "Nabla Engine".
// For conditions of distribution and use, see copyright notice in nabla.h
#ifndef _NBL_BUILTIN_HLSL_BXDF_BRDF_SPECULAR_COMMON_INCLUDED_
#define _NBL_BUILTIN_HLSL_BXDF_BRDF_SPECULAR_COMMON_INCLUDED_

#include <nbl/builtin/hlsl/bxdf/common.hlsl>
#include <nbl/builtin/hlsl/bxdf/ndf.hlsl>

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

// Here anisotropic concepts would also satisfy the isotropic ones along with more requirements
// so this should work for both isotropic and anisotropic.
template <class fresnel_t, class ndf_t, class Sample, class Interaction, class MicrofacetCache>
    NBL_REQUIRES(frensel::frensel<frensel_t> && ndf::ndf<ndf_t> && light_sample<Sample> &&
        surface_interactions::isotropic<Interaction> && isotropic_microfacet_cache<MicrofacetCache>)
struct CookTorrance
{
    using cache_t = MicrofacetCache;
    using scalar_t = typename cache_t::scalar_t;
    using matrix_t = typename cache_t::matrix_t;
    using vector_t = typename cache_t::vector_t;
    using spectrum_t = vector_t;
    using pdf_t = scalar_t;
    using q_pdf_t = quotient_and_pdf<spectrum_t, pdf_t>;
    using sample_t = Sample;
    using interaction_t = Interaction;

    fresnel_t fresnel;
    ndf::ndf_traits<ndf_t> ndf;
};

template <class IncomingRayDirInfo, class fresnel_t, class ndf_t>
    NBL_REQUIRES(frensel::frensel<frensel_t> && ndf::ndf<ndf_t> && ray_dir_info::basic<IncomingRayDirInfo>)
struct IsotropicCookTorrance : CookTorrance<fresnel_t, ndf_t, LightSample<IncomingRayDirInfo>, surface_interactions::Isotropic<IncomingRayDirInfo>, IsotropicMicrofacetCache>
{
    using base_t = CookTorrance<fresnel_t, ndf_t, LightSample<IncomingRayDirInfo>, surface_interactions::Isotropic<IncomingRayDirInfo>, IsotropicMicrofacetCache>;

    typename base_t::spectrum_t eval(
        NBL_CONST_REF_ARG(typename base_t::sample_t) s,
        NBL_CONST_REF_ARG(typename base_t::interaction_t) interaction,
        NBL_CONST_REF_ARG(typename base_t::cache_t) cache)
    {
        if (interaction.NdotV > FLT_MIN)
        {
            typename base_t::scalar_t NG = base_t::ndf.ndf.D(cache.NdotH2);
            if (base_t::ndf.ndf.a2 > FLT_MIN)
                NG *= base_t::ndf.G2(interaction.NdotV2, s.NdotL2);

            const typename base_t::vector_t fr = base_t::fresnel(cache.VdotH);

            return fr * base_t::ndf.dHdL(NG, max(interaction.NdotV, 0.0f));
        }
        else
            return typename base_t::vector_t(0.0, 0.0, 0.0);
    }

    typename base_t::sample_t generate(
        NBL_CONST_REF_ARG(typename base_t::interaction_t) interaction,
        NBL_REF_ARG(typename base_t::vector_t) u,
        NBL_REF_ARG(typename base_t::cache_t) cache)
    {
        const typename base_t::vector_t localH = base_t::ndf.ndf.generateH(interaction, u, cache);

        const typename base_t::vector_t localV = interaction.getTangentSpaceV();
        const typename base_t::vector_t localL = math::reflect(localV, localH, cache.VdotH);

        return typename base_t::sample_t::createTangentSpace(localV, IncomingRayDirInfo::create(localL), interaction.getTangentFrame());
    }

    typename base_t::q_pdf_t quotient_and_pdf(
        NBL_CONST_REF_ARG(typename base_t::sample_t) s,
        NBL_CONST_REF_ARG(typename base_t::interaction_t) interaction,
        NBL_CONST_REF_ARG(typename base_t::cache_t) cache)
    {
        typename base_t::vector_t q = 0.0f;
        if (s.NdotL > FLT_MIN && interaction.NdotV > FLT_MIN)
        {
            const typename base_t::vector_t reflectance = base_t::fresnel(cache.VdotH);

            typename base_t::scalar_t G2_over_G1 = base_t::ndf.G2_over_G1(interaction.NdotV2, s.NdotL2);
            q = reflectance * G2_over_G1;
        }
        typename base_t::pdf_t pdf = base_t::ndf.VNDF(cache.NdotH2, interaction.NdotV2, max(interaction.NdotV, 0.0f));

        return typename base_t::q_pdf_t::create(q, pdf);
    }
};

template <class IncomingRayDirInfo, class fresnel_t, class ndf_t>
    NBL_REQUIRES(frensel::frensel<frensel_t> && ndf::ndf<ndf_t> && ray_dir_info::basic<IncomingRayDirInfo>)
struct AnisotropicCookTorrance : CookTorrance<fresnel_t, ndf_t, LightSample<IncomingRayDirInfo>, surface_interactions::Anisotropic<IncomingRayDirInfo>, AnisotropicMicrofacetCache>
{
    using base_t = CookTorrance<fresnel_t, ndf_t, LightSample<IncomingRayDirInfo>, surface_interactions::Anisotropic<IncomingRayDirInfo>, AnisotropicMicrofacetCache>;

    typename base_t::spectrum_t eval(
        NBL_CONST_REF_ARG(typename base_t::sample_t) s,
        NBL_CONST_REF_ARG(typename base_t::interaction_t) interaction,
        NBL_CONST_REF_ARG(typename base_t::cache_t) cache)
    {
        if (interaction.NdotV > FLT_MIN)
        {
            const typename base_t::scalar_t TdotL2 = s.TdotL * s.TdotL;
            const typename base_t::scalar_t BdotL2 = s.BdotL * s.BdotL;

            const typename base_t::scalar_t TdotV2 = interaction.TdotV * interaction.TdotV;
            const typename base_t::scalar_t BdotV2 = interaction.BdotV * interaction.BdotV;

            typename base_t::scalar_t NG = base_t::ndf.ndf.D(cache.TdotH2, cache.BdotH2, cache.NdotH2);
            if (base_t::ndf.ndf.ax > FLT_MIN || base_t::ndf.ndf.ay > FLT_MIN)
                NG *= base_t::ndf.G2(TdotV2, BdotV2, interaction.NdotV2, TdotL2, BdotL2, s.NdotL2);

            const typename base_t::vector_t fr = base_t::fresnel(cache.VdotH);

            return fr * base_t::ndf.dHdL(NG, max(interaction.NdotV, 0.0f));
        }
        else
            return base_t::vector_t(0.0, 0.0, 0.0);
    }

    typename base_t::sample_t generate(
        NBL_CONST_REF_ARG(typename base_t::interaction_t) interaction,
        NBL_REF_ARG(typename base_t::vector_t) u,
        NBL_REF_ARG(typename base_t::cache_t) cache)
    {
        const typename base_t::vector_t localH = base_t::ndf.ndf.generateH(interaction, u, cache);

        const typename base_t::vector_t localV = interaction.getTangentSpaceV();
        const typename base_t::vector_t localL = math::reflect(localV, localH, cache.VdotH);

        return base_t::sample_t::createTangentSpace(localV, IncomingRayDirInfo::create(localL), interaction.getTangentFrame());
    }

    typename base_t::q_pdf_t quotient_and_pdf(
        NBL_CONST_REF_ARG(typename base_t::sample_t) s,
        NBL_CONST_REF_ARG(typename base_t::interaction_t) interaction,
        NBL_CONST_REF_ARG(typename base_t::cache_t) cache)
    {
        const typename base_t::scalar_t TdotV2 = interaction.TdotV * interaction.TdotV;
        const typename base_t::scalar_t BdotV2 = interaction.BdotV * interaction.BdotV;

        typename base_t::vector_t q = 0.0f;
        if (s.NdotL > FLT_MIN && interaction.NdotV > FLT_MIN)
        {
            const typename base_t::scalar_t TdotL2 = s.TdotL * s.TdotL;
            const typename base_t::scalar_t BdotL2 = s.BdotL * s.BdotL;

            const typename base_t::vector_t reflectance = base_t::fresnel(cache.VdotH);

            typename base_t::scalar_t G2_over_G1 = base_t::ndf.G2_over_G1(TdotV2, BdotV2, interaction.NdotV2, TdotL2, BdotL2, s.NdotL2);
            q = reflectance * G2_over_G1;
        }

        const typename base_t::scalar_t TdotH2 = cache.TdotH * cache.TdotH;
        const typename base_t::scalar_t BdotH2 = cache.BdotH * cache.BdotH;

        typename base_t::pdf_t pdf = base_t::ndf.VNDF(
            TdotH2, BdotH2, cache.NdotH2,
            TdotV2, BdotV2, interaction.NdotV2,
            max(interaction.NdotV, 0.0f)
        );

        return typename base_t::q_pdf_t::create(q, pdf);
    }
};

}
}
}
}
}

#endif