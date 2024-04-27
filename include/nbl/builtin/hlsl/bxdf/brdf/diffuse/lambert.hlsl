// Copyright (C) 2018-2024 - DevSH Graphics Programming Sp. z O.O.
// This file is part of the "Nabla Engine".
// For conditions of distribution and use, see copyright notice in nabla.h
#ifndef _NBL_BUILTIN_HLSL_BXDF_BRDF_DIFFUSE_LAMBERT_INCLUDED_
#define _NBL_BUILTIN_HLSL_BXDF_BRDF_DIFFUSE_LAMBERT_INCLUDED_

#include <nbl/builtin/hlsl/bxdf/common.hlsl>
#include <nbl/builtin/hlsl/sampling/cos_weighted.hlsl>
#include <nbl/builtin/hlsl/concepts.hlsl>
#include <nbl/builtin/hlsl/math/constants.hlsl>

namespace nbl
{
namespace hlsl
{
namespace bxdf
{
namespace brdf
{
namespace diffuse
{

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar lambertian()
{
    return math::ReciprocalPi<Scalar>::value;
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar lambertian_cos_eval_rec_pi_factored_out_wo_clamps(Scalar maxNdotL)
{
   return maxNdotL;
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)                    
Scalar lambertian_cos_eval_rec_pi_factored_out(Scalar NdotL)
{
   return lambertian_cos_eval_rec_pi_factored_out_wo_clamps(max(NdotL,0.0f));
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar lambertian_cos_eval_wo_clamps(Scalar maxNdotL)
{
   return lambertian_cos_eval_rec_pi_factored_out_wo_clamps(maxNdotL)*lambertian();
}

template<class IncomingRayDirInfo>
    NBL_REQUIRES(ray_dir_info::basic<IncomingRayDirInfo>)
Scalar lambertian_cos_eval(NBL_CONST_REF_ARG(LightSample<IncomingRayDirInfo>) _sample)
{
    return lambertian_cos_eval_rec_pi_factored_out(_sample.NdotL)*lambertian();
}

template<class IncomingRayDirInfo>
    NBL_REQUIRES(ray_dir_info::basic<IncomingRayDirInfo>)
LightSample<IncomingRayDirInfo> lambertian_cos_generate_wo_clamps(
    NBL_CONST_REF_ARG(vector<typename IncomingRayDirInfo::scalar_t, 3>) tangentSpaceV,
    NBL_CONST_REF_ARG(matrix<typename IncomingRayDirInfo::scalar_t, 3, 3>) m,
    NBL_CONST_REF_ARG(vector<typename IncomingRayDirInfo::scalar_t, 2>) u)
{
    vector<typename IncomingRayDirInfo::scalar_t, 3> L = sampling::projected_hemisphere_generate(u);

    return LightSample<IncomingRayDirInfo>::createTangentSpace(tangentSpaceV,IncomingRayDirInfo::create(L),m);
}

template<class IncomingRayDirInfo>
    NBL_REQUIRES(ray_dir_info::basic<IncomingRayDirInfo>)
LightSample<IncomingRayDirInfo> lambertian_cos_generate(
                        NBL_CONST_REF_ARG(surface_interactions::Anisotropic<IncomingRayDirInfo>) interaction,
                        NBL_CONST_REF_ARG(vector<typename IncomingRayDirInfo::scalar_t, 2>) u)
{
    return lambertian_cos_generate_wo_clamps(interaction.getTangentSpaceV(),interaction.getTangentFrame(),u);
}


template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar lambertian_pdf_wo_clamps(Scalar maxNdotL)
{
    return sampling::projected_hemisphere_pdf(maxNdotL);
}

template<class IncomingRayDirInfo>
    NBL_REQUIRES(ray_dir_info::basic<IncomingRayDirInfo>)
Scalar lambertian_pdf(NBL_CONST_REF_ARG(LightSample<IncomingRayDirInfo>) s, NBL_CONST_REF_ARG(surface_interactions::Isotropic<IncomingRayDirInfo>) i)
{
    return lambertian_pdf_wo_clamps(max(s.NdotL,0.0f));
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
quotient_and_pdf<Scalar> lambertian_cos_quotient_and_pdf_wo_clamps(Scalar maxNdotL)
{
    Scalar pdf;
    Scalar q = sampling::projected_hemisphere_quotient_and_pdf(pdf, maxNdotL);

    return quotient_and_pdf_scalar::create(q, pdf);
}

template<class IncomingRayDirInfo>
    NBL_REQUIRES(ray_dir_info::basic<IncomingRayDirInfo>)
quotient_and_pdf<typename IncomingRayDirInfo::scalar_t> lambertian_cos_quotient_and_pdf(in LightSample<IncomingRayDirInfo> s)
{
    typename IncomingRayDirInfo::scalar_t pdf;
    typename IncomingRayDirInfo::scalar_t q = sampling::projected_hemisphere_quotient_and_pdf(pdf, max(s.NdotL,0.0f));

    return quotient_and_pdf_scalar::create(q, pdf);
}


template<class IncomingRayDirInfo>
    NBL_REQUIRES(ray_dir_info::basic<IncomingRayDirInfo>)
struct Lambertian
{
    using scalar_t = IncomingRayDirInfo::scalar_t;
    using spectrum_t = scalar_t;
    using pdf_t = scalar_t;
    using q_pdf_t = quotient_and_pdf<spectrum_t, pdf_t>;
    using sample_t = LightSample<IncomingRayDirInfo>;
    using interaction_t = surface_interactions::Isotropic<IncomingRayDirInfo>;

    static Lambertian create()
    {
        Lambertian lambertian;
        return lambertian;
    }

    spectrum_t eval(NBL_CONST_REF_ARG(sample_t) s, NBL_CONST_REF_ARG(interaction_t) interaction)
    {
        return math::ReciprocalPi<scalar_t>::value * max(s.NdotL, scalar_t(0.0));
    }

    static sample_t generate(NBL_CONST_REF_ARG(surface_interactions::Anisotropic<IncomingRayDirInfo>) interaction, NBL_REF_ARG(vector<scalar_t, 3>) u)
    {
        vector<scalar_t, 3> L = sampling::projected_hemisphere_generate(u.xy);

        const vector<scalar_t, 3> tangentSpaceV = interaction.getTangentSpaceV();
        const matrix<scalar_t, 3, 3> tangentFrame = interaction.getTangentFrame();

        return LightSample<IncomingRayDirInfo>::createTangentSpace(tangentSpaceV, IncomingRayDirInfo::create(L), tangentFrame);
    }

    q_pdf_t quotient_and_pdf(NBL_CONST_REF_ARG(sample_t) s, NBL_CONST_REF_ARG(interaction_t) interaction)
    {
        scalar_t pdf;
        scalar_t q = sampling::projected_hemisphere_quotient_and_pdf(pdf, max(s.NdotL, scalar_t(0.0));

        return quotient_and_pdf<scalar_t, scalar_t>::create(q, pdf);
    }
};

}
}
}
}
}

#endif