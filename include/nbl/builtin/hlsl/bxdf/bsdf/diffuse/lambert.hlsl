// Copyright (C) 2018-2024 - DevSH Graphics Programming Sp. z O.O.
// This file is part of the "Nabla Engine".
// For conditions of distribution and use, see copyright notice in nabla.h

#ifndef _NBL_BUILTIN_HLSL_BXDF_BSDF_DIFFUSE_LAMBERT_INCLUDED_
#define _NBL_BUILTIN_HLSL_BXDF_BSDF_DIFFUSE_LAMBERT_INCLUDED_

#include <nbl/builtin/hlsl/bxdf/common.hlsl>
#include <nbl/builtin/hlsl/sampling/cos_weighted.hlsl>

namespace nbl
{
namespace hlsl
{
namespace bxdf
{
namespace bsdf
{
namespace diffuse
{

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar lambertian()
{
    return math::ReciprocalPi<Scalar>::value * 0.5;
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar lambertian_cos_eval_rec_2pi_factored_out_wo_clamps(Scalar absNdotL)
{
    return absNdotL;
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar lambertian_cos_eval_rec_2pi_factored_out(Scalar NdotL)
{
    return lambertian_cos_eval_rec_2pi_factored_out_wo_clamps(abs(NdotL));
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar lambertian_cos_eval_wo_clamps(Scalar absNdotL)
{
    return lambertian_cos_eval_rec_2pi_factored_out_wo_clamps(absNdotL) * lambertian();
}

template<class IncomingRayDirInfo>
    NBL_REQUIRES(ray_dir_info::basic<IncomingRayDirInfo>)
typename IncomingRayDirInfo::scalar_t lambertian_cos_eval(in LightSample<IncomingRayDirInfo> _sample)
{
    return lambertian_cos_eval_rec_2pi_factored_out(_sample.NdotL) * lambertian();
}

template<class IncomingRayDirInfo>
    NBL_REQUIRES(ray_dir_info::basic<IncomingRayDirInfo>)
LightSample<IncomingRayDirInfo> lambertian_cos_generate_wo_clamps(
    NBL_CONST_REF_ARG(vector<typename IncomingRayDirInfo::scalar_t, 3>) tangentSpaceV,
    NBL_CONST_REF_ARG(matrix<typename IncomingRayDirInfo::scalar_t, 3, 3>) m,
    NBL_REF_ARG(vector<typename IncomingRayDirInfo::scalar_t, 3>) u)
{
    vector<typename IncomingRayDirInfo::scalar_t, 3> L = sampling::projected_sphere_generate(u);

    return LightSample<IncomingRayDirInfo>::createTangentSpace(tangentSpaceV, IncomingRayDirInfo::create(L), m);
}

template<class IncomingRayDirInfo>
    NBL_REQUIRES(ray_dir_info::basic<IncomingRayDirInfo>)
LightSample<IncomingRayDirInfo> lambertian_cos_generate(
    NBL_CONST_REF_ARG(surface_interactions::Anisotropic<IncomingRayDirInfo>) interaction,
    NBL_CONST_REF_ARG(vector<typename IncomingRayDirInfo::scalar_t, 3>) u)
{
    return lambertian_cos_generate_wo_clamps(interaction.getTangentSpaceV(), interaction.getTangentFrame(), u);
}

template<class IncomingRayDirInfo>
    NBL_REQUIRES(ray_dir_info::basic<IncomingRayDirInfo>)
quotient_and_pdf<typename IncomingRayDirInfo::scalar_t, typename IncomingRayDirInfo::scalar_t> lambertian_cos_quotient_and_pdf_wo_clamps(typename IncomingRayDirInfo::scalar_t absNdotL)
{
    using scalar_t = typename IncomingRayDirInfo::scalar_t;
    scalar_t pdf;
    scalar_t q = sampling::projected_sphere_quotient_and_pdf(pdf, absNdotL);
    return quotient_and_pdf<scalar_t, scalar_t>::create(q, pdf);
}

template<class IncomingRayDirInfo>
    NBL_REQUIRES(ray_dir_info::basic<IncomingRayDirInfo>)
quotient_and_pdf<typename IncomingRayDirInfo::scalar_t, typename IncomingRayDirInfo::scalar_t> lambertian_cos_quotient_and_pdf(NBL_CONST_REF_ARG(LightSample<IncomingRayDirInfo>) s)
{
    using scalar_t = typename IncomingRayDirInfo::scalar_t;
    scalar_t pdf;
    scalar_t q = lambertian_cos_quotient_and_pdf_wo_clamps(pdf, abs(s.NdotL));
    return quotient_and_pdf<scalar_t, scalar_t>::create(q, pdf);
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar lambertian_pdf_wo_clamps(Scalar absNdotL)
{
    return sampling::projected_sphere_pdf(absNdotL);
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
        return 0.5f * math::ReciprocalPi<scalar_t> * max(s.NdotL, 0.0f);
    }

    static sample_t generate(NBL_CONST_REF_ARG(surface_interactions::Anisotropic<IncomingRayDirInfo>) interaction, NBL_REF_ARG(vector<scalar_t, 3>) u)
    {
        vector<scalar_t, 3> L = sampling::projected_sphere_generate(u);

        const vector<scalar_t, 3> tangentSpaceV = interaction.getTangentSpaceV();
        const matrix<scalar_t, 3, 3> tangentFrame = interaction.getTangentFrame();

        return LightSample<IncomingRayDirInfo>::createTangentSpace(tangentSpaceV, IncomingRayDirInfo::create(L), tangentFrame);
    }

    q_pdf_t quotient_and_pdf(NBL_CONST_REF_ARG(sample_t) s, NBL_CONST_REF_ARG(base_t::interaction_t) interaction)
    {
        pdf_t pdf;
        pdf_t q = sampling::projected_sphere_quotient_and_pdf(pdf, abs(s.NdotL));

        return quotient_and_pdf<scalar_t, scalar_t>::create(q, pdf);
    }
};

}
}
}
}
}

#endif