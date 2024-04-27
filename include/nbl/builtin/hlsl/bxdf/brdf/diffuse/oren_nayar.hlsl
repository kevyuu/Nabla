// Copyright (C) 2018-2024 - DevSH Graphics Programming Sp. z O.O.
// This file is part of the "Nabla Engine".
// For conditions of distribution and use, see copyright notice in nabla.h
#ifndef _NBL_BUILTIN_HLSL_BXDF_BRDF_DIFFUSE_OREN_NAYAR_INCLUDED_
#define _NBL_BUILTIN_HLSL_BXDF_BRDF_DIFFUSE_OREN_NAYAR_INCLUDED_

#include <nbl/builtin/hlsl/bxdf/brdf/diffuse/lambert.hlsl>

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
Scalar oren_nayar_cos_rec_pi_factored_out_wo_clamps(Scalar _a2, Scalar VdotL, Scalar maxNdotL, Scalar maxNdotV)
{
    // theta - polar angles
    // phi - azimuth angles
    Scalar a2 = _a2*0.5; //todo read about this&
    vector<Scalar, 2> AB = vector<Scalar, 2>(1.0, 0.0) + vector<Scalar, 2>(-0.5, 0.45) * vector<Scalar, 2>(a2, a2)/vector<Scalar, 2>(a2+0.33, a2+0.09);
    Scalar C = 1.0 / max(maxNdotL, maxNdotV);

    // should be equal to cos(phi)*sin(theta_i)*sin(theta_o)
    // where `phi` is the angle in the tangent plane to N, between L and V
    // and `theta_i` is the sine of the angle between L and N, similarily for `theta_o` but with V
    Scalar cos_phi_sin_theta = max(VdotL-maxNdotL*maxNdotV,0.0f);
    
    return (AB.x + AB.y * cos_phi_sin_theta * C);
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar oren_nayar_cos_eval_wo_clamps(Scalar a2, Scalar VdotL, Scalar maxNdotL, Scalar maxNdotV)
{
    return maxNdotL*math::RECIPROCAL_PI*oren_nayar_cos_rec_pi_factored_out_wo_clamps(a2,VdotL,maxNdotL,maxNdotV);
}

template<class IncomingRayDirInfo>
    NBL_REQUIRES(ray_dir_info::basic<IncomingRayDirInfo>)
typename IncomingRayDirInfo::scalar_t oren_nayar_cos_eval(
    NBL_CONST_REF_ARG(LightSample<IncomingRayDirInfo>) _sample,
    NBL_CONST_REF_ARG(surface_interactions::Isotropic<IncomingRayDirInfo>) inter,
    typename IncomingRayDirInfo::scalar_t a2)
{
    return oren_nayar_cos_eval_wo_clamps(a2, _sample.VdotL, max(_sample.NdotL,0.0f), max(inter.NdotV,0.0f));
}


template<class IncomingRayDirInfo>
    NBL_REQUIRES(ray_dir_info::basic<IncomingRayDirInfo>)
LightSample<IncomingRayDirInfo> oren_nayar_cos_generate_wo_clamps(
    NBL_CONST_REF_ARG(vector<typename IncomingRayDirInfo::scalar_t>, 3>) tangentSpaceV,
    NBL_CONST_REF_ARG(matrix<typename IncomingRayDirInfo::scalar_t>, 3, 3>) m,
    NBL_CONST_REF_ARG(vector<typename IncomingRayDirInfo::scalar_t>, 2>) u)
{
    // until we find something better
    return lambertian_cos_generate_wo_clamps<IncomingRayDirInfo>(tangentSpaceV, m, u);
}

template<class IncomingRayDirInfo>
    NBL_REQUIRES(ray_dir_info::basic<IncomingRayDirInfo>)
LightSample<IncomingRayDirInfo> oren_nayar_cos_generate(
    NBL_CONST_REF_ARG(surface_interactions::Anisotropic<IncomingRayDirInfo>) interaction,
    NBL_CONST_REF_ARG(vector<typename IncomingRayDirInfo::scalar_t>, 2>),
    typename IncomingRayDirInfo::scalar_t a2)
{
    return oren_nayar_cos_generate_wo_clamps<IncomingRayDirInfo>(getTangentSpaceV(interaction),getTangentFrame(interaction),u);
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar oren_nayar_pdf_wo_clamps(Scalar maxNdotL)
{
    return lambertian_pdf_wo_clamps(maxNdotL);
}

template<class IncomingRayDirInfo>
    NBL_REQUIRES(ray_dir_info::basic<IncomingRayDirInfo>)
typename IncomingRayDirInfo::scalar_t oren_nayar_pdf(
    NBL_CONST_REF_ARG(LightSample<IncomingRayDirInfo>) s,
    NBL_CONST_REF_ARG(surface_interactions::Isotropic<IncomingRayDirInfo>) i)
{
    return lambertian_pdf<IncomingRayDirInfo>(s, i);
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
quotient_and_pdf<Scalar, Scalar> oren_nayar_cos_quotient_and_pdf_wo_clamps(Scalar a2, Scalar VdotL, Scalar maxNdotL, Scalar maxNdotV)
{
    Scalar pdf = oren_nayar_pdf_wo_clamps(maxNdotL);
    return quotient_and_pdf<Scalar, Scalar>::create(
        oren_nayar_cos_rec_pi_factored_out_wo_clamps(a2,VdotL,maxNdotL,maxNdotV),
        pdf
    );
}

template<class IncomingRayDirInfo>
    NBL_REQUIRES(ray_dir_info::basic<IncomingRayDirInfo>)
quotient_and_pdf<Scalar, Scalar> oren_nayar_cos_quotient_and_pdf(
    NBL_CONST_REF_ARG(LightSample<IncomingRayDirInfo>) s,
    NBL_CONST_REF_ARG(surface_interactions::Isotropic<IncomingRayDirInfo>) interaction,
    typename IncomingRayDirInfo::scalar_t a2)
{
    return oren_nayar_cos_quotient_and_pdf_wo_clamps(a2,dot(interaction.V.getDirection(),s.L), max(s.NdotL,0.0f), max(interaction.NdotV,0.0f));
}


template<class IncomingRayDirInfo>
    NBL_REQUIRES(ray_dir_info::basic<IncomingRayDirInfo>)
struct OrenNayar
{
    using scalar_t = IncomingRayDirInfo::scalar_t;
    using spectrum_t = scalar_t;
    using pdf_t = scalar_t;
    using q_pdf_t = quotient_and_pdf<spectrum_t, pdf_t>;
    using sample_t = LightSample<IncomingRayDirInfo>;
    using interaction_t = surface_interactions::Isotropic<IncomingRayDirInfo>;

    scalar_t a2;
                            
    static OrenNayar<IncomingRayDirInfo> create(scalar_t _a2)
    {
        OrenNayar<IncomingRayDirInfo> orennayar;
        orennayar.a2 = _a2;
        return orennayar;
    }

    static scalar_t quotient(scalar_t _a2, scalar_t VdotL, scalar_t maxNdotL, scalar_t maxNdotV)
    {
        // theta - polar angles
        // phi - azimuth angles
        scalar_t a2 = _a2 * 0.5; //todo read about this
        vector<scalar_t, 2> AB = vector<scalar_t, 2>(1.0, 0.0) + vector<scalar_t, 2>(-0.5, 0.45) * vector<scalar_t, 2>(a2, a2) / vector<scalar_t, 2>(a2 + 0.33, a2 + 0.09);
        scalar_t C = 1.0 / max(maxNdotL, maxNdotV);

        // should be equal to cos(phi)*sin(theta_i)*sin(theta_o)
        // where `phi` is the angle in the tangent plane to N, between L and V
        // and `theta_i` is the sine of the angle between L and N, similarily for `theta_o` but with V
        scalar_t cos_phi_sin_theta = max(VdotL - maxNdotL * maxNdotV, 0.0f);

        return (AB.x + AB.y * cos_phi_sin_theta * C);
    }
                    
    scalar_t quotient(NBL_CONST_REF_ARG(sample_t) s, NBL_CONST_REF_ARG(interaction_t) interaction, NBL_REF_ARG(scalar_t) maxNdotL)
    {
        /*out*/maxNdotL = max(s.NdotL, 0.0f);
        const scalar_t maxNdotV = max(interaction.NdotV, 0.0f);

        return quotient(a2, s.VdotL, maxNdotL, maxNdotV);
    }

    spectrum_t eval(NBL_CONST_REF_ARG(sample_t) s, NBL_CONST_REF_ARG(interaction_t) interaction)
    {
        scalar_t maxNdotL;
        scalar_t q = quotient(s, interaction, maxNdotL);
        return math::RECIPROCAL_PI * maxNdotL * q;
    }

    typename base_t::sample_t generate(
        NBL_CONST_REF_ARG(surface_interactions::Anisotropic<IncomingRayDirInfo>) interaction,
        NBL_REF_ARG(vector<scalar_t, 3>) u)
    {
        return Lambertian<IncomingRayDirInfo>::generate(interaction, u);
    }

    q_pdf_t quotient_and_pdf(NBL_CONST_REF_ARG(sample_t) s, NBL_CONST_REF_ARG(interaction_t) interaction)
    {
        scalar_t maxNdotL;
        scalar_t q = quotient(s, interaction, maxNdotL);

        return quotient_and_pdf_scalar::create(q, maxNdotL * math::ReciprocalPi<scalar_t>::value);
    }
};

}
}
}
}
}

#endif