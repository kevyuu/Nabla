// Copyright (C) 2018-2022 - DevSH Graphics Programming Sp. z O.O.
// This file is part of the "Nabla Engine".
// For conditions of distribution and use, see copyright notice in nabla.h

#ifndef _NBL_BUILTIN_HLSL_BXDF_NDF_INCLUDED_
#define _NBL_BUILTIN_HLSL_BXDF_NDF_INCLUDED_

#include <nbl/builtin/hlsl/bxdf/common.hlsl>
#include <nbl/builtin/hlsl/concepts.hlsl>

namespace nbl
{
namespace hlsl
{
namespace bxdf
{
namespace ndf
{

// general path
template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar microfacet_to_light_measure_transform(
    NBL_CONST_REF_ARG(Scalar) NDFcos, NBL_CONST_REF_ARG(Scalar) absNdotV, in bool transmitted, 
    NBL_CONST_REF_ARG(Scalar) VdotH, NBL_CONST_REF_ARG(Scalar) LdotH, NBL_CONST_REF_ARG(Scalar) VdotHLdotH, NBL_CONST_REF_ARG(Scalar) orientedEta)
{
    Scalar denominator = absNdotV;
    if (transmitted)
    {
        const Scalar VdotH_etaLdotH = (VdotH+orientedEta*LdotH);
        // VdotHLdotH is negative under transmission, so thats denominator is negative
        denominator *= -VdotH_etaLdotH*VdotH_etaLdotH;
    }
    return NDFcos*(transmitted ? VdotHLdotH:0.25)/denominator;
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar microfacet_to_light_measure_transform(NBL_CONST_REF_ARG(Scalar) NDFcos, NBL_CONST_REF_ARG(Scalar) maxNdotV)
{
    return 0.25*NDFcos/maxNdotV;
}


namespace ggx
{

// specialized factorizations for GGX
template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar microfacet_to_light_measure_transform(NBL_CONST_REF_ARG(Scalar) NDFcos_already_in_reflective_dL_measure, NBL_CONST_REF_ARG(Scalar) absNdotL, in bool transmitted, NBL_CONST_REF_ARG(Scalar) VdotH, NBL_CONST_REF_ARG(Scalar) LdotH, NBL_CONST_REF_ARG(Scalar) VdotHLdotH, NBL_CONST_REF_ARG(Scalar) orientedEta)
{
    Scalar factor = absNdotL;
    if (transmitted)
    {
        const Scalar VdotH_etaLdotH = (VdotH+orientedEta*LdotH);
        // VdotHLdotH is negative under transmission, so thats denominator is negative
        factor *= -4.0*VdotHLdotH/(VdotH_etaLdotH*VdotH_etaLdotH);
    }
    return NDFcos_already_in_reflective_dL_measure*factor;
}

Scalar microfacet_to_light_measure_transform(NBL_CONST_REF_ARG(Scalar) NDFcos_already_in_reflective_dL_measure, NBL_CONST_REF_ARG(Scalar) maxNdotL)
{
    return NDFcos_already_in_reflective_dL_measure*maxNdotL;
}

}


NBL_CONCEPT_TYPE_PARAMS(typename T)
NBL_CONCEPT_SIGNATURE(ndf, T t)
NBL_CONCEPT_BODY
(
    { T::scalr_t() } -> concepts::scalar;
    { T::ray_dir_info_t() } -> ray_dir_info::basic;

    // Note that generation is always anisotropic,
    // hence both interaction and microfacet cache must be anisotropic ones.
    { T::generateH(ray_dir_info_t(), vector<T::scalar_t, 3>, AnisotropicMicrofacetCache<T::scalar_t>()) } -> concepts::same_as<vector<T::scalar_t, 3>>

    // Isotropic NDF evaluators:
    { T::D(T::scalar_t()) } -> concepts::scalar;
    { T::Lambda(T::scalar_t()) } -> concepts::scalar;

    // Anisotropic NDF evaluators:
    { T::D(T::scalar_t(), T::scalar_t(), T::scalar_t()) } -> concepts::scalar;
    { T::Lambda(T::scalar_t(), T::scalar_t(), T::scalar_t()) } -> concepts::scalar;
)



// forward declaration so we can explicitly specialize, e.g. for GGX where optimized forms of the functions provided by the trait exist
template<class NDF>
    NBL_REQUIRES(ndf<NDF>)
struct ndf_traits;

namespace impl
{
    template<class NDF>
    struct ndf_traits
    {
        using scalar_t = typename NDF::scalar_t;

        scalar_t G1(NBL_CONST_REF_ARG(scalar_t) NdotX2)
        {
            return scalar_t(1) / (scalar_t(1) + ndf.Lambda(NdotX2));
        }
        scalar_t G1(NBL_CONST_REF_ARG(scalar_t) TdotX2, NBL_CONST_REF_ARG(scalar_t) BdotX2, NBL_CONST_REF_ARG(scalar_t) NdotX2)
        {
            return scalar_t(1) / (scalar_t(1) + ndf.Lambda(TdotX2, BdotX2, NdotX2));
        }

        scalar_t G2(NBL_CONST_REF_ARG(scalar_t) NdotV2, NBL_CONST_REF_ARG(scalar_t) NdotL2) 
        { 
            return scalar_t(1) / (scalar_t(1) + ndf.Lambda(NdotV2) + ndf.Lambda(NdotL2)); 
        }
        scalar_t G2(NBL_CONST_REF_ARG(scalar_t) TdotV2, NBL_CONST_REF_ARG(scalar_t) BdotV2, NBL_CONST_REF_ARG(scalar_t) NdotV2, NBL_CONST_REF_ARG(scalar_t) TdotL2, NBL_CONST_REF_ARG(scalar_t) BdotL2, NBL_CONST_REF_ARG(scalar_t) NdotL2) 
        { 
            return  scalar_t(1) / (scalar_t(1) + ndf.Lambda(TdotV2, BdotV2, NdotV2) + ndf.Lambda(TdotL2, BdotL2, NdotL2));
        }

        static scalar_t G2_over_G1_common(in scalar_t lambdaV, in scalar_t lambdaL)
        {
            const scalar_t lambdaV_plus_one = lambdaV + scalar_t(1);
            return lambdaV_plus_one / (lambdaL + lambdaV_plus_one);
        }
        scalar_t G2_over_G1(NBL_CONST_REF_ARG(scalar_t) NdotV2, NBL_CONST_REF_ARG(scalar_t) NdotL2) 
        {
            return G2_over_G1_common(ndf.Lambda(NdotV2), ndf.Lambda(NdotL2));
        }
        scalar_t G2_over_G1(
            NBL_CONST_REF_ARG(scalar_t) TdotV2, NBL_CONST_REF_ARG(scalar_t) BdotV2, NBL_CONST_REF_ARG(scalar_t) NdotV2,
            NBL_CONST_REF_ARG(scalar_t) TdotL2, NBL_CONST_REF_ARG(scalar_t) BdotL2, NBL_CONST_REF_ARG(scalar_t) NdotL2)
        {
            return G2_over_G1_common(
                ndf.Lambda(TdotV2, BdotV2, NdotV2), 
                ndf.Lambda(TdotL2, BdotL2, NdotL2)
            );
        }

        //
        // dHdL functions
        //
        // For BRDFs only:
        scalar_t dHdL(NBL_CONST_REF_ARG(scalar_t) NDFcos, NBL_CONST_REF_ARG(scalar_t) maxNdotV)
        {
            return microfacet_to_light_measure_transform(NDFcos, maxNdotV);
        }
        // For all BxDFs:
        scalar_t dHdL(
            NBL_CONST_REF_ARG(scalar_t) NDFcos, NBL_CONST_REF_ARG(scalar_t) absNdotV, in bool transmitted, NBL_CONST_REF_ARG(scalar_t) VdotH,
            NBL_CONST_REF_ARG(scalar_t) LdotH, NBL_CONST_REF_ARG(scalar_t) VdotHLdotH, NBL_CONST_REF_ARG(scalar_t) orientedEta)
        {
            return microfacet_to_light_measure_transform(NDFcos, absNdotV, transmitted, VdotH, LdotH, VdotHLdotH, orientedEta);
        }

        //
        // VNDF functions
        //
        // Statics:
        static scalar_t VNDF_static(in scalar_t d, in scalar_t lambda_V, NBL_CONST_REF_ARG(scalar_t) maxNdotV, out scalar_t onePlusLambda_V)
        {
            onePlusLambda_V = scalar_t(1) + lambda_V;

            return microfacet_to_light_measure_transform(d / onePlusLambda_V, maxNdotV);
        }
        static scalar_t VNDF_static(
            in scalar_t d, in scalar_t lambda_V, NBL_CONST_REF_ARG(scalar_t) absNdotV, in bool transmitted, NBL_CONST_REF_ARG(scalar_t) VdotH, 
            NBL_CONST_REF_ARG(scalar_t) LdotH, NBL_CONST_REF_ARG(scalar_t) VdotHLdotH, NBL_CONST_REF_ARG(scalar_t) orientedEta, 
            NBL_CONST_REF_ARG(scalar_t) reflectance, NBL_REF_ARG(scalar_t) onePlusLambda_V)
        {
            onePlusLambda_V = scalar_t(1) + lambda_V;

            return microfacet_to_light_measure_transform((transmitted ? (1.0 - reflectance) : reflectance) * d / onePlusLambda_V, absNdotV, transmitted, VdotH, LdotH, VdotHLdotH, orientedEta);
        }
        static scalar_t VNDF_static(in scalar_t d, in scalar_t G1_over_2NdotV)
        {
            return d * 0.5 * G1_over_2NdotV;
        }

        static scalar_t VNDF_fromLambda_impl(in scalar_t d, in scalar_t lambda, NBL_CONST_REF_ARG(scalar_t) maxNdotV)
        {
            scalar_t dummy;
            return VNDF_static(d, lambda, maxNdotV, dummy);
        }
        static scalar_t VNDF_fromG1_over_2NdotV_impl(in scalar_t d, in scalar_t G1_over_2NdotV)
        {
            return VNDF_static(d, G1_over_2NdotV);
        }

        

        // VNDF isotropic variants
        scalar_t VNDF(NBL_CONST_REF_ARG(scalar_t) NdotH2, NBL_CONST_REF_ARG(scalar_t) NdotV2, NBL_CONST_REF_ARG(scalar_t) maxNdotV)
        {
            const scalar_t d = ndf.D(NdotH2);
            const scalar_t lambda = ndf.Lambda(NdotV2);
            return VNDF_fromLambda_impl(d, lambda, maxNdotV);
        }
        scalar_t VNDF(NBL_CONST_REF_ARG(scalar_t) NdotH2, NBL_CONST_REF_ARG(scalar_t) G1_over_2NdotV)
        {
            const scalar_t d = ndf.D(NdotH2);
            return VNDF_fromG1_over_2NdotV_impl(d, G1_over_2NdotV);
        }
        scalar_t VNDF(
            NBL_CONST_REF_ARG(scalar_t) NdotH2, NBL_CONST_REF_ARG(scalar_t) NdotV2, 
            NBL_CONST_REF_ARG(scalar_t) absNdotV, in bool transmitted, NBL_CONST_REF_ARG(scalar_t) VdotH, NBL_CONST_REF_ARG(scalar_t) LdotH, NBL_CONST_REF_ARG(scalar_t) VdotHLdotH, NBL_CONST_REF_ARG(scalar_t) orientedEta, NBL_CONST_REF_ARG(scalar_t) reflectance, NBL_REF_ARG(scalar_t) onePlusLambda_V)
        {
            const scalar_t d = ndf.D(NdotH2);
            const scalar_t lambda = ndf.Lambda(NdotV2);

            return VNDF_static(d, lambda, absNdotV, transmitted, VdotH, LdotH, VdotHLdotH, orientedEta, reflectance, onePlusLambda_V);
        }

        // VNDF anisotropic variants
        scalar_t VNDF(
            NBL_CONST_REF_ARG(scalar_t) TdotH2, NBL_CONST_REF_ARG(scalar_t) BdotH2, NBL_CONST_REF_ARG(scalar_t) NdotH2,
            NBL_CONST_REF_ARG(scalar_t) TdotV2, NBL_CONST_REF_ARG(scalar_t) BdotV2, NBL_CONST_REF_ARG(scalar_t) NdotV2,
            NBL_CONST_REF_ARG(scalar_t) maxNdotV)
        {
            const scalar_t d = ndf.D(TdotH2, BdotH2, NdotH2);
            const scalar_t lambda = ndf.Lambda(TdotV2, BdotV2, NdotV2);
            return VNDF_fromLambda_impl(d, lambda, maxNdotV);
        }
        scalar_t VNDF(
            NBL_CONST_REF_ARG(scalar_t) TdotH2, NBL_CONST_REF_ARG(scalar_t) BdotH2, NBL_CONST_REF_ARG(scalar_t) NdotH2, 
            NBL_CONST_REF_ARG(scalar_t) G1_over_2NdotV)
        {
            const scalar_t d = ndf.D(TdotH2, BdotH2, NdotH2);
            return VNDF_fromG1_over_2NdotV_impl(d, G1_over_2NdotV);
        }
        scalar_t VNDF(
            NBL_CONST_REF_ARG(scalar_t) TdotH2, NBL_CONST_REF_ARG(scalar_t) BdotH2, NBL_CONST_REF_ARG(scalar_t) NdotH2,
            NBL_CONST_REF_ARG(scalar_t) TdotV2, NBL_CONST_REF_ARG(scalar_t) BdotV2, NBL_CONST_REF_ARG(scalar_t) NdotV2,
            NBL_CONST_REF_ARG(scalar_t) absNdotV, in bool transmitted, NBL_CONST_REF_ARG(scalar_t) VdotH, NBL_CONST_REF_ARG(scalar_t) LdotH, NBL_CONST_REF_ARG(scalar_t) VdotHLdotH, NBL_CONST_REF_ARG(scalar_t) orientedEta, NBL_CONST_REF_ARG(scalar_t) reflectance, NBL_REF_ARG(scalar_t) onePlusLambda_V)
        {
            const scalar_t d = ndf.D(TdotH2, BdotH2, NdotH2);
            const scalar_t lambda = ndf.Lambda(TdotV2, BdotV2, NdotV2);

            return VNDF_static(d, lambda, absNdotV, transmitted, VdotH, LdotH, VdotHLdotH, orientedEta, reflectance, onePlusLambda_V);
        }

        NDF ndf;
    };
}

// default specialization
template<class NDF>
struct ndf_traits : impl::ndf_traits<NDF> {};

}
}
}
}

#include <nbl/builtin/hlsl/bxdf/ndf/beckmann.hlsl>

#endif