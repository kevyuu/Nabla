// Copyright (C) 2018-2024 - DevSH Graphics Programming Sp. z O.O.
// This file is part of the "Nabla Engine".
// For conditions of distribution and use, see copyright notice in nabla.h

#ifndef _NBL_BUILTIN_HLSL_BXDF_GEOM_SMITH_COMMON_INCLUDED_
#define _NBL_BUILTIN_HLSL_BXDF_GEOM_SMITH_COMMON_INCLUDED_


#include <nbl/builtin/hlsl/bxdf/ndf.hlsl>

namespace nbl
{
namespace hlsl
{
namespace bxdf
{
namespace geom_smith
{

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar G1(Scalar lambda)
{
    return 1.0 / (1.0 + lambda);
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar G2(Scalar lambda_V, Scalar lambda_L)
{
    return 1.0 / (1.0 + lambda_V + lambda_L);
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar VNDF_pdf_wo_clamps(Scalar ndf, Scalar lambda_V, Scalar maxNdotV, NBL_REF_ARG(Scalar) onePlusLambda_V)
{
    onePlusLambda_V = 1.0+lambda_V;

    return ndf::microfacet_to_light_measure_transform(ndf/onePlusLambda_V,maxNdotV);
}
template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar VNDF_pdf_wo_clamps(Scalar ndf, Scalar lambda_V, Scalar absNdotV, bool transmitted, Scalar VdotH, Scalar LdotH, Scalar VdotHLdotH, Scalar orientedEta, Scalar reflectance, NBL_REF_ARG(Scalar) onePlusLambda_V)
{
    onePlusLambda_V = 1.0+lambda_V;

    return ndf::microfacet_to_light_measure_transform((transmitted ? (1.0-reflectance):reflectance)*ndf/onePlusLambda_V,absNdotV,transmitted,VdotH,LdotH,VdotHLdotH,orientedEta);
}

// for when you know the NDF and the uncorrelated smith masking function
template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar VNDF_pdf_wo_clamps(Scalar ndf, Scalar G1_over_2NdotV)
{
    return ndf*0.5*G1_over_2NdotV;
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar FVNDF_pdf_wo_clamps(Scalar fresnel_ndf, Scalar G1_over_2NdotV, Scalar absNdotV, bool transmitted, Scalar VdotH, Scalar LdotH, Scalar VdotHLdotH, Scalar orientedEta)
{
    Scalar FNG = fresnel_ndf * G1_over_2NdotV;
    Scalar factor = 0.5;
    if (transmitted)
    {
        const Scalar VdotH_etaLdotH = (VdotH + orientedEta * LdotH);
        // VdotHLdotH is negative under transmission, so this factor is negative
        factor *= -2.0 * VdotHLdotH / (VdotH_etaLdotH * VdotH_etaLdotH);
    }
    return FNG * factor;
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar VNDF_pdf_wo_clamps(Scalar ndf, Scalar G1_over_2NdotV, Scalar absNdotV, bool transmitted, Scalar VdotH, Scalar LdotH, Scalar VdotHLdotH, Scalar orientedEta, Scalar reflectance)
{
    Scalar FN = (transmitted ? (1.0 - reflectance) : reflectance) * ndf;
    
    return FVNDF_pdf_wo_clamps(FN, G1_over_2NdotV, absNdotV, transmitted, VdotH, LdotH, VdotHLdotH, orientedEta);
}

	
}
}
}
}


#endif