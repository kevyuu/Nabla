// Copyright (C) 2018-2024 - DevSH Graphics Programming Sp. z O.O.
// This file is part of the "Nabla Engine".
// For conditions of distribution and use, see copyright notice in nabla.h

#ifndef _NBL_BUILTIN_HLSL_BXDF_NDF_BECKMANN_INCLUDED_
#define _NBL_BUILTIN_HLSL_BXDF_NDF_BECKMANN_INCLUDED_

#include <nbl/builtin/hlsl/math/constants.hlsl>
#include <nbl/builtin/hlsl/math/functions.hlsl>
#include <nbl/builtin/hlsl/bxdf/common.hlsl>
#include <nbl/builtin/hlsl/bxdf/ndf.hlsl>


namespace nbl
{
namespace hlsl
{
namespace bxdf
{
namespace ndf
{

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar beckmann(NBL_CONST_REF_ARG(Scalar) a2, NBL_CONST_REF_ARG(Scalar) NdotH2)
{
    Scalar nom = exp( (NdotH2-1.0)/(a2*NdotH2) ); // exp(x) == exp2(x/log(2)) ?
    Scalar denom = a2*NdotH2*NdotH2;

    return math::ReciprocalPi<Scalar>::value * nom/denom;
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar beckmann(
    NBL_CONST_REF_ARG(Scalar) ax, NBL_CONST_REF_ARG(Scalar) ay, NBL_CONST_REF_ARG(Scalar) ax2, NBL_CONST_REF_ARG(Scalar) ay2, 
    NBL_CONST_REF_ARG(Scalar) TdotH2, NBL_CONST_REF_ARG(Scalar) BdotH2, NBL_CONST_REF_ARG(Scalar) NdotH2)
{
    Scalar nom = exp(-(TdotH2/ax2+BdotH2/ay2)/NdotH2);
    Scalar denom = ax * ay * NdotH2 * NdotH2;

    return math::ReciprocalPi<Scalar>::value * nom / denom;
}

template <typename Scalar, typename RayDirInfo>
    NBL_REQUIRES(concepts::scalar<Scalar> && ray_dir_info::basic<RayDirInfo>)
struct IsotropicBeckmann
{
    using scalar_t = Scalar;
    using ray_dir_info_t = RayDirInfo;

    Scalar a;
    Scalar a2;

    static IsotropicBeckmann create(Scalar _a)
    {
        IsotropicBeckmann b;
        b.a = _a;
        b.a2 = _a*_a;
        return b;
    }
    static IsotropicBeckmann create(Scalar _a, Scalar _a2)
    {
        IsotropicBeckmann b;
        b.a  = _a;
        b.a2 = _a2;
        return b;
    }


    static Scalar C2(Scalar NdotX2, Scalar _a2)
    {
        return NdotX2 / (_a2 * (1.0 - NdotX2));
    }
    static Scalar C2(Scalar TdotX2, Scalar BdotX2, Scalar NdotX2, Scalar _ax2, Scalar _ay2)
    {
        return NdotX2 / (TdotX2 * _ax2 + BdotX2 * _ay2);
    }

    template <class IncomingRayDirInfo>
    static vector<Scalar, 3> generateH_impl(NBL_CONST_REF_ARG(surface_interactions::Anisotropic<IncomingRayDirInfo>) interaction, NBL_REF_ARG(vector<Scalar, 3>) u, NBL_REF_ARG(AnisotropicMicrofacetCache<Scalar>) cache, Scalar _ax, Scalar _ay)
    {
        const vector<Scalar, 3> localV = interaction.getTangentSpaceV();

        //stretch
        vector<Scalar, 3> V = normalize(vector<Scalar, 3>(_ax * localV.x, _ay * localV.y, localV.z));

        vector<Scalar, 2> slope;
        if (V.z > 0.9999)//V.z=NdotV=cosTheta in tangent space
        {
            Scalar r = sqrt(-log(1.0 - u.x));
            Scalar sinPhi = sin(2.0 * math::Pi<Scalar>::value * u.y);
            Scalar cosPhi = cos(2.0 * math::Pi<Scalar>::value * u.y);
            slope = vector<Scalar, 2>(r, r) * vector<Scalar, 2>(cosPhi, sinPhi);
        }
        else
        {
            Scalar cosTheta = V.z;
            Scalar sinTheta = sqrt(1.0 - cosTheta * cosTheta);
            Scalar tanTheta = sinTheta / cosTheta;
            Scalar cotTheta = 1.0 / tanTheta;

            Scalar a = -1.0;
            Scalar c = math::erf(cosTheta);
            Scalar sample_x = max(u.x, 1.0e-6f);
            Scalar theta = acos(cosTheta);
            Scalar fit = 1.0 + theta * (-0.876 + theta * (0.4265 - 0.0594 * theta));
            Scalar b = c - (1.0 + c) * pow(1.0 - sample_x, fit);

            Scalar normalization = 1.0 / (1.0 + c + math::SQRT_RECIPROCAL_PI * tanTheta * exp(-cosTheta * cosTheta));

            const int ITER_THRESHOLD = 10;
            const Scalar MAX_ACCEPTABLE_ERR = 1.0e-5;
            int it = 0;
            Scalar value = 1000.0;
            while (++it<ITER_THRESHOLD && abs(value)>MAX_ACCEPTABLE_ERR)
            {
                if (!(b >= a && b <= c))
                    b = 0.5 * (a + c);

                Scalar invErf = math::erfInv(b);
                value = normalization * (1.0 + b + math::SQRT_RECIPROCAL_PI * tanTheta * exp(-invErf * invErf)) - sample_x;
                Scalar derivative = normalization * (1.0 - invErf * cosTheta);

                if (value > 0.0)
                    c = b;
                else
                    a = b;

                b -= value / derivative;
            }
            // TODO: investigate if we can replace these two erf^-1 calls with a box muller transform
            slope.x = math::erfInv(b);
            slope.y = math::erfInv(2.0f * max(u.y, 1.0e-6f) - 1.0f);
        }

        Scalar sinTheta = sqrt(1.0f - V.z * V.z);
        Scalar cosPhi = sinTheta == 0.0f ? 1.0f : clamp(V.x / sinTheta, -1.0f, 1.0f);
        Scalar sinPhi = sinTheta == 0.0f ? 0.0f : clamp(V.y / sinTheta, -1.0f, 1.0f);
        //rotate
        Scalar tmp = cosPhi * slope.x - sinPhi * slope.y;
        slope.y = sinPhi * slope.x + cosPhi * slope.y;
        slope.x = tmp;

        //unstretch
        slope = vector<Scalar, 2>(_ax, _ay) * slope;

        const vector<Scalar, 3> localH = normalize(vector<Scalar, 3>(-slope, 1.0));

        cache = AnisotropicMicrofacetCache::create(localV, localH);

        return localH;
    }

    vector<Scalar, 3> generateH(NBL_CONST_REF_ARG(surface_interactions::Anisotropic<IncomingRayDirInfo>) interaction, NBL_REF_ARG(vector<Scalar, 3>) u, NBL_REF_ARG(AnisotropicMicrofacetCache<Scalar>) cache)
    {
        return generateH_impl(interaction, u, cache, a, a);
    }

    scalar_t D(Scalar NdotH2)
    {
        Scalar nom = exp((NdotH2 - 1.0) / (a2 * NdotH2)); // exp(x) == exp2(x/log(2)) ?
        Scalar denom = a2 * NdotH2 * NdotH2;

        return math::ReciprocalPi<Scalar>::value * nom / denom;
    }

    scalar_t Lambda_impl(Scalar c2)
    {
        Scalar c = sqrt(c2);
        Scalar nom = 1.0 - 1.259 * c + 0.396 * c2;
        Scalar denom = 2.181 * c2 + 3.535 * c;
        return lerp(0.0, nom / denom, c < 1.6);
    }
    scalar_t Lambda(Scalar NdotX2)
    {
        return Lambda_impl(C2(NdotX2, a2));
    }

    // TODO what about aniso variants of D and Lambda?
    // return nan?
    // since we dont have SFINAE in HLSL, they must be defined for ndf_traits to compile with the type
};

template <typename Scalar, typename RayDirInfo>
    NBL_REQUIRES(concepts::scalar<Scalar> && ray_dir_info::basic<RayDirInfo>)
struct Beckmann : IsotropicBeckmann<Scalar, RayDirInfo>
{
    using scalar_t = Scalar;
    using ray_dir_info_t = RayDirInfo;

    static Beckmann create(Scalar _ax, Scalar _ay, Scalar _ax2, Scalar _ay2)
    {
        Beckmann b;
        b.a = _ax;
        b.ay = _ay;
        b.a2 = _ax2;
        b.ay2 = _ay2;
        return b;
    }
    static Beckmann create(Scalar _ax, Scalar _ay)
    {
        return create(_ax, _ay, _ax*_ax, _ay*_ay);
    }


    static Scalar C2(Scalar TdotX2, Scalar BdotX2, Scalar NdotX2, Scalar _ax2, Scalar _ay2)
    {
        return NdotX2 / (TdotX2 * _ax2 + BdotX2 * _ay2);
    }


    vector<Scalar, 3> generateH(NBL_CONST_REF_ARG(surface_interactions::Anisotropic<ray_dir_info_t>) interaction, NBL_REF_ARG(vector<Scalar, 3>) u, NBL_REF_ARG(AnisotropicMicrofacetCache<Scalar>) cache)
    {
        return generateH_impl(interaction, u, cache, a, ay);
    }


    Scalar D(Scalar TdotH2, Scalar BdotH2, Scalar NdotH2)
    {
        Scalar nom = exp(-(TdotH2 / a2 + BdotH2 / ay2) / NdotH2);
        Scalar denom = a * ay * NdotH2 * NdotH2;

        return math::ReciprocalPi<Scalar>::value * nom / denom;
    }

    Scalar Lambda(Scalar TdotX2, Scalar BdotX2, Scalar NdotX2)
    {
        return Lambda_impl(C2(TdotX2, BdotX2, NdotX2, a2, ay2));
    }

    //Scalar ax; // inherited from base as `a`
    //Scalar ax2; // inherited from base as `a2`
    Scalar ay;
    Scalar ay2;
};

}
}	
}
}

#endif