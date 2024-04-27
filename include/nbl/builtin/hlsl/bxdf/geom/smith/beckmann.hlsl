// Copyright (C) 2018-2024 - DevSH Graphics Programming Sp. z O.O.
// This file is part of the "Nabla Engine".
// For conditions of distribution and use, see copyright notice in nabla.h

#ifndef _NBL_BUILTIN_HLSL_BXDF_GEOM_SMITH_BECKMANN_INCLUDED_
#define _NBL_BUILTIN_HLSL_BXDF_GEOM_SMITH_BECKMANN_INCLUDED_



namespace nbl
{
namespace hlsl
{
namespace bxdf
{
namespace geom_smith
{
namespace beckmann
{

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar C2(Scalar NdotX2, Scalar a2)
{
    return NdotX2 / (a2 * (1.0 - NdotX2));
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar C2(Scalar TdotX2, Scalar BdotX2, Scalar NdotX2, Scalar ax2, Scalar ay2)
{
    return NdotX2/(TdotX2*ax2+BdotX2*ay2);
}

//G1 = 1/(1+_Lambda)
template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar Lambda(Scalar c2)
{
    Scalar c = sqrt(c2);
    Scalar nom = 1.0 - 1.259*c + 0.396*c2;
    Scalar denom = 2.181*c2 + 3.535*c;
    return lerp(0.0, nom/denom, c<1.6);
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar Lambda(Scalar NdotX2, Scalar a2)
{
    return Lambda(C2(NdotX2, a2));
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar Lambda(Scalar TdotX2, Scalar BdotX2, Scalar NdotX2, Scalar ax2, Scalar ay2)
{
    return Lambda(C2(TdotX2, BdotX2, NdotX2, ax2, ay2));
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar correlated(Scalar NdotV2, Scalar NdotL2, Scalar a2)
{
    Scalar c2 = C2(NdotV2, a2);
    Scalar L_v = Lambda(c2);
    c2 = C2(NdotL2, a2);
    Scalar L_l = Lambda(c2);
    return 1.0 / (1.0 + L_v + L_l);
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar correlated(Scalar TdotV2, Scalar BdotV2, Scalar NdotV2, Scalar TdotL2, Scalar BdotL2, Scalar NdotL2, Scalar ax2, Scalar ay2)
{
    Scalar c2 = C2(TdotV2, BdotV2, NdotV2, ax2, ay2);
    Scalar L_v = Lambda(c2);
    c2 = C2(TdotL2, BdotL2, NdotL2, ax2, ay2);
    Scalar L_l = Lambda(c2);
    return 1.0 / (1.0 + L_v + L_l);
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar G2_over_G1(Scalar lambdaV_plus_one, Scalar NdotL2, Scalar a2)
{
    Scalar lambdaL = Lambda(NdotL2, a2);

    return lambdaV_plus_one / (lambdaV_plus_one+lambdaL);
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar G2_over_G1(Scalar lambdaV_plus_one, Scalar TdotL2, Scalar BdotL2, Scalar NdotL2, Scalar ax2, Scalar ay2)
{
    Scalar c2 = C2(TdotL2, BdotL2, NdotL2, ax2, ay2);
	Scalar lambdaL = Lambda(c2);

    return lambdaV_plus_one / (lambdaV_plus_one + lambdaL);
}


}
}
}
}
}


#endif