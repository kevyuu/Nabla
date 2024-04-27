// Copyright (C) 2018-2024 - DevSH Graphics Programming Sp. z O.O.
// This file is part of the "Nabla Engine".
// For conditions of distribution and use, see copyright notice in nabla.h

#ifndef _NBL_BUILTIN_HLSL_BXDF_GEOM_SMITH_GGX_INCLUDED_
#define _NBL_BUILTIN_HLSL_BXDF_GEOM_SMITH_GGX_INCLUDED_

#include <nbl/builtin/hlsl/concepts.hlsl>

namespace nbl
{
namespace hlsl
{
namespace bxdf
{
namespace geom_smith
{
namespace ggx
{


template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)            
Scalar devsh_part(Scalar NdotX2, Scalar a2, Scalar one_minus_a2)
{
    return sqrt(a2+one_minus_a2*NdotX2);
}
                    
template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)           
Scalar devsh_part(Scalar TdotX2, Scalar BdotX2, Scalar NdotX2, Scalar ax2, Scalar ay2)
{
    return sqrt(TdotX2*ax2+BdotX2*ay2+NdotX2);
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)           
Scalar G1_wo_numerator(Scalar NdotX, Scalar NdotX2, Scalar a2, Scalar one_minus_a2)
{
    return 1.0 / (NdotX + devsh_part(NdotX2,a2,one_minus_a2));
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar G1_wo_numerator(Scalar NdotX, Scalar TdotX2, Scalar BdotX2, Scalar NdotX2, Scalar ax2, Scalar ay2)
{
    return 1.0 / (NdotX + devsh_part(TdotX2, BdotX2, NdotX2, ax2, ay2));
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar G1_wo_numerator(Scalar NdotX, Scalar devsh_part)
{
    return 1.0 / (NdotX + devsh_part);
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar correlated_wo_numerator(Scalar NdotV, Scalar NdotV2, Scalar NdotL, Scalar NdotL2, Scalar a2, Scalar one_minus_a2)
{
    Scalar Vterm = NdotL*devsh_part(NdotV2,a2,one_minus_a2);
    Scalar Lterm = NdotV*devsh_part(NdotL2,a2,one_minus_a2);
    return 0.5 / (Vterm + Lterm);
}
                    
/* depr
template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar correlated(Scalar NdotV, Scalar NdotV2, Scalar NdotL, Scalar NdotL2, Scalar a2, Scalar one_minus_a2)
{
    return 4.0*NdotV*NdotL*correlated_wo_numerator(NdotV, NdotV2, NdotL, NdotL2, a2, one_minus_a2);
}
*/

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar correlated_wo_numerator(Scalar NdotV, Scalar NdotV2, Scalar NdotL, Scalar NdotL2, Scalar a2)
{
    return correlated_wo_numerator(NdotV,NdotV2,NdotL,NdotL2,a2,1.0-a2);
}
                    
/* depr
template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar correlated(Scalar NdotV, Scalar NdotV2, Scalar NdotL, Scalar NdotL2, Scalar a2)
{
    return 4.0*NdotV*NdotL*correlated_wo_numerator(NdotV, NdotV2, NdotL, NdotL2, a2);
}
*/

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar correlated_wo_numerator(Scalar NdotV, Scalar TdotV2, Scalar BdotV2, Scalar NdotV2, Scalar NdotL, Scalar TdotL2, Scalar BdotL2, Scalar NdotL2, Scalar ax2, Scalar ay2)
{
    Scalar Vterm = NdotL*devsh_part(TdotV2,BdotV2,NdotV2,ax2,ay2);
    Scalar Lterm = NdotV*devsh_part(TdotL2,BdotL2,NdotL2,ax2,ay2);
    return 0.5 / (Vterm + Lterm);
}
                    
/* depr
template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar correlated(Scalar NdotV, Scalar TdotV2, Scalar BdotV2, Scalar NdotV2, Scalar NdotL, Scalar TdotL2, Scalar BdotL2, Scalar NdotL2, Scalar ax2, Scalar ay2)
{
    return 4.0*NdotV*NdotL*correlated_wo_numerator(NdotV, TdotV2, BdotV2, NdotV2, NdotL, TdotL2, BdotL2, NdotL2, ax2, ay2);
}
*/
                    
template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar G2_over_G1(Scalar NdotL, Scalar NdotL2, Scalar NdotV, Scalar NdotV2, Scalar a2, Scalar one_minus_a2)
{
    Scalar devsh_v = devsh_part(NdotV2,a2,one_minus_a2);
	Scalar G2_over_G1 = NdotL*(devsh_v + NdotV); // alternative `Vterm+NdotL*NdotV /// NdotL*NdotV could come as a parameter
	G2_over_G1 /= NdotV*devsh_part(NdotL2,a2,one_minus_a2) + NdotL*devsh_v;

    return G2_over_G1;
}
                    
template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar G2_over_G1_devsh(Scalar NdotL, Scalar NdotL2, Scalar NdotV, Scalar devsh_v, Scalar a2, Scalar one_minus_a2)
{
	Scalar G2_over_G1 = NdotL*(devsh_v + NdotV); // alternative `Vterm+NdotL*NdotV /// NdotL*NdotV could come as a parameter
	G2_over_G1 /= NdotV*devsh_part(NdotL2,a2,one_minus_a2) + NdotL*devsh_v;

    return G2_over_G1;
}
                    
template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar G2_over_G1(Scalar NdotL, Scalar TdotL2, Scalar BdotL2, Scalar NdotL2, Scalar NdotV, Scalar TdotV2, Scalar BdotV2, Scalar NdotV2, Scalar ax2, Scalar ay2)
{
    Scalar devsh_v = devsh_part(TdotV2,BdotV2,NdotV2,ax2,ay2);
	Scalar G2_over_G1 = NdotL*(devsh_v + NdotV);
	G2_over_G1 /= NdotV*devsh_part(TdotL2,BdotL2,NdotL2,ax2,ay2) + NdotL*devsh_v;

    return G2_over_G1;
}
                    
template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar G2_over_G1_devsh(Scalar NdotL, Scalar TdotL2, Scalar BdotL2, Scalar NdotL2, Scalar NdotV, Scalar devsh_v, Scalar ax2, Scalar ay2)
{
	Scalar G2_over_G1 = NdotL*(devsh_v + NdotV);
	G2_over_G1 /= NdotV*devsh_part(TdotL2,BdotL2,NdotL2,ax2,ay2) + NdotL*devsh_v;

    return G2_over_G1;
}


}
}
}
}
}


#endif