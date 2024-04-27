// Copyright (C) 2018-2024 - DevSH Graphics Programming Sp. z O.O.
// This file is part of the "Nabla Engine".
// For conditions of distribution and use, see copyright notice in nabla.h
#ifndef _NBL_BUILTIN_HLSL_SAMPLING_COS_WEIGHTED_INCLUDED_
#define _NBL_BUILTIN_HLSL_SAMPLING_COS_WEIGHTED_INCLUDED_

#include <nbl/builtin/hlsl/sampling/concentric_mapping.hlsl>

namespace nbl
{
namespace hlsl
{
namespace sampling
{

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
vector<Scalar, 3> projected_hemisphere_generate(NBL_CONST_REF_ARG(vector<Scalar, 2>) _sample)
{
    vector<Scalar, 2> p = concentricMapping(_sample*0.99999f+0.000005f);
    
    Scalar z = sqrt(max(0.0f, 1.0f - p.x*p.x - p.y*p.y));
    
    return vector<Scalar, 3>(p.x,p.y,z);
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar projected_hemisphere_pdf(Scalar L_z)
{
    return L_z * math::RECIPROCAL_PI;
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar projected_hemisphere_quotient_and_pdf(NBL_REF_ARG(Scalar) pdf, Scalar L_z)
{
	pdf = projected_hemisphere_pdf(L_z);
	return 1.0;
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar projected_hemisphere_quotient_and_pdf(NBL_REF_ARG(Scalar) pdf, NBL_CONST_REF_ARG(vector<Scalar, 3>) L)
{
	return projected_hemisphere_quotient_and_pdf(pdf,L.z);
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
vector<Scalar, 3> projected_sphere_generate(NBL_REF_ARG(vector<Scalar, 3>) _sample) // TODO, it should be `inout`, right?
{
    vector<Scalar, 3> retval = projected_hemisphere_generate(_sample.xy);
    const bool chooseLower = _sample.z>0.5f;
    retval.z = chooseLower ? (-retval.z):retval.z;
    if (chooseLower)
        _sample.z -= 0.5f;
    _sample.z *= 2.f;
    return retval;
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar projected_sphere_quotient_and_pdf(NBL_REF_ARG(Scalar) pdf, Scalar L_z)
{
    Scalar retval = projected_hemisphere_quotient_and_pdf(pdf,L_z);
    pdf *= 0.5f;
	return retval;
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar projected_sphere_quotient_and_pdf(NBL_REF_ARG(Scalar), NBL_CONST_REF_ARG(vector<Scalar, 3> L)
{
    return projected_sphere_quotient_and_pdf(pdf,L.z);
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar projected_sphere_pdf(ScalarL_z)
{
    return 0.5f*projected_hemisphere_pdf(L_z);
}
}
}
}

#endif