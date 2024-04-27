// Copyright (C) 2018-2024 - DevSH Graphics Programming Sp. z O.O.
// This file is part of the "Nabla Engine".
// For conditions of distribution and use, see copyright notice in nabla.h
#ifndef _NBL_BUILTIN_HLSL_SAMPLING_CONCENTRIC_MAPPING_INCLUDED_
#define _NBL_BUILTIN_HLSL_SAMPLING_CONCENTRIC_MAPPING_INCLUDED_

#include <nbl/builtin/hlsl/math/functions.hlsl>
#include <nbl/builtin/hlsl/concepts.hlsl>
#include <nbl/builtin/hlsl/cpp_compat.hlsl>

namespace nbl
{
namespace hlsl
{
namespace sampling
{

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
vector<Scalar, 2> concentricMapping(NBL_CONST_REF_ARG(vector<Scalar, 2>) _u)
{
    //map [0;1]^2 to [-1;1]^2
    vector<Scalar, 2> u = 2.0f * _u - 1.0f;

    vector<Scalar, 2> p;
    if (all(u == vector<Scalar, 2>(0.0,0.0)))
        p = vector<Scalar, 2>(0.0,0.0);
    else
    {
        Scalar r;
        Scalar theta;
        if (abs(u.x) > abs(u.y)) {
            r = u.x;
            theta = 0.25f * math::PI * (u.y / u.x);
        }
        else {
            r = u.y;
            theta = 0.5f * math::PI - 0.25f * math::PI * (u.x / u.y);
        }
        // TODO: use nbl_glsl_sincos, but check that theta is in [-PI,PI]
        p = r * vector<Scalar, 2>(cos(theta), sin(theta));
    }

    return p;
}

}
}
}

#endif