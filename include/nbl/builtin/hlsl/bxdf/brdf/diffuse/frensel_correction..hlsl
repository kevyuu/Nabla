// Copyright (C) 2018-2024 - DevSH Graphics Programming Sp. z O.O.
// This file is part of the "Nabla Engine".
// For conditions of distribution and use, see copyright notice in nabla.h
#ifndef _NBL_BUILTIN_HLSL_BXDF_BRDF_DIFFUSE_FRESNEL_CORRECTION_INCLUDED_
#define _NBL_BUILTIN_HLSL_BXDF_BRDF_DIFFUSE_FRESNEL_CORRECTION_INCLUDED_

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
vector<Scalar, 3> diffuseFresnelCorrectionFactor(NBL_CONST_REF_ARG(vector<Scalar, 3>) n, NBL_CONST_REF_ARG(vector<Scalar, 3>) n2)
{
    const Scalar C1 = 554.33;
    const Scalar C2 = 380.7;
    const Scalar C3 = 298.25;
    const Scalar C4 = 261.38;
    const Scalar C5 = 138.43;
    const Scalar C6 = 0.8078843897748912;
    const Scalar C7 = -1.67;
    const Scalar C8 = 0.1921156102251088;


    bool3 TIR = (n < 1.0);
    vec3 invdenum = lerp(vector<Scalar, 3>(1.0,1.0,1.0), vector<Scalar, 3>(1.0,1.0,1.0)/(n2*n2*(vec3(C1,C1,C1) - C2*n)), TIR);
    vec3 num = n*lerp(vector<Scalar, 3>(C8,C8,C8),n*C3 - C4*n2 + C5,TIR);
    num += lerp(vector<Scalar, 3>(C6,C6,C6),vector<Scalar, 3>(C7,C7,C7),TIR);
    return num*invdenum;
}

}
}
}
}
}

#endif