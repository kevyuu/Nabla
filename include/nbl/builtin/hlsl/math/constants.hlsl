// Copyright (C) 2024 - DevSH Graphics Programming Sp. z O.O.
// This file is part of the "Nabla Engine".
// For conditions of distribution and use, see copyright notice in nabla.h

#ifndef _NBL_BUILTIN_HLSL_MATH_CONSTANTS_INCLUDED_
#define _NBL_BUILTIN_HLSL_MATH_CONSTANTS_INCLUDED_

#include <nbl/builtin/hlsl/type_traits.hlsl>
#include <nbl/builtin/hlsl/concepts.hlsl>

namespace nbl
{
namespace hlsl
{
namespace math
{

template <typename Scalar>
	NBL_REQUIRES(concepts::scalar<Scalar>)
struct ReciprocalPi : integral_constant<Scalar, 0.3183098861837906715377675267450287240689192914809128974953346881> {};

template <typename Scalar>
	NBL_REQUIRES(concepts::scalar<Scalar>)
struct Pi : integral_constant<Scalar, 3.1415926535897932384626433832795028841971693993751058209749445923> {};

}
}
}

#endif