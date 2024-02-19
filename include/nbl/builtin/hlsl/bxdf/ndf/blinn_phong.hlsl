// Copyright (C) 2018-2024 - DevSH Graphics Programming Sp. z O.O.
// This file is part of the "Nabla Engine".
// For conditions of distribution and use, see copyright notice in nabla.h

#ifndef _NBL_BUILTIN_HLSL_BXDF_NDF_BLINN_PHONG_INCLUDED_
#define _NBL_BUILTIN_HLSL_BXDF_NDF_BLINN_PHONG_INCLUDED_

#include <nbl/builtin/hlsl/math/constants.hlsl>
#include <nbl/builtin/hlsl/bxdf/ndf.hlsl>
// for beckmann sampling
#include <nbl/builtin/hlsl/bxdf/ndf/beckmann.hlsl>


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
Scalar blinn_phong(Scalar NdotH, Scalar n)
{
    return isinf(n) ? FLT_INF : math::ReciprocalPi<Scalar>::value*0.5*(n+2.0) * pow(NdotH,n);
}


//ashikhmin-shirley ndf
template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar blinn_phong(Scalar NdotH, Scalar one_minus_NdotH2_rcp, Scalar TdotH2, Scalar BdotH2, Scalar nx, Scalar ny)
{
    Scalar n = (TdotH2*ny + BdotH2*nx) * one_minus_NdotH2_rcp;

    return (isinf(nx)||isinf(ny)) ?  FLT_INF : sqrt((nx + 2.0)*(ny + 2.0))*math::ReciprocalPi<Scalar>::value*0.5 * pow(NdotH,n);
}

template <typename Scalar, typename RayDirInfo>
    NBL_REQUIRES(concepts::scalar<Scalar> && ray_dir_info::basic<RayDirInfo>)
struct IsotropicBlinnPhong
{
    using scalar_t = Scalar;
    using ray_dir_info_t = RayDirInfo;
	scalar_t n;

	//conversion between alpha and Phong exponent, Walter et.al.
	static scalar_t phong_exp_to_alpha2(scalar_t _n)
	{
		return 2.0 / (_n + 2.0);
	}
	//+INF for a2==0.0
	static scalar_t alpha2_to_phong_exp(scalar_t a2)
	{
		return 2.0 / a2 - 2.0;
	}

	static IsotropicBlinnPhong create(scalar_t _n)
	{
		IsotropicBlinnPhong bp;
		bp.n = _n;
		return bp;
	}

	static vector<scalar_t, 3> generateH_impl(NBL_CONST_REF_ARG(surface_interactions::Anisotropic<ray_dir_info_t>) interaction, NBL_REF_ARG(vector<scalar_t, 3>) u, NBL_REF_ARG(AnisotropicMicrofacetCache<scalar_t>) cache, scalar_t _ax, scalar_t _ay)
	{
		IsotropicBeckmann::generateH_impl(interaction, u, cache, _ax, _ay);
	}

	vector<scalar_t, 3> generateH(NBL_CONST_REF_ARG(surface_interactions::Anisotropic<ray_dir_info_t>) interaction, NBL_REF_ARG(vector<scalar_t, 3>) u, NBL_REF_ARG(AnisotropicMicrofacetCache<scalar_t>) cache)
	{
		scalar_t a = sqrt( phong_exp_to_alpha2(n) );
		return generateH_impl(interaction, u, cache, a, a);
	}

	scalar_t D(scalar_t NdotH2)
	{
		// here doing pow(NdotH2,n*0.5) instead of pow(NdotH,n) to keep compliant to common API (parameter of D() should be NdotH2)
		return isinf(n) ? FLT_INF : math::ReciprocalPi<scalar_t>::value * 0.5 * (n + 2.0) * pow(NdotH2, n*0.5);
	}

	scalar_t Lambda(scalar_t NdotX2)
	{
		// TODO
		// eh ill probably just use beckmann's lambda
		return 0.f / 0.f;//nan
	}
};


template <typename Scalar, typename RayDirInfo>
    NBL_REQUIRES(concepts::scalar<Scalar> && ray_dir_info::basic<RayDirInfo>)
struct BlinnPhong : IsotropicBlinnPhong<Scalar, RayDirInfo>
{
    using scalar_t = Scalar;
    using ray_dir_info_t = RayDirInfo;

	//scalar_t nx; // inherited from base as `n`
	scalar_t ny;

	static BlinnPhong create(scalar_t _nx, scalar_t _ny)
	{
		BlinnPhong bp;
		bp.n = _nx;
		bp.ny = _ny;
		return bp;
	}

	vector<scalar_t, 3> generateH(NBL_CONST_REF_ARG(surface_interactions::Anisotropic<ray_dir_info_t>) interaction, NBL_REF_ARG(vector<scalar_t, 3>) u, NBL_REF_ARG(AnisotropicMicrofacetCache<scalar_t>) cache)
	{
		scalar_t ax = sqrt(phong_exp_to_alpha2(n));
		scalar_t ay = sqrt(phong_exp_to_alpha2(ny));
		return generateH_impl(interaction, u, cache, ax, ay);
	}

	scalar_t D(scalar_t TdotH2, scalar_t BdotH2, scalar_t NdotH2)
	{
		scalar_t aniso_n = (TdotH2 * ny + BdotH2 * n) / (1.0 - NdotH2);

		return (isinf(n) || isinf(ny)) ? FLT_INF : sqrt((n + 2.0) * (ny + 2.0)) * math::ReciprocalPi<scalar_t>::value * 0.5 * pow(NdotH2, aniso_n*0.5);
	}

	scalar_t Lambda(scalar_t TdotH2, scalar_t BdotH2, scalar_t NdotX2)
	{
		// TODO
		// eh ill probably just use beckmann's lambda
		return 0.f / 0.f; //nan
	}
};

}
}
}
}

#endif