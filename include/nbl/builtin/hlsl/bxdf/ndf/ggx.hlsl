// Copyright (C) 2018-2024 - DevSH Graphics Programming Sp. z O.O.
// This file is part of the "Nabla Engine".
// For conditions of distribution and use, see copyright notice in nabla.h

#ifndef _NBL_BUILTIN_HLSL_BXDF_NDF_GGX_INCLUDED_
#define _NBL_BUILTIN_HLSL_BXDF_NDF_GGX_INCLUDED_

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
namespace ggx
{


template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar trowbridge_reitz(NBL_CONST_REF_ARG(Scalar) a2, NBL_CONST_REF_ARG(Scalar) NdotH2)
{
    Scalar denom = NdotH2 * (a2 - 1.0) + 1.0;
    return a2* math::ReciprocalPi<Scalar>::value / (denom*denom);
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar burley_aniso(Scalar anisotropy, Scalar a2, Scalar TdotH, Scalar BdotH, Scalar NdotH)
{
	Scalar antiAniso = 1.0-anisotropy;
	Scalar atab = a2*antiAniso;
	Scalar anisoTdotH = antiAniso*TdotH;
	Scalar anisoNdotH = antiAniso*NdotH;
	Scalar w2 = antiAniso/(BdotH*BdotH+anisoTdotH*anisoTdotH+anisoNdotH*anisoNdotH*a2);
	return w2*w2*atab * math::ReciprocalPi<Scalar>::value;
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar aniso(NBL_CONST_REF_ARG(Scalar) TdotH2, NBL_CONST_REF_ARG(Scalar) BdotH2, NBL_CONST_REF_ARG(Scalar) NdotH2, NBL_CONST_REF_ARG(Scalar) ax, NBL_CONST_REF_ARG(Scalar) ay, NBL_CONST_REF_ARG(Scalar) ax2, NBL_CONST_REF_ARG(Scalar) ay2)
{
	Scalar a2 = ax*ay;
	Scalar denom = TdotH2/ax2 + BdotH2/ay2 + NdotH2;
	return math::ReciprocalPi<Scalar>::value / (a2 * denom * denom);
}


template <typename Scalar, typename IncommingRayDirInfo>
	NBL_REQUIRES(concepts::scalar<Scalar> && ray_dir_info::basic<IncommingRayDirInfo>)
struct IsotropicGGX
{
	using scalar_t = Scalar;
	using ray_dir_info_t = IncommingRayDirInfo;

	Scalar a;
	Scalar a2;

	static IsotropicGGX create(Scalar _a)
	{
		IsotropicGGX b;
		b.a = _a;
		b.a2 = _a * _a;
		return b;
	}
	static IsotropicGGX create(Scalar _a, Scalar _a2)
	{
		IsotropicGGX b;
		b.a = _a;
		b.a2 = _a2;
		return b;
	}

	static vector<Scalar,  3> generateH_impl(NBL_CONST_REF_ARG(surface_interactions::Anisotropic<IncommingRayDirInfo>) interaction, NBL_REF_ARG(vector<Scalar, 3>) u, out AnisotropicMicrofacetCache cache, NBL_CONST_REF_ARG(Scalar) _ax, NBL_CONST_REF_ARG(Scalar) _ay)
	{
		const vector<Scalar,  3> localV = interaction.getTangentSpaceV();

		vector<Scalar,  3> V = normalize(vector<Scalar, 3>(_ax * localV.x, _ay * localV.y, localV.z));//stretch view vector so that we're sampling as if roughness=1.0

		Scalar lensq = V.x * V.x + V.y * V.y;
		vector<Scalar,  3> T1 = lensq > 0.0 ? vector<Scalar, 3>(-V.y, V.x, 0.0) * rsqrt(lensq) : vector<Scalar, 3>(1.0, 0.0, 0.0);
		vector<Scalar,  3> T2 = cross(V, T1);

		Scalar r = sqrt(u.x);
		Scalar phi = 2.0 * math::PI * u.y;
		Scalar t1 = r * cos(phi);
		Scalar t2 = r * sin(phi);
		Scalar s = 0.5 * (1.0 + V.z);
		t2 = (1.0 - s) * sqrt(1.0 - t1 * t1) + s * t2;

		//reprojection onto hemisphere
		//TODO try it wothout the& max(), not sure if -t1*t1-t2*t2>-1.0
		vector<Scalar,  3> H = t1 * T1 + t2 * T2 + sqrt(max(0.0, 1.0 - t1 * t1 - t2 * t2)) * V;
		//unstretch
		const vector<Scalar,  3> localH = normalize(vector<Scalar, 3>(_ax * H.x, _ay * H.y, H.z));

		return localH;
	}

	vector<Scalar,  3> generateH(NBL_CONST_REF_ARG(surface_interactions::Anisotropic<IncommingRayDirInfo>) interaction, NBL_REF_ARG(vector<Scalar, 3>) u, NBL_REF_ARG(AnisotropicMicrofacetCache) cache)
	{
		return generateH_impl(interaction, u, cache, a, a);
	}

	Scalar D(NBL_CONST_REF_ARG(Scalar) NdotH2)
	{
		Scalar denom = NdotH2 * (a2 - 1.0) + 1.0;
		return a2 * math::ReciprocalPi<Scalar>::value / (denom * denom);
	}

	Scalar Lambda(NBL_CONST_REF_ARG(Scalar) NdotX2)
	{
		// TODO
		// ggx is a little special...
		return 0.f / 0.f;//nan
	}
};


template <typename Scalar, typename IncommingRayDirInfo>
	NBL_REQUIRES(concepts::scalar<Scalar> && ray_dir_info::basic<IncommingRayDirInfo>)
struct GGX : IsotropicGGX<Scalar, IncommingRayDirInfo>
{
	Scalar ay;
	Scalar ay2;

	static GGX create(Scalar _ax, Scalar _ay, Scalar _ax2, Scalar _ay2)
	{
		GGX b;
		b.a = _ax;
		b.ay = _ay;
		b.a2 = _ax2;
		b.ay2 = _ay2;
		return b;
	}
	static GGX create(Scalar _ax, Scalar _ay)
	{
		return create(_ax, _ay, _ax*_ax, _ay*_ay);
	}

	vector<Scalar,  3> generateH(NBL_CONST_REF_ARG(surface_interactions::Anisotropic<IncommingRayDirInfo>) interaction, NBL_REF_ARG(vector<Scalar, 3>) u, NBL_REF_ARG(AnisotropicMicrofacetCache<Scalar>) cache)
	{
		return generateH_impl(interaction, u, cache, a, ay);
	}

	Scalar D(NBL_CONST_REF_ARG(Scalar) TdotH2, NBL_CONST_REF_ARG(Scalar) BdotH2, NBL_CONST_REF_ARG(Scalar) NdotH2)
	{
		Scalar aa = a * ay;
		Scalar denom = TdotH2 / a2 + BdotH2 / ay2 + NdotH2;
		return math::ReciprocalPi<Scalar>::value / (aa * denom * denom);
	}

	Scalar Lambda(NBL_CONST_REF_ARG(Scalar) TdotH2, NBL_CONST_REF_ARG(Scalar) BdotH2, NBL_CONST_REF_ARG(Scalar) NdotX2)
	{
		// TODO
		// ggx is a little special...
		return 0.f / 0.f; //nan
	}
};

}
}
}
}
}



#endif