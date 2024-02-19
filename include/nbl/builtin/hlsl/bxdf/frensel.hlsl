// Copyright (C) 2018-2024 - DevSH Graphics Programming Sp. z O.O.
// This file is part of the "Nabla Engine".
// For conditions of distribution and use, see copyright notice in nabla.h

// TODO: rename Eta to something more fitting to let it be known that reciprocal Eta convention is used (ior_dest/ior_src)

#ifndef _NBL_BUILTIN_HLSL_BXDF_FRESNEL_INCLUDED_
#define _NBL_BUILTIN_HLSL_BXDF_FRESNEL_INCLUDED_

#include <nbl/builtin/hlsl/math/functions.hlsl>
#include <nbl/builtin/hlsl/bxdf/common.hlsl>
#include <nbl/builtin/hlsl/concepts.hlsl>

namespace nbl
{
namespace hlsl
{
namespace bxdf
{
namespace fresnel
{

// only works for implied IoR==1.333....
template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
vector<Scalar, 3> schlick(NBL_CONST_REF_ARG(vector<Scalar, 3>) F0, NBL_CONST_REF_ARG(Scalar) VdotH)
{
    Scalar x = 1.0 - VdotH;
    return F0 + (1.0 - F0) * x*x*x*x*x;
}

// TODO: provide a `nbl_glsl_fresnel_conductor_impl` that take `Eta` and `EtaLen2` directly
// conductors, only works for `CosTheta>=0`
template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
vector<Scalar, 3> conductor(NBL_CONST_REF_ARG(vector<Scalar, 3>) Eta, NBL_CONST_REF_ARG(vector<Scalar, 3>) Etak, NBL_CONST_REF_ARG(Scalar) CosTheta)
{  
   const Scalar CosTheta2 = CosTheta*CosTheta;
   const Scalar SinTheta2 = 1.0 - CosTheta2;

   const vector<Scalar, 3> EtaLen2 = Eta*Eta+Etak*Etak;
   const vector<Scalar, 3> etaCosTwice = Eta*CosTheta*2.0;

   const vector<Scalar, 3> rs_common = EtaLen2 + (CosTheta2).xxx;
   const vector<Scalar, 3> rs2 = (rs_common - etaCosTwice)/(rs_common + etaCosTwice);

   const vector<Scalar, 3> rp_common = EtaLen2*CosTheta2 + (1.0).xxx;
   const vector<Scalar, 3> rp2 = (rp_common - etaCosTwice)/(rp_common + etaCosTwice);
   
   return (rs2 + rp2)*0.5;
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
vector<Scalar, 3> conductor_impl(NBL_CONST_REF_ARG(vector<Scalar, 3>) Eta, NBL_CONST_REF_ARG(vector<Scalar, 3>) EtaLen2, NBL_CONST_REF_ARG(Scalar) CosTheta)
{  
   const Scalar CosTheta2 = CosTheta*CosTheta;
   const Scalar SinTheta2 = 1.0 - CosTheta2;

   const vector<Scalar, 3> etaCosTwice = Eta*CosTheta*2.0;

   const vector<Scalar, 3> rs_common = EtaLen2 + (CosTheta2).xxx;
   const vector<Scalar, 3> rs2 = (rs_common - etaCosTwice)/(rs_common + etaCosTwice);

   const vector<Scalar, 3> rp_common = EtaLen2*CosTheta2 + (1.0).xxx;
   const vector<Scalar, 3> rp2 = (rp_common - etaCosTwice)/(rp_common + etaCosTwice);
   
   return (rs2 + rp2)*0.5;
}


// dielectrics
template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
vector<Scalar, 3> dielectric_common(NBL_CONST_REF_ARG(Scalar) orientedEta2, NBL_CONST_REF_ARG(Scalar) AbsCosTheta)
{
    const Scalar SinTheta2 = 1.0-AbsCosTheta*AbsCosTheta;

    // the max() clamping can handle TIR when orientedEta2<1.0
    const Scalar t0 = sqrt(max(orientedEta2-SinTheta2,0.0));
    const Scalar rs = (AbsCosTheta - t0) / (AbsCosTheta + t0);

    const Scalar t2 = orientedEta2 * AbsCosTheta;
    const Scalar rp = (t0 - t2) / (t0 + t2);

    return (rs * rs + rp * rp) * 0.5;
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
vector<Scalar, 3> dielectric_common(NBL_CONST_REF_ARG(vector<Scalar, 3>) orientedEta2, NBL_CONST_REF_ARG(Scalar) AbsCosTheta)
{
   const Scalar SinTheta2 = 1.0-AbsCosTheta*AbsCosTheta;

   // the max() clamping can handle TIR when orientedEta2<1.0
   const vector<Scalar, 3> t0 = sqrt(max(vector<Scalar, 3>(orientedEta2)-SinTheta2, (0.0).xxx));
   const vector<Scalar, 3> rs = ((AbsCosTheta).xxx - t0) / ((AbsCosTheta).xxx + t0);

   const vector<Scalar, 3> t2 = orientedEta2*AbsCosTheta;
   const vector<Scalar, 3> rp = (t0 - t2) / (t0 + t2);

   return (rs*rs + rp*rp)*0.5;
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
vector<Scalar, 3> dielectric_frontface_only(NBL_CONST_REF_ARG(vector<Scalar, 3>) Eta, NBL_CONST_REF_ARG(Scalar) CosTheta)
{
    return dielectric_common(Eta*Eta,CosTheta);
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
vector<Scalar, 3> dielectric(NBL_CONST_REF_ARG(vector<Scalar, 3>) Eta, NBL_CONST_REF_ARG(Scalar) CosTheta)
{
    vector<Scalar, 3> orientedEta,rcpOrientedEta;
    math::getOrientedEtas(orientedEta,rcpOrientedEta,CosTheta,Eta);
    return dielectric_common(orientedEta*orientedEta,abs(CosTheta));
}


// gets the sum of all R, T R T, T R^3 T, T R^5 T, ... paths
template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
vector<Scalar, 3> thindielectric_infinite_scatter(NBL_CONST_REF_ARG(vector<Scalar, 3>) singleInterfaceReflectance)
{
    const vector<Scalar, 3> doubleInterfaceReflectance = singleInterfaceReflectance*singleInterfaceReflectance;
        
    return lerp(
        (singleInterfaceReflectance-doubleInterfaceReflectance)/(vector<Scalar, 3>(1.0,1.0,1.0)-doubleInterfaceReflectance)*2.0,
        vector<Scalar, 3>(1.0,1.0,1.0),
        (doubleInterfaceReflectance > vector<Scalar, 3>(0.9999,0.9999,0.9999))
    );
}

template <typename Scalar>
    NBL_REQUIRES(concepts::scalar<Scalar>)
Scalar thindielectric_infinite_scatter(NBL_CONST_REF_ARG(Scalar) singleInterfaceReflectance)
{
    const Scalar doubleInterfaceReflectance = singleInterfaceReflectance*singleInterfaceReflectance;
        
    return doubleInterfaceReflectance>0.9999 ? 1.0:((singleInterfaceReflectance-doubleInterfaceReflectance)/(1.0-doubleInterfaceReflectance)*2.0);
}


NBL_CONCEPT_TYPE_PARAMS(typename T)
NBL_CONCEPT_SIGNATURE(frensel, T t)
NBL_CONCEPT_BODY
(
    { T::scalar_t } -> concepts::scalar;
    { T::spectrum_t } -> spectral_of<T::scalar_t>;

    { T::operator()(T::scalar_t()) } -> concepts::same_as<T::spectrum_t>;
)


template <typename Scalar, typename Spectrum>
    NBL_REQUIRES(concepts::scalar<Scalar> && spectral_of<Spectrum, Scalar>)
struct FresnelSchlick
{
    using scalar_t = Scalar;
    using spectrum_t = Spectrum;

    spectrum_t F0;

    static FresnelSchlick create(NBL_CONST_REF_ARG(spectrum_t) _F0)
    {
        FresnelSchlick fs;
        fs.F0 = _F0;
        return fs;
    }

    spectrum_t operator()(NBL_CONST_REF_ARG(Scalar) cosTheta)
    {
        Scalar x = 1.0 - cosTheta;
        return F0 + (1.0 - F0) * x*x*x*x*x;
    }
};


using FresnelSchlickScalar  = FresnelSchlick<float, float>;
using FresnelSchlickRGB     = FresnelSchlick<float, float3>;


template <typename Scalar, typename Spectrum>
    NBL_REQUIRES(concepts::scalar<Scalar> && spectral_of<Spectrum, Scalar>)
struct FresnelConductor
{

    using scalar_t = Scalar;
    using spectrum_t = Spectrum;

    Spectrum eta;
    Spectrum etak;

    static FresnelConductor<Spectrum> create(NBL_CONST_REF_ARG(Spectrum) _eta, NBL_CONST_REF_ARG(Spectrum) _etak)
    {
        FresnelConductor<Spectrum> f;

        f.eta     = _eta;
        f.etak    = _etak;

        return f;
    }

    Spectrum operator()(NBL_CONST_REF_ARG(Scalar) cosTheta)
    {
       const Scalar CosTheta2 = cosTheta*cosTheta;
       const Scalar SinTheta2 = 1.0 - CosTheta2;

       const Spectrum EtaLen2 = eta*eta + etak*etak;
       const Spectrum etaCosTwice = eta*cosTheta*2.0;

       const Spectrum rs_common = EtaLen2 + (CosTheta2).xxx;
       const Spectrum rs2 = (rs_common - etaCosTwice)/(rs_common + etaCosTwice);

       const Spectrum rp_common = EtaLen2*CosTheta2 + Spectrum(1);
       const Spectrum rp2 = (rp_common - etaCosTwice)/(rp_common + etaCosTwice);
   
       return (rs2 + rp2)*0.5;
    }
};


using FresnelConductorScalar    = FresnelConductor<float, float>;
using FresnelConductorRGB       = FresnelConductor<float, float3>;


template <typename Scalar, typename Spectrum>
    NBL_REQUIRES(concepts::scalar<Scalar> && spectral_of<Spectrum, Scalar>)
struct FresnelDielectric
{
    using scalar_t = Scalar;
    using spectrum_t = Spectrum;

    Spectrum eta;

    static FresnelDielectric<Spectrum> create(in Spectrum _eta)
    {
        FresnelDielectric<Spectrum> f;

        f.eta = _eta;

        return f;
    }

    Spectrum operator()(NBL_CONST_REF_ARG(Scalar) cosTheta)
    {
        Spectrum orientedEta, rcpOrientedEta;
        math::getOrientedEtas<Spectrum>(orientedEta, rcpOrientedEta, cosTheta, eta);

        const Scalar AbsCosTheta = abs(cosTheta);
        const Spectrum orientedEta2 = orientedEta * orientedEta;

        const Scalar SinTheta2 = 1.0-AbsCosTheta*AbsCosTheta;

        // the max() clamping can handle TIR when orientedEta2<1.0
        const Spectrum t0 = sqrt(max(orientedEta2-SinTheta2, Spectrum(0)));
        const Spectrum rs = (AbsCosTheta - t0) / (AbsCosTheta + t0);

        const Spectrum t2 = orientedEta2*AbsCosTheta;
        const Spectrum rp = (t0 - t2) / (t0 + t2);

        return (rs*rs + rp*rp)*0.5;
    }
};

using FresnelDielectricScalar = FresnelDielectric<float, float>;
using FresnelDielectricRGB = FresnelDielectric<float, float3>;

}
}
}
}


#endif