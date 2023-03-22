// Copyright (C) 2018-2022 - DevSH Graphics Programming Sp. z O.O.
// This file is part of the "Nabla Engine".
// For conditions of distribution and use, see copyright notice in nabla.h
#ifndef _NBL_BUILTIN_HLSL_BXDF_COMMON_INCLUDED_
#define _NBL_BUILTIN_HLSL_BXDF_COMMON_INCLUDED_

#include <nbl/builtin/hlsl/limits/numeric.hlsl>
//#include <nbl/builtin/hlsl/numbers.hlsl> // ???
#include <nbl/builtin/hlsl/math/functions.hlsl>

namespace nbl
{
namespace hlsl
{
namespace bxdf
{


namespace ray_dir_info
{

// no ray-differentials, nothing
struct Basic
{
  static Basic create(float3 dir)
  {
    Basic retval;
    retval.direction = dir;
    return retval;
  }

  // `tform` must be orthonormal
  void transform(float3x3 tform)
  {
    direction = mul(tform, direction); // TODO (make sure) is this correct multiplication order in HLSL? Original GLSL code is `tform * direction`
  }

  // `tform` must be orthonormal
  static Basic transform(Basic r, float3x3 tform)
  {
    r.transform(tform);
    return r;
  }

  float3 getDirection() {return direction;}

  Basic transmit()
  {
    Basic retval;
    retval.direction = -direction;
    return retval;
  }
  
  Basic reflect(const float3 N, const float directionDotN)
  {
    Basic retval;
    retval.direction = math::reflect(direction,N,directionDotN);
    return retval;
  }

  float3 direction;
};
// more to come!

}


namespace surface_interactions
{

template<class RayDirInfo>
struct Isotropic
{
  // WARNING: Changed since GLSL, now arguments need to be normalized!
  static Isotropic<RayDirInfo> create(const RayDirInfo normalizedV, const float3 normalizedN)
  {
    Isotropic<RayDirInfo> retval;
    retval.V = normalizedV;
    retval.N = normalizedN;

    retval.NdotV = dot(retval.N,retval.V.getDirection());
    retval.NdotV_squared = retval.NdotV*retval.NdotV;

    return retval;
  }

  RayDirInfo V;
  float3 N;
  float NdotV;
  float NdotV2; // old NdotV_squared
};

template<class RayDirInfo>
struct Anisotropic : Isotropic<RayDirInfo>
{
  // WARNING: Changed since GLSL, now arguments need to be normalized!
  static Anisotropic<RayDirInfo> create(
    const Isotropic<RayDirInfo> isotropic,
    const float3 normalizedT,
    const float normalizedB
  )
  {
    Anisotropic<RayDirInfo> retval;
    (Isotropic<RayDirInfo>) retval = isotropic;
    retval.T = normalizedT;
    retval.B = normalizedB;
    
    const float3 V = retval.getDirection();
    retval.TdotV = dot(V,retval.T);
    retval.BdotV = dot(V,retval.B);

    return retval;
  }
  static Anisotropic<RayDirInfo> create(const Isotropic<RayDirInfo> isotropic, const float3 normalizedT)
  {
    return create(isotropic,normalizedT,cross(isotropic.N,normalizedT));
  }
  static Anisotropic<RayDirInfo> create(const Isotropic<RayDirInfo> isotropic)
  {
    float2x3 TB = math::frisvad(isotropic.N);
    return create(isotropic,TB[0],TB[1]);
  }

  float3 getTangentSpaceV() {return float3(TdotV,BdotV,Isotropic<RayDirInfo>::NdotV);}
  // WARNING: its the transpose of the old GLSL function return value!
  float3x3 getTangentFrame() {return float3x3(T,B,Isotropic<RayDirInfo>::N);}

  float3 T;
  float3 B;
  float3 TdotV;
  float3 BdotV;
};

}


template<class RayDirInfo>
struct LightSample
{
  static LightSample<RayDirInfo> createTangentSpace(
    const float3 tangentSpaceV,
    const RayDirInfo tangentSpaceL,
    const float3x3 tangentFrame // WARNING: its the transpose of the old GLSL function return value!
  )
  {
    LightSample<RayDirInfo> retval;
    
    retval.L = RayDirInfo::transform(tangentSpaceL,tangentFrame);
    retval.VdotL = dot(tangentSpaceV,tangentSpaceL);

    retval.TdotL = tangentSpaceL.x;
    retval.BdotL = tangentSpaceL.y;
    retval.NdotL = tangentSpaceL.z;
    retval.NdotL2 = retval.NdotL*retval.NdotL;
    
    return retval;
  }
  static LightSample<RayDirInfo> create(const RayDirInfo L, const float VdotL, const float3 N)
  {
    LightSample<RayDirInfo> retval;
    
    retval.L = L;
    retval.VdotL = VdotL;

    retval.TdotL = nbl::hlsl::numeric_limits<float>::nan();
    retval.BdotL = nbl::hlsl::numeric_limits<float>::nan();
    retval.NdotL = dot(N,L);
    retval.NdotL2 = retval.NdotL*retval.NdotL;
    
    return retval;
  }
  static LightSample<RayDirInfo> create(const RayDirInfo L, const float VdotL, const float3 T, const float3 B, const float3 N)
  {
    LightSample<RayDirInfo> retval = create(L,VdotL,N);
    
    retval.TdotL = dot(T,L);
    retval.BdotL = dot(B,L);
    
    return retval;
  }
  // overloads for surface_interactions
  template<class ObserverRayDirInfo>
  static LightSample<RayDirInfo> create(const float3 L, const surface_interactions::Isotropic<ObserverRayDirInfo> interaction)
  {
    const float3 V = interaction.V.getDirection();
    const float VdotL = dot(V,L);
    return create(L,VdotL,interaction.N);
  }
  template<class ObserverRayDirInfo>
  static LightSample<RayDirInfo> create(const float3 L, const surface_interactions::Anisotropic<ObserverRayDirInfo> interaction)
  {
    const float3 V = interaction.V.getDirection();
    const float VdotL = dot(V,L);
    return create(L,VdotL,interaction.T,interaction.B,interaction.N);
  }
  //
  float3 getTangentSpaceL()
  {
    return float3(TdotL,BdotL,NdotL);
  }

  RayDirInfo L;
  float VdotL;

  float TdotL; 
  float BdotL;
  float NdotL;
  float NdotL2;
};


//
struct IsotropicMicrofacetCache
{
  // always valid because its specialized for the reflective case
  static IsotropicMicrofacetCache createForReflection(const float NdotV, const float NdotL, const float VdotL, out float LplusV_rcpLen)
  {
    LplusV_rcpLen = rsqrt(2.0+2.0*VdotL);

    IsotropicMicrofacetCache retval;
    
    retval.VdotH = LplusV_rcpLen*VdotL+LplusV_rcpLen;
    retval.LdotH = retval.VdotH;
    retval.NdotH = (NdotL+NdotV)*LplusV_rcpLen;
    retval.NdotH2 = retval.NdotH*retval.NdotH;
    
    return retval;
  }
  static IsotropicMicrofacetCache createForReflection(const float NdotV, const float NdotL, const float VdotL)
  {
    float dummy;
    return createForReflection(NdotV,NdotL,VdotL,dummy);
  }
  template<class ObserverRayDirInfo, class IncomingRayDirInfo>
  static IsotropicMicrofacetCache createForReflection(
    const surface_interactions::Isotropic<ObserverRayDirInfo> interaction, 
    const LightSample<IncomingRayDirInfo> _sample)
  {
    return createForReflection(interaction.NdotV,_sample.NdotL,_sample.VdotL);
  }
  // transmissive cases need to be checked if the path is valid before usage
  static bool compute(
    out IsotropicMicrofacetCache retval,
    const bool transmitted, const float3 V, const float3 L,
    const float3 N, const float NdotL, const float VdotL,
    const float orientedEta, const float rcpOrientedEta, out float3 H
  )
  {
    // TODO: can we optimize?
    H = math::computeMicrofacetNormal(transmitted,V,L,orientedEta);
    retval.NdotH = dot(N,H);
    
    // not coming from the medium (reflected) OR
    // exiting at the macro scale AND ( (not L outside the cone of possible directions given IoR with constraint VdotH*LdotH<0.0) OR (microfacet not facing toward the macrosurface, i.e. non heightfield profile of microsurface) ) 
    const bool valid = !transmitted || (VdotL<=-min(orientedEta,rcpOrientedEta) && retval.NdotH>nbl::hlsl::numeric_limits<float>::min());
    if (valid)
    {
      // TODO: can we optimize?
      retval.VdotH = dot(V,H);
      retval.LdotH = dot(L,H);
      retval.NdotH2 = retval.NdotH*retval.NdotH;
      return true;
    }
    return false;
  }
  template<class ObserverRayDirInfo, class IncomingRayDirInfo>
  static bool compute(
    out IsotropicMicrofacetCache retval,
    const surface_interactions::Isotropic<ObserverRayDirInfo> interaction, 
    const LightSample<IncomingRayDirInfo> _sample,
    const float eta, out float3 H
  )
  {
    const float NdotV = interaction.NdotV;
    const float NdotL = _sample.NdotL;
    const bool transmitted = math::isTransmissionPath(NdotV,NdotL);
    
    float orientedEta, rcpOrientedEta;
    const bool backside = math::getOrientedEtas(orientedEta,rcpOrientedEta,NdotV,eta);

    const float3 V = interaction.V.getDirection();
    const float3 L = _sample.L;
    const float VdotL = dot(V,L);
    return nbl_glsl_calcIsotropicMicrofacetCache(retval,transmitted,V,L,interaction.N,NdotL,VdotL,orientedEta,rcpOrientedEta,H);
  }
  template<class ObserverRayDirInfo, class IncomingRayDirInfo>
  static bool compute(
    out IsotropicMicrofacetCache retval,
    const surface_interactions::Isotropic<ObserverRayDirInfo> interaction, 
    const LightSample<IncomingRayDirInfo> _sample,
    const float eta
  )
  {
    float3 dummy;
    //float3 V = interaction.V.getDirection();
    //return nbl_glsl_calcIsotropicMicrofacetCache(retval,transmitted,V,L,interaction.N,NdotL,_sample.VdotL,orientedEta,rcpOrientedEta,dummy);
    // TODO try to use other function that doesnt need `dummy`
    //  (computing H that is then lost)
    return compute(retval, interaction, _sample, eta, dummy);
  }

  bool isValidVNDFMicrofacet(const bool is_bsdf, const bool transmission, const float VdotL, const float eta, const float rcp_eta)
  {
    return NdotH >= 0.0 && !(is_bsdf && transmission && (VdotL > -min(eta,rcp_eta)));
  }

  float VdotH;
  float LdotH;
  float NdotH;
  float NdotH2;
};

struct AnisotropicMicrofacetCache : IsotropicMicrofacetCache
{
  // always valid by construction
  static AnisotropicMicrofacetCache create(const float3 tangentSpaceV, const float3 tangentSpaceH)
  {
    AnisotropicMicrofacetCache retval;
    
    retval.VdotH = dot(tangentSpaceV,tangentSpaceH);
    retval.LdotH = retval.VdotH;
    retval.NdotH = tangentSpaceH.z;
    retval.NdotH2 = retval.NdotH*retval.NdotH;
    retval.TdotH = tangentSpaceH.x;
    retval.BdotH = tangentSpaceH.y;
    
    return retval;
  }
  static AnisotropicMicrofacetCache create(
    const float3 tangentSpaceV, 
    const float3 tangentSpaceH,
    const bool transmitted,
    const float rcpOrientedEta,
    const float rcpOrientedEta2
  )
  {
    AnisotropicMicrofacetCache retval = create(tangentSpaceV,tangentSpaceH);
    if (transmitted)
    {
      const float VdotH = retval.VdotH;
      retval.LdotH = transmitted ? math::refract_compute_NdotT(VdotH<0.0,VdotH*VdotH,rcpOrientedEta2) : VdotH;
    }
    
    return retval;
  }
  // always valid because its specialized for the reflective case
  static AnisotropicMicrofacetCache createForReflection(const float3 tangentSpaceV, const float3 tangentSpaceL, const float VdotL)
  {
    AnisotropicMicrofacetCache retval;
    
    float LplusV_rcpLen;
    (IsotropicMicrofacetCache) retval = IsotropicMicrofacetCache::createForReflection(tangentSpaceV.z,tangentSpaceL.z,VdotL,LplusV_rcpLen);
    retval.TdotH = (tangentSpaceV.x+tangentSpaceL.x)*LplusV_rcpLen;
    retval.BdotH = (tangentSpaceV.y+tangentSpaceL.y)*LplusV_rcpLen;
    
    return retval;
  }
  template<class ObserverRayDirInfo, class IncomingRayDirInfo>
  static AnisotropicMicrofacetCache createForReflection(
    const surface_interactions::Anisotropic<ObserverRayDirInfo> interaction, 
    const LightSample<IncomingRayDirInfo> _sample)
  {
    return createForReflection(interaction.getTangentSpaceV(),_sample.getTangentSpaceL(),_sample.VdotL);
  }
  // transmissive cases need to be checked if the path is valid before usage
  static bool compute(
    out AnisotropicMicrofacetCache retval,
    const bool transmitted, const float3 V, const float3 L,
    const float3 T, const float3 B, const float3 N,
    const float NdotL, const float VdotL,
    const float orientedEta, const float rcpOrientedEta, out float3 H
  )
  {
    const bool valid = IsotropicMicrofacetCache::compute(retval,transmitted,V,L,N,NdotL,VdotL,orientedEta,rcpOrientedEta,H);
    if (valid)
    {
      retval.TdotH = dot(T,H);
      retval.BdotH = dot(B,H);
    }
    return valid;
  }
  template<class ObserverRayDirInfo, class IncomingRayDirInfo>
  static bool compute(
    out AnisotropicMicrofacetCache retval,
    const surface_interactions::Anisotropic<ObserverRayDirInfo> interaction, 
    const LightSample<IncomingRayDirInfo> _sample,
    const float eta
  )
  {
    float3 H;
    const bool valid = IsotropicMicrofacetCache::compute(retval,interaction,_sample,eta,H);
    if (valid)
    {
      retval.TdotH = dot(interaction.T,H);
      retval.BdotH = dot(interaction.B,H);
    }
    return valid;
  }

  float TdotH;
  float BdotH;
};


// finally fixed the semantic F-up, value/pdf = quotient not remainder
template<typename SpectralBins>
struct quotient_and_pdf
{
  static quotient_and_pdf<SpectralBins> create(const SpectralBins _quotient, const float _pdf)
  {
    quotient_and_pdf<SpectralBins> retval;
    retval.quotient = _quotient;
    retval.pdf = _pdf;
    return retval;
  }

  SpectralBins value()
  {
    return quotient*pdf;
  }
  
  SpectralBins quotient;
  float pdf;
};

using quotient_and_pdf_scalar = quotient_and_pdf<float>;
using quotient_and_pdf_rgb = quotient_and_pdf<float3>;


// Utility class
// Making it easier to conform your BxDFs to the concept
template <class Spectrum, class Pdf, class Sample, class Interaction, class MicrofacetCache>
struct BxDFBase
{
    // BxDFs must define such typenames:
    using spectrum_t =      Spectrum;
    using pdf_t =           Pdf;
    using q_pdf_t =         quotient_and_pdf<spectrum_t>;

    using sample_t =        Sample;
    using interaction_t =   Interaction;
    using cache_t =         MicrofacetCache;

    /**
    * BxDFs (i.e. types derived from this base) must define following member functions:
    * 
    * spectrum_t    eval(in sample_t, in interaction_t, in cache_t);
    * 
    * sample_t      generate(in interaction_t, inout float3 u);
    * 
    * q_pdf_t       quotient_and_pdf(in sample_t, in interaction_t, in cache_t);
    */
};

}
}
}

#endif
