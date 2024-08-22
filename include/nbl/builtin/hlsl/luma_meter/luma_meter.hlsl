// Copyright (C) 2018-2024 - DevSH Graphics Programming Sp. z O.O.
// This file is part of the "Nabla Engine".
// For conditions of distribution and use, see copyright notice in nabla.h

#ifndef _NBL_BUILTIN_HLSL_LUMA_METER_INCLUDED_
#define _NBL_BUILTIN_HLSL_LUMA_METER_INCLUDED_

#include "nbl/builtin/hlsl/glsl_compat/core.hlsl"
#include "nbl/builtin/hlsl/glsl_compat/subgroup_basic.hlsl"
#include "nbl/builtin/hlsl/workgroup/basic.hlsl"
#include "nbl/builtin/hlsl/workgroup/arithmetic.hlsl"
#include "nbl/builtin/hlsl/type_traits.hlsl"
#include "nbl/builtin/hlsl/math/morton.hlsl"
#include "nbl/builtin/hlsl/luma_meter/common.hlsl"

namespace nbl
{
namespace hlsl
{
namespace luma_meter
{

template<uint32_t GroupSize, typename ValueAccessor, typename SharedAccessor, typename TexAccessor>
struct geom_meter {
    using float_t = typename SharedAccessor::type;
    using float_t2 = typename conditional<is_same_v<float_t, float32_t>, float32_t2, float16_t2>::type;
    using float_t3 = typename conditional<is_same_v<float_t, float32_t>, float32_t3, float16_t3>::type;
    using this_t = geom_meter<GroupSize, ValueAccessor, SharedAccessor, TexAccessor>;

    static this_t create(float_t2 lumaMinMax)
    {
        this_t retval;
        retval.lumaMinMax = lumaMinMax;
        return retval;
    }

    float_t reduction(float_t value, NBL_REF_ARG(SharedAccessor) sdata)
    {
        return workgroup::reduction < plus < float_t >, GroupSize >::
            template __call <SharedAccessor>(value, sdata);
    }

    float_t computeLumaLog2(
        NBL_CONST_REF_ARG(MeteringWindow) window,
        NBL_REF_ARG(TexAccessor) tex,
        float_t2 shiftedCoord
    )
    {
        float_t2 uvPos = shiftedCoord * window.meteringWindowScale + window.meteringWindowOffset;
        float_t3 color = tex.get(uvPos);
        float_t luma = TexAccessor::toXYZ(color);

        luma = clamp(luma, lumaMinMax.x, lumaMinMax.y);

        return max(log2(luma), log2(lumaMinMax.x));
    }

    void sampleLuma(
        NBL_CONST_REF_ARG(MeteringWindow) window,
        NBL_REF_ARG(ValueAccessor) val,
        NBL_REF_ARG(TexAccessor) tex,
        NBL_REF_ARG(SharedAccessor) sdata,
        float_t2 tileOffset
    )
    {
        uint32_t tid = workgroup::SubgroupContiguousIndex();
        uint32_t2 coord = {
            morton2d_decode_x(tid),
            morton2d_decode_y(tid)
        };

        float_t luma = 0.0f;
        luma = computeLumaLog2(window, tex, tileOffset + (float32_t2)(coord));
        float_t lumaSum = reduction(luma, sdata);

        if (tid == GroupSize - 1) {
            uint32_t3 workGroupCount = glsl::gl_NumWorkGroups();
            uint32_t fixedPointBitsLeft = 32 - uint32_t(ceil(log2(workGroupCount.x * workGroupCount.y * workGroupCount.z))) + glsl::gl_SubgroupSizeLog2();

            uint32_t lumaSumBitPattern = uint32_t(clamp((lumaSum - log2(lumaMinMax.x)) * (log2(lumaMinMax.y) - log2(lumaMinMax.x)), 0.f, float32_t((1 << fixedPointBitsLeft) - 1)));

            uint32_t3 workgroupSize = glsl::gl_WorkGroupSize();
            uint32_t workgroupIndex = dot(uint32_t3(workgroupSize.y * workgroupSize.z, workgroupSize.z, 1), glsl::gl_WorkGroupID());

            val.atomicAdd(workgroupIndex & ((1 << glsl::gl_SubgroupSizeLog2()) - 1), lumaSumBitPattern);
        }
    }

    float_t2 lumaMinMax;
};
}
}
}

#endif