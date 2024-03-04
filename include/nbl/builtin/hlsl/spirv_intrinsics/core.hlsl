// Copyright (C) 2023 - DevSH Graphics Programming Sp. z O.O.
// This file is part of the "Nabla Engine".
// For conditions of distribution and use, see copyright notice in nabla.h
#ifndef _NBL_BUILTIN_HLSL_SPIRV_INTRINSICS_CORE_INCLUDED_
#define _NBL_BUILTIN_HLSL_SPIRV_INTRINSICS_CORE_INCLUDED_


#ifdef __HLSL_VERSION // TODO: AnastZIuk fix public search paths so we don't choke
#include "spirv/unified1/spirv.hpp"
//#include "spirv/unified1/spirv.hpp" TODO: also GLSL.450.std for the extended instruction set
#endif


namespace nbl 
{
namespace hlsl
{
#ifdef __HLSL_VERSION
namespace spirv
{
[[vk::ext_builtin_input(spv::BuiltInHelperInvocation)]]
static const bool HelperInvocation;

[[vk::ext_builtin_input(spv::BuiltInNumWorkgroups)]]
static const uint32_t3 NumWorkGroups;
// TODO: Doesn't work, find out why and file issue on DXC!
//[[vk::ext_builtin_input(spv::BuiltInWorkgroupSize)]]
//static const uint32_t3 WorkgroupSize;
[[vk::ext_builtin_input(spv::BuiltInWorkgroupId)]]
static const uint32_t3 WorkgroupId;
[[vk::ext_builtin_input(spv::BuiltInLocalInvocationId)]]
static const uint32_t3 LocalInvocationId;
[[vk::ext_builtin_input(spv::BuiltInGlobalInvocationId)]]
static const uint32_t3 GlobalInvocationId;
[[vk::ext_builtin_input(spv::BuiltInLocalInvocationIndex)]]
static const uint32_t LocalInvocationIndex;

template<typename T>
T atomicAdd([[vk::ext_reference]] T ptr, uint32_t memoryScope, uint32_t memorySemantics, T value);
template<>
[[vk::ext_instruction(spv::OpAtomicIAdd)]]
int32_t atomicAdd([[vk::ext_reference]] int32_t ptr, uint32_t memoryScope, uint32_t memorySemantics, int32_t value);
template<>
[[vk::ext_instruction( spv::OpAtomicIAdd )]]
uint32_t atomicAdd([[vk::ext_reference]] uint32_t ptr, uint32_t memoryScope, uint32_t memorySemantics, uint32_t value);

template<typename T>
T atomicAnd([[vk::ext_reference]] T ptr, uint32_t memoryScope, uint32_t memorySemantics, T value);
template<>
[[vk::ext_instruction( spv::OpAtomicAnd )]]
int32_t atomicAnd([[vk::ext_reference]] int32_t ptr, uint32_t memoryScope, uint32_t memorySemantics, int32_t value);
template<>
[[vk::ext_instruction( spv::OpAtomicAnd )]]
uint32_t atomicAnd([[vk::ext_reference]] uint32_t ptr, uint32_t memoryScope, uint32_t memorySemantics, uint32_t value);

template<typename T>
T atomicOr([[vk::ext_reference]] T ptr, uint32_t memoryScope, uint32_t memorySemantics, T value);
template<>
[[vk::ext_instruction( spv::OpAtomicOr )]]
int32_t atomicOr([[vk::ext_reference]] int32_t ptr, uint32_t memoryScope, uint32_t memorySemantics, int32_t value);
template<>
[[vk::ext_instruction( spv::OpAtomicOr )]]
uint32_t atomicOr([[vk::ext_reference]] uint32_t ptr, uint32_t memoryScope, uint32_t memorySemantics, uint32_t value);

template<typename T>
T atomicXor([[vk::ext_reference]] T ptr, uint32_t memoryScope, uint32_t memorySemantics, T value);
template<>
[[vk::ext_instruction( spv::OpAtomicXor )]]
int32_t atomicXor([[vk::ext_reference]] int32_t ptr, uint32_t memoryScope, uint32_t memorySemantics, int32_t value);
template<>
[[vk::ext_instruction( spv::OpAtomicXor )]]
uint32_t atomicXor([[vk::ext_reference]] uint32_t ptr, uint32_t memoryScope, uint32_t memorySemantics, uint32_t value);

template<typename T>
T atomicMin([[vk::ext_reference]] T ptr, uint32_t memoryScope, uint32_t memorySemantics, T value);
template<>
[[vk::ext_instruction( spv::OpAtomicSMin )]]
int32_t atomicMin([[vk::ext_reference]] int32_t ptr, uint32_t memoryScope, uint32_t memorySemantics, int32_t value);
template<>
[[vk::ext_instruction( spv::OpAtomicSMin )]]
uint32_t atomicMin([[vk::ext_reference]] uint32_t ptr, uint32_t memoryScope, uint32_t memorySemantics, uint32_t value);

template<typename T>
T atomicMax([[vk::ext_reference]] T ptr, uint32_t memoryScope, uint32_t memorySemantics, T value);
template<>
[[vk::ext_instruction( spv::OpAtomicSMax )]]
int32_t atomicMax([[vk::ext_reference]] int32_t ptr, uint32_t memoryScope, uint32_t memorySemantics, int32_t value);
template<>
[[vk::ext_instruction( spv::OpAtomicSMax )]]
uint32_t atomicMax([[vk::ext_reference]] uint32_t ptr, uint32_t memoryScope, uint32_t memorySemantics, uint32_t value);

template<typename T>
T atomicExchange([[vk::ext_reference]] T ptr, uint32_t memoryScope, uint32_t memorySemantics, T value);
template<>
[[vk::ext_instruction( spv::OpAtomicExchange )]]
int32_t atomicExchange([[vk::ext_reference]] int32_t ptr, uint32_t memoryScope, uint32_t memorySemantics, int32_t value);
template<>
[[vk::ext_instruction( spv::OpAtomicExchange )]]
uint32_t atomicExchange([[vk::ext_reference]] uint32_t ptr, uint32_t memoryScope, uint32_t memorySemantics, uint32_t value);
template<>
[[vk::ext_instruction( spv::OpAtomicExchange )]]
float32_t atomicExchange([[vk::ext_reference]] float32_t ptr, uint32_t memoryScope, uint32_t memorySemantics, float32_t value);


template<typename T>
T atomicCompSwap([[vk::ext_reference]] T ptr, uint32_t memoryScope, uint32_t memSemanticsEqual, uint32_t memSemanticsUnequal, T value, T comparator);
template<>
[[vk::ext_instruction( spv::OpAtomicCompareExchange )]]
int32_t atomicCompSwap([[vk::ext_reference]] int32_t ptr, uint32_t memoryScope, uint32_t memSemanticsEqual, uint32_t memSemanticsUnequal, int32_t value, int32_t comparator);
template<>
[[vk::ext_instruction( spv::OpAtomicCompareExchange )]]
uint32_t atomicCompSwap([[vk::ext_reference]] uint32_t ptr, uint32_t memoryScope, uint32_t memSemanticsEqual, uint32_t memSemanticsUnequal, uint32_t value, uint32_t comparator);



// Memory Semantics link here: https://registry.khronos.org/SPIR-V/specs/unified1/SPIRV.html#Memory_Semantics_-id-

// https://registry.khronos.org/SPIR-V/specs/unified1/SPIRV.html#_memory_semantics_id
// By providing memory semantics None we do both control and memory barrier as is done in GLSL

[[vk::ext_instruction( spv::OpControlBarrier )]]
void controlBarrier(uint32_t executionScope, uint32_t memoryScope, uint32_t memorySemantics);

[[vk::ext_instruction( spv::OpMemoryBarrier )]]
void memoryBarrier(uint32_t memoryScope, uint32_t memorySemantics);


// Add specializations if you need to emit a `ext_capability` (this means that the instruction needs to forward through an `impl::` struct and so on)
template<class T, class U>
[[vk::ext_instruction(spv::OpBitcast)]]
T bitcast(U);

template<typename Unsigned>
[[vk::ext_instruction( spv::OpBitFieldUExtract )]]
Unsigned bitFieldUExtract( Unsigned val, uint32_t offsetBits, uint32_t numBits );

template<typename Signed>
[[vk::ext_instruction( spv::OpBitFieldSExtract )]]
Signed bitFieldSExtract( Signed val, uint32_t offsetBits, uint32_t numBits );

}
#endif
}
}

#endif