#define _NBL_VIDEO_I_GPU_ACCELERATION_STRUCTURE_CPP_
#include "nbl/video/ILogicalDevice.h"
#include "nbl/video/IPhysicalDevice.h"

namespace nbl::video
{

template<class BufferType>
bool IGPUAccelerationStructure::BuildInfo<BufferType>::invalid(const IGPUAccelerationStructure* const src, const IGPUAccelerationStructure* const dst)
{
    if (!dst)
        return true;

	if (isUpdate)
	{
		// https://registry.khronos.org/vulkan/specs/1.3-extensions/html/vkspec.html#VUID-vkBuildAccelerationStructuresKHR-srcAccelerationStructure-04629
		// https://registry.khronos.org/vulkan/specs/1.3-extensions/html/vkspec.html#VUID-vkBuildAccelerationStructuresKHR-pInfos-04630
		if (!src || src->getOriginDevice()!=device)
			return true;
	}

    const auto device = dst->getOriginDevice();
	if constexpr (std::is_same_v<BufferType,asset::ICPUBuffer>)
	{
		// https://registry.khronos.org/vulkan/specs/1.3-extensions/html/vkspec.html#VUID-vkBuildAccelerationStructuresKHR-accelerationStructureHostCommands-03581
		if (!device->getEnabledFeatures().accelerationStructureHostCommands)
			return true;
		// https://registry.khronos.org/vulkan/specs/1.3-extensions/html/vkspec.html#VUID-vkBuildAccelerationStructuresKHR-pInfos-03722
		if (device->invalidAccelerationStructureForHostOperations(dst))
			return true;
		// https://registry.khronos.org/vulkan/specs/1.3-extensions/html/vkspec.html#VUID-vkBuildAccelerationStructuresKHR-pInfos-03723
		if (isUpdate && device->invalidAccelerationStructureForHostOperations(src))
			return true;
	}

    return false;
}


template<class BufferType>
inline uint32_t IGPUBottomLevelAccelerationStructure::BuildInfo<BufferType>::valid(const BuildRangeInfo* const buildRangeInfos) const
{
	if (IGPUAccelerationStructure::BuildInfo<BufferType>::invalid(srcAS,dstAS))
		return false;

	// destination and scratch
	uint32_t retval = 2u;
	if (isUpdate) // source
		retval++;

	#ifdef _NBL_DEBUG
	const auto& bufferUsages = getOriginDevice()->getPhysicalDevice()->getBufferFormatUsages();
	for (auto i=0u; i<geometries.size(); i++)
	{
		const auto& geometry = geometries[i];
		if (geometry.isAABB)
		{
			// https://registry.khronos.org/vulkan/specs/1.3-extensions/html/vkspec.html#VUID-vkBuildAccelerationStructuresKHR-pInfos-03774
			if (!geometry.aabbs.data.isValid())
				return false;
		}
		else
		{
			if (!bufferUsages[geometry.triangles.vertexFormat].accelerationStructureVertex)
				return false;
			// https://registry.khronos.org/vulkan/specs/1.3-extensions/html/vkspec.html#VUID-vkBuildAccelerationStructuresKHR-pInfos-03771
			if (!geometry.indexData.isValid())
				return false;
			// https://registry.khronos.org/vulkan/specs/1.3-extensions/html/vkspec.html#VUID-vkBuildAccelerationStructuresKHR-pInfos-03772
			if (geometry.indexType!=asset::EIT_UNKNOWN && !geometry.indexData.isValid())
				return false;
			// https://registry.khronos.org/vulkan/specs/1.3-extensions/html/vkspec.html#VUID-vkBuildAccelerationStructuresKHR-pInfos-03773
			if (geometry.transformData.buffer && !geometry.transformData.isValid())
				return false;
		}
	}
	#endif

	constexpr uint32_t MaxBuffersPerGeometry = dstAS->getCreationFlags().hasFlags(IAccelerationStructure::CREATE_FLAGS::MOTION_BIT) ? 4u:3u;
	retval += MaxBuffersPerGeometry*geometries.size();
	return retval;
}

template<class BufferType>
inline uint32_t IGPUTopLevelAccelerationStructure::BuildInfo<BufferType>::valid(const uint32_t* const instanceCounts)
{
	if (IGPUAccelerationStructure::BuildInfo<BufferType>::invalid(srcAS,dstAS))
		return false;

	// destination and scratch
	uint32_t retval = 2u;
	if (isUpdate) // source
		retval++;

	constexpr bool HostBuild = std::is_same_v<BufferType,asset::ICPUBuffer>;

	const auto device = dst->getOriginDevice();
	const auto physDev = device->getPhysicalDevice();
	for (auto i=0u; i<geometries.size(); i++)
	{
		const auto& geometry = geometries[i];
		// https://registry.khronos.org/vulkan/specs/1.3-extensions/html/vkspec.html#VUID-vkBuildAccelerationStructuresKHR-pInfos-03778
		if (!geometry.instanceData.isValid())
			return false;
		// https://registry.khronos.org/vulkan/specs/1.3-extensions/html/vkspec.html#VUID-vkBuildAccelerationStructuresKHR-pInfos-03801
		if (instanceCounts && instanceCounts[i]>physDev->getLimits().maxAccelerationStructureInstanceCount) // TODO: review
			return false;
		#ifdef _NBL_DEBUG
		/* TODO: with `EXT_private_data
		// https://registry.khronos.org/vulkan/specs/1.3-extensions/html/vkspec.html#VUID-vkBuildAccelerationStructuresKHR-pInfos-03724
		if constexpr (HostBuild)
		{
			for (auto 
			if (device->invalidAccelerationStructureForHostOperations(getAccelerationStructureFromReference(geometry.instanceData.blas)))
				return false;
		}
		*/
		#endif
	}

	if (totalInstanceCount>m_maxInstanceCount)
		return false;

	retval += geometries.size();
	return retval;
}

}