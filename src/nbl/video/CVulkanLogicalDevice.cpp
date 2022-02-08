#include "CVulkanLogicalDevice.h"

#include "nbl/video/CVulkanPhysicalDevice.h"
#include "nbl/video/CVulkanQueryPool.h"
#include "nbl/video/CVulkanCommandBuffer.h"
#include "nbl/video/CVulkanEvent.h"
#include "nbl/video/surface/CSurfaceVulkan.h"

namespace nbl::video
{

core::smart_refctd_ptr<ISwapchain> CVulkanLogicalDevice::createSwapchain(ISwapchain::SCreationParams&& params)
{
    constexpr uint32_t MAX_SWAPCHAIN_IMAGE_COUNT = 100u;

    if (params.surface->getAPIType() != EAT_VULKAN)
        return nullptr;

#ifdef _NBL_PLATFORM_WINDOWS_
    // Todo(achal): not sure yet, how would I handle multiple platforms without making
    // this function templated
    VkSurfaceKHR vk_surface = static_cast<const CSurfaceVulkanWin32*>(params.surface.get())->getInternalObject();
#endif

    VkPresentModeKHR vkPresentMode;
    if((params.presentMode & ISurface::E_PRESENT_MODE::EPM_IMMEDIATE) == ISurface::E_PRESENT_MODE::EPM_IMMEDIATE)
        vkPresentMode = VK_PRESENT_MODE_IMMEDIATE_KHR;
    else if((params.presentMode & ISurface::E_PRESENT_MODE::EPM_MAILBOX) == ISurface::E_PRESENT_MODE::EPM_MAILBOX)
        vkPresentMode = VK_PRESENT_MODE_MAILBOX_KHR;
    else if((params.presentMode & ISurface::E_PRESENT_MODE::EPM_FIFO) == ISurface::E_PRESENT_MODE::EPM_FIFO)
        vkPresentMode = VK_PRESENT_MODE_FIFO_KHR;
    else if((params.presentMode & ISurface::E_PRESENT_MODE::EPM_FIFO_RELAXED) == ISurface::E_PRESENT_MODE::EPM_FIFO_RELAXED)
        vkPresentMode = VK_PRESENT_MODE_FIFO_RELAXED_KHR;

    VkSwapchainCreateInfoKHR vk_createInfo = { VK_STRUCTURE_TYPE_SWAPCHAIN_CREATE_INFO_KHR };
#ifdef _NBL_PLATFORM_WINDOWS_
    vk_createInfo.surface = vk_surface;
#endif        
    vk_createInfo.minImageCount = params.minImageCount;
    vk_createInfo.imageFormat = getVkFormatFromFormat(params.surfaceFormat.format);
    vk_createInfo.imageColorSpace = getVkColorSpaceKHRFromColorSpace(params.surfaceFormat.colorSpace);
    vk_createInfo.imageExtent = { params.width, params.height };
    vk_createInfo.imageArrayLayers = params.arrayLayers;
    vk_createInfo.imageUsage = static_cast<VkImageUsageFlags>(params.imageUsage);
    vk_createInfo.imageSharingMode = static_cast<VkSharingMode>(params.imageSharingMode);
    vk_createInfo.queueFamilyIndexCount = params.queueFamilyIndexCount;
    vk_createInfo.pQueueFamilyIndices = params.queueFamilyIndices;
    vk_createInfo.preTransform = static_cast<VkSurfaceTransformFlagBitsKHR>(params.preTransform);
    vk_createInfo.compositeAlpha = static_cast<VkCompositeAlphaFlagBitsKHR>(params.compositeAlpha);
    vk_createInfo.presentMode = vkPresentMode;
    vk_createInfo.clipped = VK_FALSE;
    vk_createInfo.oldSwapchain = VK_NULL_HANDLE;
    if (params.oldSwapchain && (params.oldSwapchain->getAPIType() == EAT_VULKAN))
        vk_createInfo.oldSwapchain = IBackendObject::device_compatibility_cast<CVulkanSwapchain*>(params.oldSwapchain.get(), this)->getInternalObject();

    VkSwapchainKHR vk_swapchain;
    if (m_devf.vk.vkCreateSwapchainKHR(m_vkdev, &vk_createInfo, nullptr, &vk_swapchain) != VK_SUCCESS)
        return nullptr;

    uint32_t imageCount;
    VkResult retval = m_devf.vk.vkGetSwapchainImagesKHR(m_vkdev, vk_swapchain, &imageCount, nullptr);
    if ((retval != VK_SUCCESS) && (retval != VK_INCOMPLETE))
        return nullptr;

    assert(imageCount <= MAX_SWAPCHAIN_IMAGE_COUNT);

    VkImage vk_images[MAX_SWAPCHAIN_IMAGE_COUNT];
    retval = m_devf.vk.vkGetSwapchainImagesKHR(m_vkdev, vk_swapchain, &imageCount, vk_images);
    if ((retval != VK_SUCCESS) && (retval != VK_INCOMPLETE))
        return nullptr;

    ISwapchain::images_array_t images = core::make_refctd_dynamic_array<ISwapchain::images_array_t>(imageCount);

    uint32_t i = 0u;
    for (auto& image : (*images))
    {
        CVulkanForeignImage::SCreationParams creationParams;
        creationParams.arrayLayers = params.arrayLayers;
        creationParams.extent = { params.width, params.height, 1u };
        creationParams.flags = static_cast<CVulkanForeignImage::E_CREATE_FLAGS>(0); // Todo(achal)
        creationParams.format = params.surfaceFormat.format;
        creationParams.mipLevels = 1u;
        creationParams.samples = CVulkanImage::ESCF_1_BIT; // Todo(achal)
        creationParams.type = CVulkanImage::ET_2D;

        image = core::make_smart_refctd_ptr<CVulkanForeignImage>(
            core::smart_refctd_ptr<CVulkanLogicalDevice>(this), std::move(creationParams),
            vk_images[i++]);
    }

    return core::make_smart_refctd_ptr<CVulkanSwapchain>(
        core::smart_refctd_ptr<CVulkanLogicalDevice>(this), std::move(params),
        std::move(images), vk_swapchain);
}

core::smart_refctd_ptr<IGPUEvent> CVulkanLogicalDevice::createEvent(IGPUEvent::E_CREATE_FLAGS flags)
{
    VkEventCreateInfo vk_createInfo = { VK_STRUCTURE_TYPE_EVENT_CREATE_INFO };
    vk_createInfo.pNext = nullptr;
    vk_createInfo.flags = static_cast<VkEventCreateFlags>(flags);

    VkEvent vk_event;
    if (m_devf.vk.vkCreateEvent(m_vkdev, &vk_createInfo, nullptr, &vk_event) == VK_SUCCESS)
        return core::make_smart_refctd_ptr<CVulkanEvent>(
            core::smart_refctd_ptr<const CVulkanLogicalDevice>(this), flags, vk_event);
    else
        return nullptr;
};

IGPUEvent::E_STATUS CVulkanLogicalDevice::getEventStatus(const IGPUEvent* _event)
{
    if (!_event || _event->getAPIType() != EAT_VULKAN)
        return IGPUEvent::E_STATUS::ES_FAILURE;

    VkEvent vk_event = IBackendObject::device_compatibility_cast<const CVulkanEvent*>(_event, this)->getInternalObject();
    VkResult retval = m_devf.vk.vkGetEventStatus(m_vkdev, vk_event);
    switch (retval)
    {
    case VK_EVENT_SET:
        return IGPUEvent::ES_SET;
    case VK_EVENT_RESET:
        return IGPUEvent::ES_RESET;
    default:
        return IGPUEvent::ES_FAILURE;
    }
}

IGPUEvent::E_STATUS CVulkanLogicalDevice::resetEvent(IGPUEvent* _event)
{
    if (!_event || _event->getAPIType() != EAT_VULKAN)
        return IGPUEvent::E_STATUS::ES_FAILURE;

    VkEvent vk_event = IBackendObject::device_compatibility_cast<const CVulkanEvent*>(_event, this)->getInternalObject();
    if (m_devf.vk.vkResetEvent(m_vkdev, vk_event) == VK_SUCCESS)
        return IGPUEvent::ES_RESET;
    else
        return IGPUEvent::ES_FAILURE;
}

IGPUEvent::E_STATUS CVulkanLogicalDevice::setEvent(IGPUEvent* _event)
{
    if (!_event || _event->getAPIType() != EAT_VULKAN)
        return IGPUEvent::E_STATUS::ES_FAILURE;

    VkEvent vk_event = IBackendObject::device_compatibility_cast<const CVulkanEvent*>(_event, this)->getInternalObject();
    if (m_devf.vk.vkSetEvent(m_vkdev, vk_event) == VK_SUCCESS)
        return IGPUEvent::ES_SET;
    else
        return IGPUEvent::ES_FAILURE;
}

core::smart_refctd_ptr<IDriverMemoryAllocation> CVulkanLogicalDevice::allocateDeviceLocalMemory(
    const IDriverMemoryBacked::SDriverMemoryRequirements& additionalReqs)
{
    IDriverMemoryBacked::SDriverMemoryRequirements memoryReqs = getDeviceLocalGPUMemoryReqs();
    memoryReqs.vulkanReqs.alignment = core::max(memoryReqs.vulkanReqs.alignment, additionalReqs.vulkanReqs.alignment);
    memoryReqs.vulkanReqs.size = core::max(memoryReqs.vulkanReqs.size, additionalReqs.vulkanReqs.size);
    memoryReqs.vulkanReqs.memoryTypeBits &= additionalReqs.vulkanReqs.memoryTypeBits;
    memoryReqs.mappingCapability = additionalReqs.mappingCapability;
    memoryReqs.memoryHeapLocation = IDriverMemoryAllocation::ESMT_DEVICE_LOCAL;
    memoryReqs.prefersDedicatedAllocation = additionalReqs.prefersDedicatedAllocation;
    memoryReqs.requiresDedicatedAllocation = additionalReqs.requiresDedicatedAllocation;

    return allocateGPUMemory(additionalReqs);
}

core::smart_refctd_ptr<IDriverMemoryAllocation> CVulkanLogicalDevice::allocateSpilloverMemory(
    const IDriverMemoryBacked::SDriverMemoryRequirements& additionalReqs)
{
    if (additionalReqs.memoryHeapLocation == IDriverMemoryAllocation::ESMT_DEVICE_LOCAL)
        return nullptr;

    IDriverMemoryBacked::SDriverMemoryRequirements memoryReqs = getSpilloverGPUMemoryReqs();
    memoryReqs.vulkanReqs.alignment = core::max(memoryReqs.vulkanReqs.alignment, additionalReqs.vulkanReqs.alignment);
    memoryReqs.vulkanReqs.size = core::max(memoryReqs.vulkanReqs.size, additionalReqs.vulkanReqs.size);
    memoryReqs.vulkanReqs.memoryTypeBits &= additionalReqs.vulkanReqs.memoryTypeBits;
    memoryReqs.mappingCapability = additionalReqs.mappingCapability;
    memoryReqs.memoryHeapLocation = additionalReqs.memoryHeapLocation;
    memoryReqs.prefersDedicatedAllocation = additionalReqs.prefersDedicatedAllocation;
    memoryReqs.requiresDedicatedAllocation = additionalReqs.requiresDedicatedAllocation;

    return allocateGPUMemory(memoryReqs);
}

core::smart_refctd_ptr<IDriverMemoryAllocation> CVulkanLogicalDevice::allocateUpStreamingMemory(
    const IDriverMemoryBacked::SDriverMemoryRequirements& additionalReqs)
{
    IDriverMemoryBacked::SDriverMemoryRequirements memoryReqs = getUpStreamingMemoryReqs();
    memoryReqs.vulkanReqs.alignment = core::max(memoryReqs.vulkanReqs.alignment, additionalReqs.vulkanReqs.alignment);
    memoryReqs.vulkanReqs.size = core::max(memoryReqs.vulkanReqs.size, additionalReqs.vulkanReqs.size);
    memoryReqs.vulkanReqs.memoryTypeBits &= additionalReqs.vulkanReqs.memoryTypeBits;
    memoryReqs.mappingCapability = additionalReqs.mappingCapability;
    memoryReqs.memoryHeapLocation = additionalReqs.memoryHeapLocation;
    memoryReqs.prefersDedicatedAllocation = additionalReqs.prefersDedicatedAllocation;
    memoryReqs.requiresDedicatedAllocation = additionalReqs.requiresDedicatedAllocation;

    return allocateGPUMemory(memoryReqs);
}

core::smart_refctd_ptr<IDriverMemoryAllocation> CVulkanLogicalDevice::allocateDownStreamingMemory(
    const IDriverMemoryBacked::SDriverMemoryRequirements& additionalReqs)
{
    IDriverMemoryBacked::SDriverMemoryRequirements memoryReqs = getDownStreamingMemoryReqs();
    memoryReqs.vulkanReqs.alignment = core::max(memoryReqs.vulkanReqs.alignment, additionalReqs.vulkanReqs.alignment);
    memoryReqs.vulkanReqs.size = core::max(memoryReqs.vulkanReqs.size, additionalReqs.vulkanReqs.size);
    memoryReqs.vulkanReqs.memoryTypeBits &= additionalReqs.vulkanReqs.memoryTypeBits;
    memoryReqs.mappingCapability = additionalReqs.mappingCapability;
    memoryReqs.memoryHeapLocation = additionalReqs.memoryHeapLocation;
    memoryReqs.prefersDedicatedAllocation = additionalReqs.prefersDedicatedAllocation;
    memoryReqs.requiresDedicatedAllocation = additionalReqs.requiresDedicatedAllocation;

    return allocateGPUMemory(memoryReqs);
}

core::smart_refctd_ptr<IDriverMemoryAllocation> CVulkanLogicalDevice::allocateCPUSideGPUVisibleMemory(
    const IDriverMemoryBacked::SDriverMemoryRequirements& additionalReqs)
{
    if (additionalReqs.memoryHeapLocation != IDriverMemoryAllocation::ESMT_NOT_DEVICE_LOCAL)
        return nullptr;

    IDriverMemoryBacked::SDriverMemoryRequirements memoryReqs = getCPUSideGPUVisibleGPUMemoryReqs();
    memoryReqs.vulkanReqs.alignment = core::max(memoryReqs.vulkanReqs.alignment, additionalReqs.vulkanReqs.alignment);
    memoryReqs.vulkanReqs.size = core::max(memoryReqs.vulkanReqs.size, additionalReqs.vulkanReqs.size);
    memoryReqs.vulkanReqs.memoryTypeBits &= additionalReqs.vulkanReqs.memoryTypeBits;
    memoryReqs.mappingCapability = additionalReqs.mappingCapability;
    memoryReqs.memoryHeapLocation = additionalReqs.memoryHeapLocation;
    memoryReqs.prefersDedicatedAllocation = additionalReqs.prefersDedicatedAllocation;
    memoryReqs.requiresDedicatedAllocation = additionalReqs.requiresDedicatedAllocation;

    return allocateGPUMemory(memoryReqs);
}

core::smart_refctd_ptr<IDriverMemoryAllocation> CVulkanLogicalDevice::allocateGPUMemory(
    const IDriverMemoryBacked::SDriverMemoryRequirements& reqs,
    core::bitflag<IDriverMemoryAllocation::E_MEMORY_ALLOCATE_FLAGS> allocateFlags)
{
    VkMemoryPropertyFlags desiredMemoryProperties = static_cast<VkMemoryPropertyFlags>(0u);

    if (reqs.memoryHeapLocation == IDriverMemoryAllocation::ESMT_DEVICE_LOCAL)
        desiredMemoryProperties |= VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT;

    if ((reqs.mappingCapability & IDriverMemoryAllocation::EMCF_CAN_MAP_FOR_READ) ||
        (reqs.mappingCapability & IDriverMemoryAllocation::EMCF_CAN_MAP_FOR_WRITE))
        desiredMemoryProperties |= VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT;

    if (reqs.mappingCapability & IDriverMemoryAllocation::EMCF_COHERENT)
        desiredMemoryProperties |= VK_MEMORY_PROPERTY_HOST_COHERENT_BIT;

    if (reqs.mappingCapability & IDriverMemoryAllocation::EMCF_CACHED)
        desiredMemoryProperties |= VK_MEMORY_PROPERTY_HOST_CACHED_BIT;

    const IPhysicalDevice::SMemoryProperties& memoryProperties = m_physicalDevice->getMemoryProperties();

    uint32_t compatibleMemoryTypeCount = 0u;
    uint32_t compatibleMemoryTypeIndices[VK_MAX_MEMORY_TYPES];

    for (uint32_t i = 0u; i < memoryProperties.memoryTypeCount; ++i)
    {
        const bool memoryTypeSupportedForResource = (reqs.vulkanReqs.memoryTypeBits & (1 << i));

        const bool memoryHasDesirableProperties = (memoryProperties.memoryTypes[i].propertyFlags
            & desiredMemoryProperties) == desiredMemoryProperties;

        if (memoryTypeSupportedForResource && memoryHasDesirableProperties)
            compatibleMemoryTypeIndices[compatibleMemoryTypeCount++] = i;
    }

    for (uint32_t i = 0u; i < compatibleMemoryTypeCount; ++i)
    {
        // Todo(achal): Make use of requiresDedicatedAllocation and prefersDedicatedAllocation
        VkMemoryAllocateFlagsInfo vk_allocateFlags = { VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_FLAGS_INFO, nullptr };
        if (allocateFlags.hasValue(IDriverMemoryAllocation::EMAF_DEVICE_MASK_BIT))
            vk_allocateFlags.flags |= VK_MEMORY_ALLOCATE_DEVICE_MASK_BIT;
        else if(allocateFlags.hasValue(IDriverMemoryAllocation::EMAF_DEVICE_ADDRESS_BIT))
            vk_allocateFlags.flags |= VK_MEMORY_ALLOCATE_DEVICE_ADDRESS_BIT;
        else if (allocateFlags.hasValue(IDriverMemoryAllocation::EMAF_DEVICE_ADDRESS_CAPTURE_REPLAY_BIT))
            vk_allocateFlags.flags |= VK_MEMORY_ALLOCATE_DEVICE_ADDRESS_CAPTURE_REPLAY_BIT;
        vk_allocateFlags.deviceMask = 0u; // unused

        VkMemoryAllocateInfo vk_allocateInfo = { VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO };
        vk_allocateInfo.pNext = &vk_allocateFlags; // No extensions for now
        vk_allocateInfo.allocationSize = reqs.vulkanReqs.size;
        vk_allocateInfo.memoryTypeIndex = compatibleMemoryTypeIndices[i];

        VkDeviceMemory vk_deviceMemory;
        if (m_devf.vk.vkAllocateMemory(m_vkdev, &vk_allocateInfo, nullptr, &vk_deviceMemory) == VK_SUCCESS)
        {
            // Todo(achal): Change dedicate to not always be false
            return core::make_smart_refctd_ptr<CVulkanMemoryAllocation>(this, reqs.vulkanReqs.size, false, vk_deviceMemory, allocateFlags);
        }
    }

    return nullptr;
}

bool CVulkanLogicalDevice::createCommandBuffers_impl(IGPUCommandPool* cmdPool, IGPUCommandBuffer::E_LEVEL level,
    uint32_t count, core::smart_refctd_ptr<IGPUCommandBuffer>* outCmdBufs)
{
    constexpr uint32_t MAX_COMMAND_BUFFER_COUNT = 1000u;

    if (cmdPool->getAPIType() != EAT_VULKAN)
        return false;

    auto vulkanCommandPool = IBackendObject::device_compatibility_cast<CVulkanCommandPool*>(cmdPool, this)->getInternalObject();

    assert(count <= MAX_COMMAND_BUFFER_COUNT);
    VkCommandBuffer vk_commandBuffers[MAX_COMMAND_BUFFER_COUNT];

    VkCommandBufferAllocateInfo vk_allocateInfo = { VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO };
    vk_allocateInfo.pNext = nullptr; // this must be NULL
    vk_allocateInfo.commandPool = vulkanCommandPool;
    vk_allocateInfo.level = static_cast<VkCommandBufferLevel>(level);
    vk_allocateInfo.commandBufferCount = count;

    if (m_devf.vk.vkAllocateCommandBuffers(m_vkdev, &vk_allocateInfo, vk_commandBuffers) == VK_SUCCESS)
    {
        for (uint32_t i = 0u; i < count; ++i)
        {
            outCmdBufs[i] = core::make_smart_refctd_ptr<CVulkanCommandBuffer>(
                core::smart_refctd_ptr<ILogicalDevice>(this), level, vk_commandBuffers[i],
                core::smart_refctd_ptr<IGPUCommandPool>(cmdPool));
        }

        return true;
    }
    else
    {
        return false;
    }
}

core::smart_refctd_ptr<IGPUImage> CVulkanLogicalDevice::createGPUImage(asset::IImage::SCreationParams&& params)
{
    VkImageCreateInfo vk_createInfo = { VK_STRUCTURE_TYPE_IMAGE_CREATE_INFO };
    vk_createInfo.pNext = nullptr; // there are a lot of extensions
    vk_createInfo.flags = static_cast<VkImageCreateFlags>(params.flags);
    vk_createInfo.imageType = static_cast<VkImageType>(params.type);
    vk_createInfo.format = getVkFormatFromFormat(params.format);
    vk_createInfo.extent = { params.extent.width, params.extent.height, params.extent.depth };
    vk_createInfo.mipLevels = params.mipLevels;
    vk_createInfo.arrayLayers = params.arrayLayers;
    vk_createInfo.samples = static_cast<VkSampleCountFlagBits>(params.samples);
    vk_createInfo.tiling = static_cast<VkImageTiling>(params.tiling);
    vk_createInfo.usage = static_cast<VkImageUsageFlags>(params.usage.value);
    vk_createInfo.sharingMode = static_cast<VkSharingMode>(params.sharingMode);
    vk_createInfo.queueFamilyIndexCount = params.queueFamilyIndexCount;
    vk_createInfo.pQueueFamilyIndices = params.queueFamilyIndices;
    vk_createInfo.initialLayout = static_cast<VkImageLayout>(params.initialLayout);

    VkImage vk_image;
    if (m_devf.vk.vkCreateImage(m_vkdev, &vk_createInfo, nullptr, &vk_image) == VK_SUCCESS)
    {
        VkImageMemoryRequirementsInfo2 vk_memReqsInfo = { VK_STRUCTURE_TYPE_IMAGE_MEMORY_REQUIREMENTS_INFO_2 };
        vk_memReqsInfo.pNext = nullptr;
        vk_memReqsInfo.image = vk_image;

        VkMemoryDedicatedRequirements vk_memDedReqs = { VK_STRUCTURE_TYPE_MEMORY_DEDICATED_REQUIREMENTS };
        VkMemoryRequirements2 vk_memReqs = { VK_STRUCTURE_TYPE_MEMORY_REQUIREMENTS_2 };
        vk_memReqs.pNext = &vk_memDedReqs;

        m_devf.vk.vkGetImageMemoryRequirements2(m_vkdev, &vk_memReqsInfo, &vk_memReqs);

        IDriverMemoryBacked::SDriverMemoryRequirements imageMemReqs = {};
        imageMemReqs.vulkanReqs.alignment = vk_memReqs.memoryRequirements.alignment;
        imageMemReqs.vulkanReqs.size = vk_memReqs.memoryRequirements.size;
        imageMemReqs.vulkanReqs.memoryTypeBits = vk_memReqs.memoryRequirements.memoryTypeBits;
        imageMemReqs.memoryHeapLocation = 0u; // doesn't matter, would get overwritten during memory allocation for this resource anyway
        imageMemReqs.mappingCapability = 0u; // doesn't matter, would get overwritten during memory allocation for this resource anyway
        imageMemReqs.prefersDedicatedAllocation = vk_memDedReqs.prefersDedicatedAllocation;
        imageMemReqs.requiresDedicatedAllocation = vk_memDedReqs.requiresDedicatedAllocation;

        return core::make_smart_refctd_ptr<CVulkanImage>(
            core::smart_refctd_ptr<CVulkanLogicalDevice>(this), std::move(params),
            vk_image, imageMemReqs);
    }
    else
    {
        return nullptr;
    }
}

core::smart_refctd_ptr<IGPUGraphicsPipeline> CVulkanLogicalDevice::createGPUGraphicsPipeline_impl(
    IGPUPipelineCache* pipelineCache,
    IGPUGraphicsPipeline::SCreationParams&& params)
{
    core::smart_refctd_ptr<IGPUGraphicsPipeline> result;
    if (createGPUGraphicsPipelines_impl(pipelineCache, { &params, &params + 1 }, &result))
        return result;
    else
        return nullptr;
}

bool CVulkanLogicalDevice::createGPUGraphicsPipelines_impl(
    IGPUPipelineCache* pipelineCache,
    core::SRange<const IGPUGraphicsPipeline::SCreationParams> params,
    core::smart_refctd_ptr<IGPUGraphicsPipeline>* output)
{
    IGPUGraphicsPipeline::SCreationParams* creationParams = const_cast<IGPUGraphicsPipeline::SCreationParams*>(params.begin());

    VkPipelineCache vk_pipelineCache = VK_NULL_HANDLE;
    if (pipelineCache && pipelineCache->getAPIType() == EAT_VULKAN)
        vk_pipelineCache = IBackendObject::device_compatibility_cast<const CVulkanPipelineCache*>(pipelineCache, this)->getInternalObject();

    // Shader stages
    uint32_t shaderStageCount_total = 0u;
    core::vector<VkPipelineShaderStageCreateInfo> vk_shaderStages(params.size() * IGPURenderpassIndependentPipeline::SHADER_STAGE_COUNT);
    uint32_t specInfoCount_total = 0u;
    core::vector<VkSpecializationInfo> vk_specInfos(vk_shaderStages.size());
    constexpr uint32_t MAX_MAP_ENTRIES_PER_SHADER = 100u;
    uint32_t mapEntryCount_total = 0u;
    core::vector<VkSpecializationMapEntry> vk_mapEntries(vk_specInfos.size()*MAX_MAP_ENTRIES_PER_SHADER);

    // Vertex input
    uint32_t vertexBindingDescriptionCount_total = 0u;
    core::vector<VkVertexInputBindingDescription> vk_vertexBindingDescriptions(params.size() * asset::SVertexInputParams::MAX_ATTR_BUF_BINDING_COUNT);
    uint32_t vertexAttribDescriptionCount_total = 0u;
    core::vector<VkVertexInputAttributeDescription> vk_vertexAttribDescriptions(params.size() * asset::SVertexInputParams::MAX_VERTEX_ATTRIB_COUNT);
    core::vector<VkPipelineVertexInputStateCreateInfo> vk_vertexInputStates(params.size());

    // Input Assembly
    core::vector<VkPipelineInputAssemblyStateCreateInfo> vk_inputAssemblyStates(params.size());

    core::vector<VkPipelineViewportStateCreateInfo> vk_viewportStates(params.size());

    core::vector<VkPipelineRasterizationStateCreateInfo> vk_rasterizationStates(params.size());

    core::vector<VkPipelineMultisampleStateCreateInfo> vk_multisampleStates(params.size());

    core::vector<VkStencilOpState> vk_stencilFrontStates(params.size());
    core::vector<VkStencilOpState> vk_stencilBackStates(params.size());
    core::vector<VkPipelineDepthStencilStateCreateInfo> vk_depthStencilStates(params.size());

    uint32_t colorBlendAttachmentCount_total = 0u;
    core::vector<VkPipelineColorBlendAttachmentState> vk_colorBlendAttachmentStates(params.size() * asset::SBlendParams::MAX_COLOR_ATTACHMENT_COUNT);
    core::vector<VkPipelineColorBlendStateCreateInfo> vk_colorBlendStates(params.size());

    constexpr uint32_t DYNAMIC_STATE_COUNT = 2u;
    VkDynamicState vk_dynamicStates[DYNAMIC_STATE_COUNT] = { VK_DYNAMIC_STATE_VIEWPORT , VK_DYNAMIC_STATE_SCISSOR };
    VkPipelineDynamicStateCreateInfo vk_dynamicStateCreateInfo = { VK_STRUCTURE_TYPE_PIPELINE_DYNAMIC_STATE_CREATE_INFO };
    vk_dynamicStateCreateInfo.pNext = nullptr;
    vk_dynamicStateCreateInfo.flags = 0u;
    vk_dynamicStateCreateInfo.dynamicStateCount = DYNAMIC_STATE_COUNT;
    vk_dynamicStateCreateInfo.pDynamicStates = vk_dynamicStates;

    core::vector<VkGraphicsPipelineCreateInfo> vk_createInfos(params.size());
    for (size_t i = 0ull; i < params.size(); ++i)
    {
        const auto& rpIndie = creationParams[i].renderpassIndependent;

        vk_createInfos[i].sType = VK_STRUCTURE_TYPE_GRAPHICS_PIPELINE_CREATE_INFO;
        vk_createInfos[i].pNext = nullptr;
        vk_createInfos[i].flags = static_cast<VkPipelineCreateFlags>(creationParams[i].createFlags.value);

        uint32_t shaderStageCount = 0u;
        for (uint32_t ss = 0u; ss < IGPURenderpassIndependentPipeline::SHADER_STAGE_COUNT; ++ss)
        {
            const IGPUSpecializedShader* shader = rpIndie->getShaderAtIndex(ss);
            if (!shader || shader->getAPIType() != EAT_VULKAN)
                continue;

            const auto* vulkanSpecShader = IBackendObject::device_compatibility_cast<const CVulkanSpecializedShader*>(shader, this);

            auto& vk_shaderStage = vk_shaderStages[shaderStageCount_total + shaderStageCount];

            vk_shaderStage.sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
            vk_shaderStage.pNext = nullptr;
            vk_shaderStage.flags = 0u;
            vk_shaderStage.stage = static_cast<VkShaderStageFlagBits>(shader->getStage());
            vk_shaderStage.module = vulkanSpecShader->getInternalObject();
            vk_shaderStage.pName = "main";

            const auto& shaderSpecInfo = vulkanSpecShader->getSpecInfo();

            if (shaderSpecInfo.m_backingBuffer && shaderSpecInfo.m_entries)
            {
                for (uint32_t me = 0u; me < shaderSpecInfo.m_entries->size(); ++me)
                {
                    const auto entry = shaderSpecInfo.m_entries->begin() + me;

                    vk_mapEntries[mapEntryCount_total + me].constantID = entry->specConstID;
                    vk_mapEntries[mapEntryCount_total + me].offset = entry->offset;
                    vk_mapEntries[mapEntryCount_total + me].size = entry->size;
                }

                vk_specInfos[specInfoCount_total].mapEntryCount = static_cast<uint32_t>(shaderSpecInfo.m_entries->size());
                vk_specInfos[specInfoCount_total].pMapEntries = vk_mapEntries.data() + mapEntryCount_total;
                mapEntryCount_total += vk_specInfos[specInfoCount_total].mapEntryCount;
                vk_specInfos[specInfoCount_total].dataSize = shaderSpecInfo.m_backingBuffer->getSize();
                vk_specInfos[specInfoCount_total].pData = shaderSpecInfo.m_backingBuffer->getPointer();

                vk_shaderStage.pSpecializationInfo = vk_specInfos.data() + specInfoCount_total++;
            }
            else
            {
                vk_shaderStage.pSpecializationInfo = nullptr;
            }

            ++shaderStageCount;
        }
        vk_createInfos[i].stageCount = shaderStageCount;
        vk_createInfos[i].pStages = vk_shaderStages.data() + shaderStageCount_total;
        shaderStageCount_total += shaderStageCount;

        // Vertex Input
        {
            vk_vertexInputStates[i].sType = VK_STRUCTURE_TYPE_PIPELINE_VERTEX_INPUT_STATE_CREATE_INFO;
            vk_vertexInputStates[i].pNext = nullptr;
            vk_vertexInputStates[i].flags = 0u;

            const auto& vertexInputParams = rpIndie->getVertexInputParams();

            // Fill up vertex binding descriptions
            uint32_t offset = vertexBindingDescriptionCount_total;
            uint32_t vertexBindingDescriptionCount = 0u;

            for (uint32_t b = 0u; b < asset::SVertexInputParams::MAX_ATTR_BUF_BINDING_COUNT; ++b)
            {
                if (vertexInputParams.enabledBindingFlags & (1 << b))
                {
                    auto& bndDesc = vk_vertexBindingDescriptions[offset + vertexBindingDescriptionCount++];

                    bndDesc.binding = b;
                    bndDesc.stride = vertexInputParams.bindings[b].stride;
                    bndDesc.inputRate = static_cast<VkVertexInputRate>(vertexInputParams.bindings[b].inputRate);
                }
            }
            vk_vertexInputStates[i].vertexBindingDescriptionCount = vertexBindingDescriptionCount;
            vk_vertexInputStates[i].pVertexBindingDescriptions = vk_vertexBindingDescriptions.data() + offset;
            vertexBindingDescriptionCount_total += vertexBindingDescriptionCount;

            // Fill up vertex attribute descriptions
            offset = vertexAttribDescriptionCount_total;
            uint32_t vertexAttribDescriptionCount = 0u;

            for (uint32_t l = 0u; l < asset::SVertexInputParams::MAX_VERTEX_ATTRIB_COUNT; ++l)
            {
                if (vertexInputParams.enabledAttribFlags & (1 << l))
                {
                    auto& attribDesc = vk_vertexAttribDescriptions[offset + vertexAttribDescriptionCount++];

                    attribDesc.location = l;
                    attribDesc.binding = vertexInputParams.attributes[l].binding;
                    attribDesc.format = getVkFormatFromFormat(static_cast<asset::E_FORMAT>(vertexInputParams.attributes[l].format));
                    attribDesc.offset = vertexInputParams.attributes[l].relativeOffset;
                }
            }
            vk_vertexInputStates[i].vertexAttributeDescriptionCount = vertexAttribDescriptionCount;
            vk_vertexInputStates[i].pVertexAttributeDescriptions = vk_vertexAttribDescriptions.data() + offset;
            vertexAttribDescriptionCount_total += vertexAttribDescriptionCount;
        }
        vk_createInfos[i].pVertexInputState = &vk_vertexInputStates[i];

        // Input Assembly
        {
            const auto& primAssParams = rpIndie->getPrimitiveAssemblyParams();

            vk_inputAssemblyStates[i].sType = VK_STRUCTURE_TYPE_PIPELINE_INPUT_ASSEMBLY_STATE_CREATE_INFO;
            vk_inputAssemblyStates[i].pNext = nullptr;
            vk_inputAssemblyStates[i].flags = 0u; // reserved for future use by Vulkan
            vk_inputAssemblyStates[i].topology = static_cast<VkPrimitiveTopology>(primAssParams.primitiveType);
            vk_inputAssemblyStates[i].primitiveRestartEnable = primAssParams.primitiveRestartEnable;
        }
        vk_createInfos[i].pInputAssemblyState = &vk_inputAssemblyStates[i];

        // Tesselation
        vk_createInfos[i].pTessellationState = nullptr;

        // Viewport State
        {
            const uint32_t viewportCount = rpIndie->getRasterizationParams().viewportCount;

            vk_viewportStates[i].sType = VK_STRUCTURE_TYPE_PIPELINE_VIEWPORT_STATE_CREATE_INFO;
            vk_viewportStates[i].pNext = nullptr;
            vk_viewportStates[i].flags = 0u;
            vk_viewportStates[i].viewportCount = viewportCount;
            vk_viewportStates[i].pViewports = nullptr; // ignored
            vk_viewportStates[i].scissorCount = viewportCount; // must be identical to viewport count unless VK_DYNAMIC_STATE_VIEWPORT_WITH_COUNT_EXT or VK_DYNAMIC_STATE_SCISSOR_WITH_COUNT_EXT are used
            vk_viewportStates[i].pScissors = nullptr; // ignored
        }
        vk_createInfos[i].pViewportState = &vk_viewportStates[i];

        // Rasterization
        {
            const auto& rasterizationParams = rpIndie->getRasterizationParams();

            vk_rasterizationStates[i].sType = VK_STRUCTURE_TYPE_PIPELINE_RASTERIZATION_STATE_CREATE_INFO;
            vk_rasterizationStates[i].pNext = nullptr;
            vk_rasterizationStates[i].flags = 0u;
            vk_rasterizationStates[i].depthClampEnable = rasterizationParams.depthClampEnable;
            vk_rasterizationStates[i].rasterizerDiscardEnable = rasterizationParams.rasterizerDiscard;
            vk_rasterizationStates[i].polygonMode = static_cast<VkPolygonMode>(rasterizationParams.polygonMode);
            vk_rasterizationStates[i].cullMode = static_cast<VkCullModeFlags>(rasterizationParams.faceCullingMode);
            vk_rasterizationStates[i].frontFace = rasterizationParams.frontFaceIsCCW ? VK_FRONT_FACE_COUNTER_CLOCKWISE : VK_FRONT_FACE_CLOCKWISE;
            vk_rasterizationStates[i].depthBiasEnable = rasterizationParams.depthBiasEnable;
            vk_rasterizationStates[i].depthBiasConstantFactor = rasterizationParams.depthBiasConstantFactor;
            vk_rasterizationStates[i].depthBiasClamp = 0.f;
            vk_rasterizationStates[i].depthBiasSlopeFactor = rasterizationParams.depthBiasSlopeFactor;
            vk_rasterizationStates[i].lineWidth = 1.f;
        }
        vk_createInfos[i].pRasterizationState = &vk_rasterizationStates[i];

        // Multisampling
        {
            const auto& rasterizationParams = rpIndie->getRasterizationParams();

            vk_multisampleStates[i].sType = VK_STRUCTURE_TYPE_PIPELINE_MULTISAMPLE_STATE_CREATE_INFO;
            vk_multisampleStates[i].pNext = nullptr;
            vk_multisampleStates[i].flags = 0u;
            vk_multisampleStates[i].rasterizationSamples = static_cast<VkSampleCountFlagBits>(rasterizationParams.rasterizationSamplesHint);            
            vk_multisampleStates[i].sampleShadingEnable = rasterizationParams.sampleShadingEnable;
            vk_multisampleStates[i].minSampleShading = rasterizationParams.minSampleShading;
            vk_multisampleStates[i].pSampleMask = rasterizationParams.sampleMask;
            vk_multisampleStates[i].alphaToCoverageEnable = rasterizationParams.alphaToCoverageEnable;
            vk_multisampleStates[i].alphaToOneEnable = rasterizationParams.alphaToOneEnable;
        }
        vk_createInfos[i].pMultisampleState = &vk_multisampleStates[i];

        // Depth-stencil
        {
            const auto& rasterParams = rpIndie->getRasterizationParams();

            // Front stencil state
            vk_stencilFrontStates[i].failOp = static_cast<VkStencilOp>(rasterParams.frontStencilOps.failOp);
            vk_stencilFrontStates[i].passOp = static_cast<VkStencilOp>(rasterParams.frontStencilOps.passOp);
            vk_stencilFrontStates[i].depthFailOp = static_cast<VkStencilOp>(rasterParams.frontStencilOps.depthFailOp);
            vk_stencilFrontStates[i].compareOp = static_cast<VkCompareOp>(rasterParams.frontStencilOps.compareOp);
            vk_stencilFrontStates[i].compareMask = 0xFFFFFFFF;
            vk_stencilFrontStates[i].writeMask = rasterParams.frontStencilOps.writeMask;
            vk_stencilFrontStates[i].reference = rasterParams.frontStencilOps.reference;

            // Back stencil state
            vk_stencilBackStates[i].failOp = static_cast<VkStencilOp>(rasterParams.backStencilOps.failOp);
            vk_stencilBackStates[i].passOp = static_cast<VkStencilOp>(rasterParams.backStencilOps.passOp);
            vk_stencilBackStates[i].depthFailOp = static_cast<VkStencilOp>(rasterParams.backStencilOps.depthFailOp);
            vk_stencilBackStates[i].compareOp = static_cast<VkCompareOp>(rasterParams.backStencilOps.compareOp);
            vk_stencilBackStates[i].compareMask = 0xFFFFFFFF;
            vk_stencilBackStates[i].writeMask = rasterParams.backStencilOps.writeMask;
            vk_stencilBackStates[i].reference = rasterParams.backStencilOps.reference;
            
            vk_depthStencilStates[i].sType = VK_STRUCTURE_TYPE_PIPELINE_DEPTH_STENCIL_STATE_CREATE_INFO;
            vk_depthStencilStates[i].pNext = nullptr;
            vk_depthStencilStates[i].flags = static_cast<VkPipelineDepthStencilStateCreateFlags>(0u);
            vk_depthStencilStates[i].depthTestEnable = rasterParams.depthTestEnable;
            vk_depthStencilStates[i].depthWriteEnable = rasterParams.depthWriteEnable;
            vk_depthStencilStates[i].depthCompareOp = static_cast<VkCompareOp>(rasterParams.depthCompareOp);
            vk_depthStencilStates[i].depthBoundsTestEnable = rasterParams.depthBoundsTestEnable;
            vk_depthStencilStates[i].stencilTestEnable = rasterParams.stencilTestEnable;
            vk_depthStencilStates[i].front = vk_stencilFrontStates[i];
            vk_depthStencilStates[i].back = vk_stencilBackStates[i];
            vk_depthStencilStates[i].minDepthBounds = 0.f;
            vk_depthStencilStates[i].maxDepthBounds = 1.f;
        }
        vk_createInfos[i].pDepthStencilState = &vk_depthStencilStates[i];

        // Color blend
        {
            const auto& blendParams = rpIndie->getBlendParams();
            
            uint32_t offset = colorBlendAttachmentCount_total;

            assert(creationParams[i].subpassIx < creationParams[i].renderpass->getCreationParameters().subpassCount);
            auto subpassDescription = creationParams[i].renderpass->getCreationParameters().subpasses[creationParams[i].subpassIx];
            uint32_t colorBlendAttachmentCount = subpassDescription.colorAttachmentCount;

            for (uint32_t as = 0u; as < colorBlendAttachmentCount; ++as)
            {
                const auto& inBlendParams = blendParams.blendParams[as];
                auto& outBlendState = vk_colorBlendAttachmentStates[offset + as];

                outBlendState.blendEnable = inBlendParams.blendEnable;
                outBlendState.srcColorBlendFactor = getVkBlendFactorFromBlendFactor(static_cast<asset::E_BLEND_FACTOR>(inBlendParams.srcColorFactor));
                outBlendState.dstColorBlendFactor = getVkBlendFactorFromBlendFactor(static_cast<asset::E_BLEND_FACTOR>(inBlendParams.dstColorFactor));
                assert(inBlendParams.colorBlendOp <= asset::EBO_MAX);
                outBlendState.colorBlendOp = getVkBlendOpFromBlendOp(static_cast<asset::E_BLEND_OP>(inBlendParams.colorBlendOp));
                outBlendState.srcAlphaBlendFactor = getVkBlendFactorFromBlendFactor(static_cast<asset::E_BLEND_FACTOR>(inBlendParams.srcAlphaFactor));
                outBlendState.dstAlphaBlendFactor = getVkBlendFactorFromBlendFactor(static_cast<asset::E_BLEND_FACTOR>(inBlendParams.dstAlphaFactor));
                assert(inBlendParams.alphaBlendOp <= asset::EBO_MAX);
                outBlendState.alphaBlendOp = getVkBlendOpFromBlendOp(static_cast<asset::E_BLEND_OP>(inBlendParams.alphaBlendOp));
                outBlendState.colorWriteMask = getVkColorComponentFlagsFromColorWriteMask(inBlendParams.colorWriteMask);
            }
            colorBlendAttachmentCount_total += colorBlendAttachmentCount;

            vk_colorBlendStates[i].sType = VK_STRUCTURE_TYPE_PIPELINE_COLOR_BLEND_STATE_CREATE_INFO;
            vk_colorBlendStates[i].pNext = nullptr;
            vk_colorBlendStates[i].flags = 0u;
            vk_colorBlendStates[i].logicOpEnable = blendParams.logicOpEnable;
            vk_colorBlendStates[i].logicOp = getVkLogicOpFromLogicOp(static_cast<asset::E_LOGIC_OP>(blendParams.logicOp));
            vk_colorBlendStates[i].attachmentCount = colorBlendAttachmentCount;
            vk_colorBlendStates[i].pAttachments = vk_colorBlendAttachmentStates.data() + offset;
            vk_colorBlendStates[i].blendConstants[0] = 0.0f;
            vk_colorBlendStates[i].blendConstants[1] = 0.0f;
            vk_colorBlendStates[i].blendConstants[2] = 0.0f;
            vk_colorBlendStates[i].blendConstants[3] = 0.0f;
        }
        vk_createInfos[i].pColorBlendState = &vk_colorBlendStates[i];

        // Dynamic state
        vk_createInfos[i].pDynamicState = &vk_dynamicStateCreateInfo;

        vk_createInfos[i].layout = IBackendObject::device_compatibility_cast<const CVulkanPipelineLayout*>(rpIndie->getLayout(), this)->getInternalObject();
        vk_createInfos[i].renderPass = IBackendObject::device_compatibility_cast<const CVulkanRenderpass*>(creationParams[i].renderpass.get(), this)->getInternalObject();
        vk_createInfos[i].subpass = creationParams[i].subpassIx;
        vk_createInfos[i].basePipelineHandle = VK_NULL_HANDLE;
        vk_createInfos[i].basePipelineIndex = 0u;
    }

    core::vector<VkPipeline> vk_pipelines(params.size());
    if (m_devf.vk.vkCreateGraphicsPipelines(m_vkdev, vk_pipelineCache,
        static_cast<uint32_t>(params.size()), vk_createInfos.data(), nullptr, vk_pipelines.data()) == VK_SUCCESS)
    {
        for (size_t i = 0ull; i < params.size(); ++i)
        {
            output[i] = core::make_smart_refctd_ptr<CVulkanGraphicsPipeline>(
                core::smart_refctd_ptr<CVulkanLogicalDevice>(this),
                std::move(creationParams[i]),
                vk_pipelines[i]);
        }
        return true;
    }
    else
    {
        return false;
    }
}

core::smart_refctd_ptr<IGPUAccelerationStructure> CVulkanLogicalDevice::createGPUAccelerationStructure_impl(IGPUAccelerationStructure::SCreationParams&& params) 
{
    auto physicalDevice = static_cast<const CVulkanPhysicalDevice*>(getPhysicalDevice());
    auto features = physicalDevice->getFeatures();
    
    if(!features.accelerationStructure)
    {
        assert(false && "device accelerationStructures is not enabled.");
        return nullptr;
    }

    VkAccelerationStructureKHR vk_as = VK_NULL_HANDLE;
    VkAccelerationStructureCreateInfoKHR vasci = {VK_STRUCTURE_TYPE_ACCELERATION_STRUCTURE_CREATE_INFO_KHR, nullptr};
    vasci.createFlags = CVulkanAccelerationStructure::getVkASCreateFlagsFromASCreateFlags(params.flags);
    vasci.type = CVulkanAccelerationStructure::getVkASTypeFromASType(params.type);
    vasci.buffer = IBackendObject::device_compatibility_cast<const CVulkanBuffer*>(params.bufferRange.buffer.get(), this)->getInternalObject();
    vasci.offset = static_cast<VkDeviceSize>(params.bufferRange.offset);
    vasci.size = static_cast<VkDeviceSize>(params.bufferRange.size);
    auto vk_res = m_devf.vk.vkCreateAccelerationStructureKHR(m_vkdev, &vasci, nullptr, &vk_as);
    if(VK_SUCCESS != vk_res)
        return nullptr;
    return core::make_smart_refctd_ptr<CVulkanAccelerationStructure>(core::smart_refctd_ptr<CVulkanLogicalDevice>(this), std::move(params), vk_as);
}

bool CVulkanLogicalDevice::buildAccelerationStructures(
    core::smart_refctd_ptr<IDeferredOperation>&& deferredOperation,
    const core::SRange<IGPUAccelerationStructure::HostBuildGeometryInfo>& pInfos,
    IGPUAccelerationStructure::BuildRangeInfo* const* ppBuildRangeInfos)
{
    auto physicalDevice = static_cast<const CVulkanPhysicalDevice*>(getPhysicalDevice());
    auto features = physicalDevice->getFeatures();
    if(!features.accelerationStructure)
    {
        assert(false && "device acceleration structures is not enabled.");
        return false;
    }


    bool ret = false;
    if(!pInfos.empty() && deferredOperation.get() != nullptr)
    {
        VkDeferredOperationKHR vk_deferredOp = IBackendObject::device_compatibility_cast<CVulkanDeferredOperation *>(deferredOperation.get(), this)->getInternalObject();
        static constexpr size_t MaxGeometryPerBuildInfoCount = 64;
        static constexpr size_t MaxBuildInfoCount = 128;
        size_t infoCount = pInfos.size();
        assert(infoCount <= MaxBuildInfoCount);
                
        // TODO: Use better container when ready for these stack allocated memories.
        VkAccelerationStructureBuildGeometryInfoKHR vk_buildGeomsInfos[MaxBuildInfoCount] = {};

        uint32_t geometryArrayOffset = 0u;
        VkAccelerationStructureGeometryKHR vk_geometries[MaxGeometryPerBuildInfoCount * MaxBuildInfoCount] = {};

        IGPUAccelerationStructure::HostBuildGeometryInfo* infos = pInfos.begin();
        for(uint32_t i = 0; i < infoCount; ++i)
        {
            uint32_t geomCount = infos[i].geometries.size();

            assert(geomCount > 0);
            assert(geomCount <= MaxGeometryPerBuildInfoCount);

            vk_buildGeomsInfos[i] = CVulkanAccelerationStructure::getVkASBuildGeomInfoFromBuildGeomInfo(m_vkdev, &m_devf, infos[i], &vk_geometries[geometryArrayOffset]);
            geometryArrayOffset += geomCount; 
        }
                
        static_assert(sizeof(IGPUAccelerationStructure::BuildRangeInfo) == sizeof(VkAccelerationStructureBuildRangeInfoKHR));
        auto buildRangeInfos = reinterpret_cast<const VkAccelerationStructureBuildRangeInfoKHR* const*>(ppBuildRangeInfos);
        VkResult vk_res = m_devf.vk.vkBuildAccelerationStructuresKHR(m_vkdev, vk_deferredOp, infoCount, vk_buildGeomsInfos, buildRangeInfos);
        if(VK_SUCCESS == vk_res)
        {
            ret = true;
        }
    }
    return ret;
}

bool CVulkanLogicalDevice::copyAccelerationStructure(core::smart_refctd_ptr<IDeferredOperation>&& deferredOperation, const IGPUAccelerationStructure::CopyInfo& copyInfo)
{
    auto physicalDevice = static_cast<const CVulkanPhysicalDevice*>(getPhysicalDevice());
    auto features = physicalDevice->getFeatures();
    if(!features.accelerationStructureHostCommands || !features.accelerationStructure)
    {
        assert(false && "device accelerationStructuresHostCommands is not enabled.");
        return false;
    }

    bool ret = false;
    if(deferredOperation.get() != nullptr)
    {
        VkDeferredOperationKHR vk_deferredOp = IBackendObject::device_compatibility_cast<CVulkanDeferredOperation *>(deferredOperation.get(), this)->getInternalObject();
        if(copyInfo.dst == nullptr || copyInfo.src == nullptr) 
        {
            assert(false && "invalid src or dst");
            return false;
        }

        VkCopyAccelerationStructureInfoKHR info = CVulkanAccelerationStructure::getVkASCopyInfo(m_vkdev, &m_devf, copyInfo);
        VkResult res = m_devf.vk.vkCopyAccelerationStructureKHR(m_vkdev, vk_deferredOp, &info);
        if(VK_SUCCESS == res)
        {
            ret = true;
        }
    }
    return ret;
}
    
bool CVulkanLogicalDevice::copyAccelerationStructureToMemory(core::smart_refctd_ptr<IDeferredOperation>&& deferredOperation, const IGPUAccelerationStructure::HostCopyToMemoryInfo& copyInfo)
{
    auto physicalDevice = static_cast<const CVulkanPhysicalDevice*>(getPhysicalDevice());
    auto features = physicalDevice->getFeatures();
    if(!features.accelerationStructureHostCommands || !features.accelerationStructure)
    {
        assert(false && "device accelerationStructuresHostCommands is not enabled.");
        return false;
    }

    bool ret = false;
    if(deferredOperation.get() != nullptr)
    {
        VkDeferredOperationKHR vk_deferredOp = IBackendObject::device_compatibility_cast<CVulkanDeferredOperation *>(deferredOperation.get(), this)->getInternalObject();

        if(copyInfo.dst.isValid() == false || copyInfo.src == nullptr) 
        {
            assert(false && "invalid src or dst");
            return false;
        }

        VkCopyAccelerationStructureToMemoryInfoKHR info = CVulkanAccelerationStructure::getVkASCopyToMemoryInfo(m_vkdev, &m_devf, copyInfo);
        VkResult res = m_devf.vk.vkCopyAccelerationStructureToMemoryKHR(m_vkdev, vk_deferredOp, &info);
        if(VK_SUCCESS == res)
        {
            ret = true;
        }
    }
    return ret;
}

bool CVulkanLogicalDevice::copyAccelerationStructureFromMemory(core::smart_refctd_ptr<IDeferredOperation>&& deferredOperation, const IGPUAccelerationStructure::HostCopyFromMemoryInfo& copyInfo)
{
    auto physicalDevice = static_cast<const CVulkanPhysicalDevice*>(getPhysicalDevice());
    auto features = physicalDevice->getFeatures();
    if(!features.accelerationStructureHostCommands || !features.accelerationStructure)
    {
        assert(false && "device accelerationStructuresHostCommands is not enabled.");
        return false;
    }

    bool ret = false;
    if(deferredOperation.get() != nullptr)
    {
        VkDeferredOperationKHR vk_deferredOp = IBackendObject::device_compatibility_cast<CVulkanDeferredOperation *>(deferredOperation.get(), this)->getInternalObject();
        if(copyInfo.dst == nullptr || copyInfo.src.isValid() == false) 
        {
            assert(false && "invalid src or dst");
            return false;
        }

        VkCopyMemoryToAccelerationStructureInfoKHR info = CVulkanAccelerationStructure::getVkASCopyFromMemoryInfo(m_vkdev, &m_devf, copyInfo);
        VkResult res = m_devf.vk.vkCopyMemoryToAccelerationStructureKHR(m_vkdev, vk_deferredOp, &info);
        if(VK_SUCCESS == res)
        {
            ret = true;
        }
    }
    return ret;
}

IGPUAccelerationStructure::BuildSizes CVulkanLogicalDevice::getAccelerationStructureBuildSizes(const IGPUAccelerationStructure::HostBuildGeometryInfo& pBuildInfo, const uint32_t* pMaxPrimitiveCounts)
{
    // TODO(Validation): Rayquery or RayTracing Pipeline must be enabled
    return getAccelerationStructureBuildSizes_impl(VK_ACCELERATION_STRUCTURE_BUILD_TYPE_HOST_KHR, pBuildInfo, pMaxPrimitiveCounts);
}

IGPUAccelerationStructure::BuildSizes CVulkanLogicalDevice::getAccelerationStructureBuildSizes(const IGPUAccelerationStructure::DeviceBuildGeometryInfo& pBuildInfo, const uint32_t* pMaxPrimitiveCounts)
{
    // TODO(Validation): Rayquery or RayTracing Pipeline must be enabled
    return getAccelerationStructureBuildSizes_impl(VK_ACCELERATION_STRUCTURE_BUILD_TYPE_DEVICE_KHR, pBuildInfo, pMaxPrimitiveCounts);
}

core::smart_refctd_ptr<IQueryPool> CVulkanLogicalDevice::createQueryPool(IQueryPool::SCreationParams&& params)
{
    VkQueryPool vk_queryPool = VK_NULL_HANDLE;
    VkQueryPoolCreateInfo vk_qpci = CVulkanQueryPool::getVkCreateInfoFromCreationParams(std::move(params));
    auto vk_res = m_devf.vk.vkCreateQueryPool(m_vkdev, &vk_qpci, nullptr, &vk_queryPool);
    if(VK_SUCCESS != vk_res)
        return nullptr;
    return core::make_smart_refctd_ptr<CVulkanQueryPool>(core::smart_refctd_ptr<CVulkanLogicalDevice>(this), std::move(params), vk_queryPool);
}

bool CVulkanLogicalDevice::getQueryPoolResults(IQueryPool* queryPool, uint32_t firstQuery, uint32_t queryCount, size_t dataSize, void * pData, uint64_t stride, core::bitflag<IQueryPool::E_QUERY_RESULTS_FLAGS> flags)
{
    bool ret = false;
    if(queryPool != nullptr)
    {
        auto vk_queryPool = IBackendObject::device_compatibility_cast<CVulkanQueryPool*>(queryPool, this)->getInternalObject();
        auto vk_queryResultsflags = CVulkanQueryPool::getVkQueryResultsFlagsFromQueryResultsFlags(flags.value);
        auto vk_res = m_devf.vk.vkGetQueryPoolResults(m_vkdev, vk_queryPool, firstQuery, queryCount, dataSize, pData, static_cast<VkDeviceSize>(stride), vk_queryResultsflags);
        if(VK_SUCCESS == vk_res)
            ret = true;
    }
    return ret;
}

}