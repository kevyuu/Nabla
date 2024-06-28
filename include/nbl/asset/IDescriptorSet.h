// Copyright (C) 2018-2023 - DevSH Graphics Programming Sp. z O.O.
// This file is part of the "Nabla Engine".
// For conditions of distribution and use, see copyright notice in nabla.h
#ifndef _NBL_ASSET_I_DESCRIPTOR_SET_H_INCLUDED_
#define _NBL_ASSET_I_DESCRIPTOR_SET_H_INCLUDED_

#include <algorithm>
#include <compare>


#include "nbl/core/declarations.h"

#include "nbl/asset/format/EFormat.h"
#include "nbl/asset/IDescriptor.h"
#include "nbl/asset/IDescriptorSetLayout.h" //for IDescriptor::E_TYPE
#include "nbl/core/SRange.h"

namespace nbl::asset
{

class IAccelerationStructure;

//! Interface class for various Descriptor Set's resources
/*
	Buffers, Images and Samplers all derive from IDescriptor
	and can be bound under different bindings in a DescriptorSet,
	so that they can be bound to the GPU API together in a single
	API call with efficiency.

	@see ICPUDescriptorSet
	@see IReferenceCounted
*/

template<typename LayoutType>
class IDescriptorSet : public virtual core::IReferenceCounted // TODO: try to remove this inheritance and see what happens
{
	public:
		using layout_t = LayoutType;
		struct SDescriptorInfo
		{
                struct SBufferInfo
                {
                    size_t offset;
                    size_t size;//in Vulkan it's called `range` but IMO it's misleading so i changed to `size`

					auto operator<=>(const SBufferInfo&) const = default;
                };
				struct SImageInfo
				{
					// If the binding is a SAMPLER, this field is ignored
					IImage::LAYOUT imageLayout;
					// If the binding is COMBINED, this will be ignored if the DS layout already has an immutable sampler specified for the binding
					core::smart_refctd_ptr<typename layout_t::sampler_type> sampler;
				};
                    
				core::smart_refctd_ptr<IDescriptor> desc;
				union SBufferImageInfo
				{
					SBufferImageInfo()
					{
						memset(&buffer, 0, core::max<size_t>(sizeof(buffer), sizeof(image)));
					};

					~SBufferImageInfo() {};

					SBufferInfo buffer;
					SImageInfo image;
				} info;

				SDescriptorInfo() {}

				template<typename BufferType>
				SDescriptorInfo(const SBufferBinding<BufferType>& binding) : desc()
				{
					desc = binding.buffer;
					info.buffer.offset = binding.offset;
					info.buffer.size = SBufferRange<BufferType>::WholeBuffer;
				}
				template<typename BufferType>
				SDescriptorInfo(const SBufferRange<BufferType>& range) : desc()
				{
					desc = range.buffer;
					info.buffer.offset = range.offset;
					info.buffer.size = range.size;
				}
				SDescriptorInfo(const SDescriptorInfo& other) : SDescriptorInfo()
				{
					operator=(other);
				}
				SDescriptorInfo(SDescriptorInfo&& other): SDescriptorInfo()
				{
					operator=(std::move(other));
				}
				~SDescriptorInfo()
				{
					if (desc and (desc->getTypeCategory()==IDescriptor::EC_IMAGE or desc->getTypeCategory() == IDescriptor::EC_SAMPLER))
						info.image.sampler = nullptr;
				}

				inline SDescriptorInfo& operator=(const SDescriptorInfo& other)
				{
					if (desc and (desc->getTypeCategory()==IDescriptor::EC_IMAGE or desc->getTypeCategory() == IDescriptor::EC_SAMPLER))
						info.image.sampler = nullptr;
					desc = other.desc;
					const auto type = desc->getTypeCategory();
					if (type!=IDescriptor::EC_IMAGE and type!=IDescriptor::EC_SAMPLER)
						info.buffer = other.info.buffer;
					else
						info.image = other.info.image;
					return *this;
				}
				inline SDescriptorInfo& operator=(SDescriptorInfo&& other)
				{
					if (desc and (desc->getTypeCategory() == IDescriptor::EC_IMAGE or desc->getTypeCategory() == IDescriptor::EC_SAMPLER))
						info.image = {nullptr,IImage::LAYOUT::UNDEFINED};
					desc = std::move(other.desc);
					if (desc)
					{
						const auto type = desc->getTypeCategory();
						if (type!=IDescriptor::EC_IMAGE and type!=IDescriptor::EC_SAMPLER)
							info.buffer = other.info.buffer;
						else
							info.image = other.info.image;
					}
					return *this;
				}

				inline bool operator!=(const SDescriptorInfo& other) const
				{
					if (desc != other.desc)
						return true;
					return info.buffer != other.info.buffer;
				}
		};

		const layout_t* getLayout() const { return m_layout.get(); }

	protected:
		IDescriptorSet(core::smart_refctd_ptr<layout_t>&& _layout) : m_layout(std::move(_layout)) {}
		virtual ~IDescriptorSet() {}

		core::smart_refctd_ptr<layout_t> m_layout;
};

}

#endif
