// Copyright (C) 2018-2020 - DevSH Graphics Programming Sp. z O.O.
// This file is part of the "Nabla Engine".
// For conditions of distribution and use, see copyright notice in nabla.h
#ifndef _NBL_EXT_MITSUBA_LOADER_C_ELEMENT_TRANSFORM_H_INCLUDED_
#define _NBL_EXT_MITSUBA_LOADER_C_ELEMENT_TRANSFORM_H_INCLUDED_

#include "nbl/ext/MitsubaLoader/IElement.h"

namespace nbl::ext::MitsubaLoader
{

class CElementTransform : public IElement
{
	public:
		CElementTransform() : IElement(""), matrix() {}
		virtual ~CElementTransform() {}

		bool addProperty(SNamedPropertyElement&& _property) override;
		bool onEndTag(asset::IAssetLoader::IAssetLoaderOverride* _override, CMitsubaMetadata* globalMetadata) override { return true; }
		IElement::Type getType() const override { return IElement::Type::TRANSFORM; }
		std::string getLogName() const override { return "transform"; }
		/*
		inline CElementTransform& operator=(const CElementTransform& other)
		{
			IElement::operator=(other);
			matrix = other.matrix;
			return *this;
		}
		*/

		core::matrix4SIMD matrix;
};

}

#endif