// Copyright (C) 2018-2020 - DevSH Graphics Programming Sp. z O.O.
// This file is part of the "Nabla Engine".
// For conditions of distribution and use, see copyright notice in nabla.h
#ifndef _NBL_EXT_MITSUBA_LOADER_C_MITSUBA_LOADER_H_INCLUDED_
#define _NBL_EXT_MITSUBA_LOADER_C_MITSUBA_LOADER_H_INCLUDED_

#include "nbl/asset/asset.h"

#include "IFileSystem.h"
#include "nbl/asset/utils/ICPUVirtualTexture.h"

#include "nbl/ext/MitsubaLoader/CSerializedLoader.h"
#include "nbl/ext/MitsubaLoader/CMitsubaMetadata.h"
#include "nbl/ext/MitsubaLoader/CElementShape.h"
#include "nbl/ext/MitsubaLoader/SContext.h"


namespace nbl::ext::MitsubaLoader
{

namespace impl
{
#define uint uint32_t
#define uvec2 uint64_t
#define mat4x3 nbl::core::matrix3x4SIMD
#define nbl_glsl_MC_material_data_t asset::material_compiler::impl::nbl_glsl_MC_material_data_t
struct vec3
{
	float x, y, z;
};
#include <nbl/builtin/glsl/ext/MitsubaLoader/instance_data_struct.glsl>
#undef uint
#undef uvec2
#undef mat4x3
#undef nbl_glsl_MC_material_data_t
}


class CElementBSDF;
class CMaterialCompilerFrontend;

class CMitsubaLoader : public asset::IRenderpassIndependentPipelineLoader
{
		friend class CMitsubaMaterialCompilerFrontend;

	public:
		using instance_data_t = impl::nbl_glsl_ext_Mitsuba_Loader_instance_data_t;
		//! Constructor
		CMitsubaLoader(asset::IAssetManager* _manager, io::IFileSystem* _fs);

		void initialize() override;

		//! Check if the file might be loaded by this class
		/** Check might look into the file.
		\param file File handle to check.
		\return True if file seems to be loadable. */
		bool isALoadableFileFormat(io::IReadFile* _file) const override;

		//! Returns an array of string literals terminated by nullptr
		const char** getAssociatedFileExtensions() const override;

		//! Returns the assets loaded by the loader
		/** Bits of the returned value correspond to each IAsset::E_TYPE
		enumeration member, and the return value cannot be 0. */
		uint64_t getSupportedAssetTypesBitfield() const override { return asset::IAsset::ET_MESH/*|asset::IAsset::ET_SCENE|asset::IAsset::ET_IMPLEMENTATION_SPECIFIC_METADATA*/; }

		//! Loads an asset from an opened file, returns nullptr in case of failure.
		asset::SAssetBundle loadAsset(io::IReadFile* _file, const asset::IAssetLoader::SAssetLoadParams& _params, asset::IAssetLoader::IAssetLoaderOverride* _override = nullptr, uint32_t _hierarchyLevel = 0u) override;
		
	protected:
		//! Destructor
		virtual ~CMitsubaLoader() = default;

		static core::smart_refctd_ptr<asset::ICPUPipelineLayout> createPipelineLayout(asset::IAssetManager* _manager, asset::ICPUVirtualTexture* _vt);

		//
		core::vector<SContext::shape_ass_type>	getMesh(SContext& ctx, uint32_t hierarchyLevel, CElementShape* shape);
		core::vector<SContext::shape_ass_type>	loadShapeGroup(SContext& ctx, uint32_t hierarchyLevel, const CElementShape::ShapeGroup* shapegroup, const core::matrix3x4SIMD& relTform);
		SContext::shape_ass_type				loadBasicShape(SContext& ctx, uint32_t hierarchyLevel, CElementShape* shape, const core::matrix3x4SIMD& relTform);
		
		void									cacheTexture(SContext& ctx, uint32_t hierarchyLevel, const CElementTexture* texture, const CMaterialCompilerFrontend::E_IMAGE_VIEW_SEMANTIC semantic);

		SContext::bsdf_type getBSDFtreeTraversal(SContext& ctx, const CElementBSDF* bsdf);
		SContext::bsdf_type genBSDFtreeTraversal(SContext& ctx, const CElementBSDF* bsdf);

		template <typename Iter>
		core::smart_refctd_ptr<asset::ICPUDescriptorSet> createDS0(const SContext& _ctx, asset::ICPUPipelineLayout* _layout, const asset::material_compiler::CGLSLBackendCommon::result_t& _compResult, Iter meshBegin, Iter meshEnd);


		// members
		io::IFileSystem* m_filesystem;
};

}
#endif