#ifndef _NBL_BUILTIN_HLSL_CPP_COMPAT_IMPL_INTRINSICS_IMPL_INCLUDED_
#define _NBL_BUILTIN_HLSL_CPP_COMPAT_IMPL_INTRINSICS_IMPL_INCLUDED_

#include <nbl/builtin/hlsl/cpp_compat/basic.h>
#include <nbl/builtin/hlsl/matrix_utils/matrix_traits.hlsl>
#include <nbl/builtin/hlsl/concepts.hlsl>
#include <nbl/builtin/hlsl/spirv_intrinsics/core.hlsl>
#include <nbl/builtin/hlsl/spirv_intrinsics/glsl.std.450.hlsl>
#include <nbl/builtin/hlsl/ieee754.hlsl>
#include <nbl/builtin/hlsl/concepts/core.hlsl>
#include <nbl/builtin/hlsl/concepts/vector.hlsl>
#include <nbl/builtin/hlsl/concepts/matrix.hlsl>
#include <nbl/builtin/hlsl/cpp_compat/promote.hlsl>
#include <nbl/builtin/hlsl/numbers.hlsl>
#ifndef __HLSL_VERSION
#include <boost/preprocessor/comparison/not_equal.hpp>
#include <boost/preprocessor/punctuation/comma_if.hpp>
#include <boost/preprocessor/seq/for_each_i.hpp>
#endif

namespace nbl
{
namespace hlsl
{
namespace cpp_compat_intrinsics_impl
{

template<typename T NBL_STRUCT_CONSTRAINABLE>
struct dot_helper;
template<typename T NBL_STRUCT_CONSTRAINABLE>
struct cross_helper;
template<typename T NBL_STRUCT_CONSTRAINABLE>
struct clamp_helper;
template<typename Integer NBL_STRUCT_CONSTRAINABLE>
struct find_msb_helper;
template<typename Integer NBL_STRUCT_CONSTRAINABLE>
struct find_lsb_helper;
template<typename T NBL_STRUCT_CONSTRAINABLE>
struct bitReverse_helper;
template<typename Matrix NBL_STRUCT_CONSTRAINABLE>
struct transpose_helper;
template<typename Vector NBL_STRUCT_CONSTRAINABLE>
struct length_helper;
template<typename Vector NBL_STRUCT_CONSTRAINABLE>
struct normalize_helper;
template<typename T NBL_STRUCT_CONSTRAINABLE>
struct max_helper;
template<typename T NBL_STRUCT_CONSTRAINABLE>
struct min_helper;
template<typename Integer NBL_STRUCT_CONSTRAINABLE>
struct bitCount_helper;
template<typename LhsT, typename RhsT NBL_STRUCT_CONSTRAINABLE>
struct mul_helper;
template<typename T NBL_STRUCT_CONSTRAINABLE>
struct determinant_helper;
template<typename T NBL_STRUCT_CONSTRAINABLE>
struct inverse_helper;
template<typename T NBL_STRUCT_CONSTRAINABLE>
struct rsqrt_helper;
template<typename T NBL_STRUCT_CONSTRAINABLE>
struct all_helper;
template<typename T NBL_STRUCT_CONSTRAINABLE>
struct any_helper;
template<typename T, uint16_t Bits NBL_STRUCT_CONSTRAINABLE>
struct bitReverseAs_helper;
template<typename T NBL_STRUCT_CONSTRAINABLE>
struct frac_helper;
template<typename T, typename U NBL_STRUCT_CONSTRAINABLE>
struct mix_helper;
template<typename T NBL_STRUCT_CONSTRAINABLE>
struct sign_helper;
template<typename T NBL_STRUCT_CONSTRAINABLE>
struct radians_helper;
template<typename T NBL_STRUCT_CONSTRAINABLE>
struct degrees_helper;
template<typename T NBL_STRUCT_CONSTRAINABLE>
struct step_helper;
template<typename T NBL_STRUCT_CONSTRAINABLE>
struct smoothStep_helper;
template<typename T NBL_STRUCT_CONSTRAINABLE>
struct faceForward_helper;
template<typename T NBL_STRUCT_CONSTRAINABLE>
struct reflect_helper;
template<typename T, typename U NBL_STRUCT_CONSTRAINABLE>
struct refract_helper;

#ifdef __HLSL_VERSION // HLSL only specializations

// it is crucial these partial specializations appear first because thats what makes the helpers match SPIR-V intrinsics first

#define DECLVAL(r,data,i,_T) BOOST_PP_COMMA_IF(BOOST_PP_NOT_EQUAL(i,0)) experimental::declval<_T>()
#define DECL_ARG(r,data,i,_T) BOOST_PP_COMMA_IF(BOOST_PP_NOT_EQUAL(i,0)) const _T arg##i
#define WRAP(r,data,i,_T) BOOST_PP_COMMA_IF(BOOST_PP_NOT_EQUAL(i,0)) _T
#define ARG(r,data,i,_T) BOOST_PP_COMMA_IF(BOOST_PP_NOT_EQUAL(i,0)) arg##i

// the template<> needs to be written ourselves
// return type is __VA_ARGS__ to protect against `,` in templated return types
#define AUTO_SPECIALIZE_TRIVIAL_CASE_HELPER(HELPER_NAME, SPIRV_FUNCTION_NAME, ARG_TYPE_LIST, ARG_TYPE_SET, ...)\
NBL_PARTIAL_REQ_TOP(is_same_v<decltype(spirv::SPIRV_FUNCTION_NAME<T>(BOOST_PP_SEQ_FOR_EACH_I(DECLVAL, _, ARG_TYPE_SET))), __VA_ARGS__ >) \
struct HELPER_NAME<BOOST_PP_SEQ_FOR_EACH_I(WRAP, _, ARG_TYPE_LIST) NBL_PARTIAL_REQ_BOT(is_same_v<decltype(spirv::SPIRV_FUNCTION_NAME<T>(BOOST_PP_SEQ_FOR_EACH_I(DECLVAL, _, ARG_TYPE_SET))), __VA_ARGS__ >) >\
{\
	using return_t = __VA_ARGS__;\
	static inline return_t __call( BOOST_PP_SEQ_FOR_EACH_I(DECL_ARG, _, ARG_TYPE_SET) )\
	{\
		return spirv::SPIRV_FUNCTION_NAME<T>( BOOST_PP_SEQ_FOR_EACH_I(ARG, _, ARG_TYPE_SET) );\
	}\
};

#define FIND_MSB_LSB_RETURN_TYPE conditional_t<is_vector_v<T>, vector<int32_t, vector_traits<T>::Dimension>, int32_t>
template<typename T> AUTO_SPECIALIZE_TRIVIAL_CASE_HELPER(find_msb_helper, findUMsb, (T), (T), FIND_MSB_LSB_RETURN_TYPE);
template<typename T> AUTO_SPECIALIZE_TRIVIAL_CASE_HELPER(find_msb_helper, findSMsb, (T), (T), FIND_MSB_LSB_RETURN_TYPE)
template<typename T> AUTO_SPECIALIZE_TRIVIAL_CASE_HELPER(find_lsb_helper, findILsb, (T), (T), FIND_MSB_LSB_RETURN_TYPE)
#undef FIND_MSB_LSB_RETURN_TYPE

template<typename T> AUTO_SPECIALIZE_TRIVIAL_CASE_HELPER(bitReverse_helper, bitReverse, (T), (T), T)
template<typename T> AUTO_SPECIALIZE_TRIVIAL_CASE_HELPER(transpose_helper, transpose, (T), (T), T)
template<typename T> AUTO_SPECIALIZE_TRIVIAL_CASE_HELPER(length_helper, length, (T), (T), typename vector_traits<T>::scalar_type)
template<typename T> AUTO_SPECIALIZE_TRIVIAL_CASE_HELPER(normalize_helper, normalize, (T), (T), T)
template<typename T> AUTO_SPECIALIZE_TRIVIAL_CASE_HELPER(rsqrt_helper, inverseSqrt, (T), (T), T)
template<typename T> AUTO_SPECIALIZE_TRIVIAL_CASE_HELPER(frac_helper, fract, (T), (T), T)
template<typename T> AUTO_SPECIALIZE_TRIVIAL_CASE_HELPER(all_helper, any, (T), (T), T)
template<typename T> AUTO_SPECIALIZE_TRIVIAL_CASE_HELPER(any_helper, any, (T), (T), T)
template<typename T> AUTO_SPECIALIZE_TRIVIAL_CASE_HELPER(sign_helper, fSign, (T), (T), T)
template<typename T> AUTO_SPECIALIZE_TRIVIAL_CASE_HELPER(sign_helper, sSign, (T), (T), T)
template<typename T> AUTO_SPECIALIZE_TRIVIAL_CASE_HELPER(radians_helper, radians, (T), (T), T)
template<typename T> AUTO_SPECIALIZE_TRIVIAL_CASE_HELPER(degrees_helper, degrees, (T), (T), T)
template<typename T> AUTO_SPECIALIZE_TRIVIAL_CASE_HELPER(max_helper, fMax, (T), (T)(T), T)
template<typename T> AUTO_SPECIALIZE_TRIVIAL_CASE_HELPER(max_helper, uMax, (T), (T)(T), T)
template<typename T> AUTO_SPECIALIZE_TRIVIAL_CASE_HELPER(max_helper, sMax, (T), (T)(T), T)
template<typename T> AUTO_SPECIALIZE_TRIVIAL_CASE_HELPER(min_helper, fMin, (T), (T)(T), T)
template<typename T> AUTO_SPECIALIZE_TRIVIAL_CASE_HELPER(min_helper, uMin, (T), (T)(T), T)
template<typename T> AUTO_SPECIALIZE_TRIVIAL_CASE_HELPER(min_helper, sMin, (T), (T)(T), T)
template<typename T> AUTO_SPECIALIZE_TRIVIAL_CASE_HELPER(step_helper, step, (T), (T)(T), T)
template<typename T> AUTO_SPECIALIZE_TRIVIAL_CASE_HELPER(reflect_helper, reflect, (T), (T)(T), T)
template<typename T> AUTO_SPECIALIZE_TRIVIAL_CASE_HELPER(clamp_helper, fClamp, (T), (T)(T)(T), T)
template<typename T> AUTO_SPECIALIZE_TRIVIAL_CASE_HELPER(clamp_helper, uClamp, (T), (T)(T)(T), T)
template<typename T> AUTO_SPECIALIZE_TRIVIAL_CASE_HELPER(clamp_helper, sClamp, (T), (T)(T)(T), T)
template<typename T> AUTO_SPECIALIZE_TRIVIAL_CASE_HELPER(smoothStep_helper, smoothStep, (T), (T)(T)(T), T)
template<typename T> AUTO_SPECIALIZE_TRIVIAL_CASE_HELPER(faceForward_helper, faceForward, (T), (T)(T)(T), T)
template<typename T, typename U> AUTO_SPECIALIZE_TRIVIAL_CASE_HELPER(refract_helper, refract, (T)(U), (T)(T)(U), T)

#define BITCOUNT_HELPER_RETRUN_TYPE conditional_t<is_vector_v<T>, vector<int32_t, vector_traits<T>::Dimension>, int32_t>
template<typename T> AUTO_SPECIALIZE_TRIVIAL_CASE_HELPER(bitCount_helper, bitCount, (T), (T), BITCOUNT_HELPER_RETRUN_TYPE)
#undef BITCOUNT_HELPER_RETRUN_TYPE

#undef DECLVAL
#undef DECL_ARG
#undef WRAP
#undef ARG
#undef AUTO_SPECIALIZE_TRIVIAL_CASE_HELPER

template<typename UInt64> NBL_PARTIAL_REQ_TOP(is_same_v<UInt64, uint64_t>)
struct find_msb_helper<UInt64 NBL_PARTIAL_REQ_BOT(is_same_v<UInt64, uint64_t>) >
{
	using return_t = int32_t;
	static return_t __call(NBL_CONST_REF_ARG(UInt64) val)
	{
		const uint32_t highBits = uint32_t(val >> 32);
		const int32_t highMsb = find_msb_helper<uint32_t>::__call(highBits);

		if (highMsb == -1)
		{
			const uint32_t lowBits = uint32_t(val);
			const int32_t lowMsb = find_msb_helper<uint32_t>::__call(lowBits);
			if (lowMsb == -1)
				return -1;

			return lowMsb;
		}

		return highMsb + 32;
	}
};
template<typename UInt64> NBL_PARTIAL_REQ_TOP(is_same_v<UInt64, uint64_t>)
struct find_lsb_helper<UInt64 NBL_PARTIAL_REQ_BOT(is_same_v<UInt64, uint64_t>) >
{
	static int32_t __call(NBL_CONST_REF_ARG(uint64_t) val)
	{
		const uint32_t lowBits = uint32_t(val);
		const int32_t lowLsb = find_lsb_helper<uint32_t>::__call(lowBits);

		if (lowLsb == -1)
		{
			const uint32_t highBits = uint32_t(val >> 32);
			const int32_t highLsb = find_lsb_helper<uint32_t>::__call(highBits);
			if (highLsb == -1)
				return -1;
			else
				return 32 + highLsb;
		}

		return lowLsb;
	}
};

template<typename SquareMatrix>
NBL_PARTIAL_REQ_TOP(concepts::Matrix<SquareMatrix>&& matrix_traits<SquareMatrix>::Square)
struct inverse_helper<SquareMatrix NBL_PARTIAL_REQ_BOT(concepts::Matrix<SquareMatrix>&& matrix_traits<SquareMatrix>::Square) >
{
	static SquareMatrix __call(NBL_CONST_REF_ARG(SquareMatrix) mat)
	{
		return spirv::matrixInverse(mat);
	}
};

template<typename T, typename U> NBL_PARTIAL_REQ_TOP(always_true<decltype(spirv::fMix<T>(experimental::declval<T>(), experimental::declval<T>(), experimental::declval<U>()))>)
struct mix_helper<T, U NBL_PARTIAL_REQ_BOT(always_true<decltype(spirv::fMix<T>(experimental::declval<T>(), experimental::declval<T>(), experimental::declval<U>()))>) >
{
	using return_t = conditional_t<is_vector_v<T>, vector<typename vector_traits<T>::scalar_type, vector_traits<T>::Dimension>, T>;
	static inline return_t __call(const T x, const T y, const U a)
	{
		T aAsT = a;
		return spirv::fMix<T>(x, y, aAsT);
	}
};

template<typename SquareMatrix> NBL_PARTIAL_REQ_TOP(matrix_traits<SquareMatrix>::Square)
struct determinant_helper<SquareMatrix NBL_PARTIAL_REQ_BOT(matrix_traits<SquareMatrix>::Square) >
{
	static typename matrix_traits<SquareMatrix>::scalar_type __call(NBL_CONST_REF_ARG(SquareMatrix) mat)
	{
		return spirv::determinant(mat);
	}
};

#else // C++ only specializations

#define DECL_ARG(r,data,i,_T) BOOST_PP_COMMA_IF(BOOST_PP_NOT_EQUAL(i,0)) const _T arg##i
#define WRAP(r,data,i,_T) BOOST_PP_COMMA_IF(BOOST_PP_NOT_EQUAL(i,0)) _T
#define ARG(r,data,i,_T) BOOST_PP_COMMA_IF(BOOST_PP_NOT_EQUAL(i,0)) arg##i

// the template<> needs to be written ourselves
// return type is __VA_ARGS__ to protect against `,` in templated return types
#define AUTO_SPECIALIZE_TRIVIAL_CASE_HELPER(HELPER_NAME, STD_FUNCTION_NAME, REQUIREMENT, ARG_TYPE_LIST, ARG_TYPE_SET, ...)\
requires REQUIREMENT \
struct HELPER_NAME<BOOST_PP_SEQ_FOR_EACH_I(WRAP, _, ARG_TYPE_LIST)>\
{\
	using return_t = __VA_ARGS__;\
	static inline return_t __call( BOOST_PP_SEQ_FOR_EACH_I(DECL_ARG, _, ARG_TYPE_SET) )\
	{\
		return std::STD_FUNCTION_NAME<BOOST_PP_SEQ_FOR_EACH_I(WRAP, _, ARG_TYPE_LIST)>( BOOST_PP_SEQ_FOR_EACH_I(ARG, _, ARG_TYPE_SET) );\
	}\
};

template<typename T> AUTO_SPECIALIZE_TRIVIAL_CASE_HELPER(clamp_helper, clamp, concepts::Scalar<T>, (T), (T)(T)(T), T)
template<typename T> AUTO_SPECIALIZE_TRIVIAL_CASE_HELPER(max_helper, max, concepts::Scalar<T>, (T), (T)(T), T)
template<typename T> AUTO_SPECIALIZE_TRIVIAL_CASE_HELPER(min_helper, min, concepts::Scalar<T>, (T), (T)(T), T)

#undef DECL_ARG
#undef WRAP
#undef ARG
#undef AUTO_SPECIALIZE_TRIVIAL_CASE_HELPER

template<typename T>
requires concepts::IntegralScalar<T>
struct bitReverse_helper<T>
{
	static inline T __call(NBL_CONST_REF_ARG(T) arg)
	{
		return glm::bitfieldReverse<T>(arg);
	}
};
template<typename Matrix>
requires concepts::Matrix<Matrix>
struct transpose_helper<Matrix>
{
	using transposed_t = typename matrix_traits<Matrix>::transposed_type;

	static transposed_t __call(NBL_CONST_REF_ARG(Matrix) m)
	{
		return reinterpret_cast<transposed_t&>(glm::transpose(reinterpret_cast<typename Matrix::Base const&>(m)));
	}
};
template<typename Vector>
requires concepts::FloatingPointVector<Vector>
struct length_helper<Vector>
{
	static inline typename vector_traits<Vector>::scalar_type __call(NBL_CONST_REF_ARG(Vector) vec)
	{
		return std::sqrt(dot_helper<Vector>::__call(vec, vec));
	}
};
template<typename Vectorial>
requires concepts::FloatingPointLikeVectorial<Vectorial>
struct normalize_helper<Vectorial>
{
	static inline Vectorial __call(NBL_CONST_REF_ARG(Vectorial) vec)
	{
		return vec / length_helper<Vectorial>::__call(vec);
	}
};

template<typename T>
requires concepts::IntegralScalar<T>
struct find_lsb_helper<T>
{
	using return_t = int32_t;
	static inline T __call(const T arg)
	{
		return glm::findLSB<T>(arg);
	}
};
template<typename Integer>
NBL_PARTIAL_REQ_TOP(concepts::IntegralScalar<Integer>)
struct find_msb_helper<Integer NBL_PARTIAL_REQ_BOT(concepts::IntegralScalar<Integer>) >
{
	using return_t = int32_t;
	static return_t __call(NBL_CONST_REF_ARG(Integer) val)
	{
		return glm::findMSB<Integer>(val);
	}
};
// TODO: implemet to be compatible with both C++ and HLSL when it works with DXC
template<typename EnumType>
requires std::is_enum_v<EnumType>
struct find_lsb_helper<EnumType>
{
	using return_t = int32_t;
	static int32_t __call(NBL_CONST_REF_ARG(EnumType) val)
	{
		using underlying_t = std::underlying_type_t<EnumType>;
		return find_lsb_helper<underlying_t>::__call(static_cast<underlying_t>(val));
	}
};
template<typename EnumType>
requires std::is_enum_v<EnumType>
struct find_msb_helper<EnumType>
{
	using return_t = int32_t;
	static return_t __call(NBL_CONST_REF_ARG(EnumType) val)
	{
		using underlying_t = std::underlying_type_t<EnumType>;
		return find_msb_helper<underlying_t>::__call(static_cast<underlying_t>(val));
	}
};

template<typename FloatingPoint>
requires concepts::FloatingPointScalar<FloatingPoint>
struct rsqrt_helper<FloatingPoint>
{
	static FloatingPoint __call(NBL_CONST_REF_ARG(FloatingPoint) x)
	{
		// TODO: https://stackoverflow.com/a/62239778
		return 1.0f / std::sqrt(x);
	}
};

template<typename T>
requires concepts::FloatingPointScalar<T>
struct frac_helper<T>
{
	using return_t = T;
	static inline return_t __call(const T x)
	{
		return x - std::floor(x);
	}
};

template<typename Integer>
requires concepts::IntegralScalar<Integer>
struct bitCount_helper<Integer>
{
	using return_t = int32_t;
	static return_t __call(NBL_CONST_REF_ARG(Integer) val)
	{
		using UnsignedInteger = typename hlsl::unsigned_integer_of_size_t<sizeof(Integer)>;
		return std::popcount(static_cast<UnsignedInteger>(val));
	}
};

template<typename SquareMatrix>
requires concepts::Matrix<SquareMatrix> && matrix_traits<SquareMatrix>::Square
struct inverse_helper<SquareMatrix>
{
	static SquareMatrix __call(NBL_CONST_REF_ARG(SquareMatrix) mat)
	{
		return reinterpret_cast<SquareMatrix&>(glm::inverse(reinterpret_cast<typename SquareMatrix::Base const&>(mat)));
	}
};

template<typename EnumT>
requires std::is_enum_v<EnumT>
struct bitCount_helper<EnumT>
{
	using return_t = int32_t;
	using underlying_t = std::underlying_type_t<EnumT>;
	static return_t __call(NBL_CONST_REF_ARG(EnumT) val)
	{
		return bitCount_helper<const underlying_t>::__call(reinterpret_cast<const underlying_t&>(val));
	}
};

template<typename T, typename U>
requires concepts::FloatingPoint<T> && (concepts::FloatingPoint<T> || concepts::Boolean<T>)
struct mix_helper<T, U>
{
	using return_t = T;
	static inline return_t __call(const T x, const T y, const U a)
	{
		return glm::mix(x, y, a);
	}
};

template<typename T>
requires concepts::FloatingPointScalar<T> || concepts::IntegralScalar<T>
struct sign_helper<T>
{
	using return_t = T;
	static inline return_t __call(const T val)
	{
		if (val < 0)
			return -1;
		if (val > 0)
			return 1;

		return 0;
	}
};

template<typename T>
requires concepts::FloatingPointScalar<T>
struct radians_helper<T>
{
	using return_t = T;
	static inline return_t __call(const T degrees)
	{
		return degrees * (numbers::pi<T> / static_cast<T>(180.0));
	}
};

template<typename T>
requires concepts::FloatingPointScalar<T>
struct degrees_helper<T>
{
	using return_t = T;
	static inline return_t __call(const T radians)
	{
		return radians * (static_cast<T>(180.0) / numbers::pi<T>);
	}
};

template<typename T>
requires concepts::FloatingPointScalar<T>
struct step_helper<T>
{
	using return_t = T;
	static inline return_t __call(const T edge, const T x)
	{
		return x < edge ? 0.0 : 1.0;
	}
};

template<typename T>
requires concepts::FloatingPointScalar<T>
struct smoothStep_helper<T>
{
	using return_t = T;
	static inline return_t __call(const T edge0, const T edge1, const T x)
	{
		T t = clamp_helper<T>::__call((x - edge0) / (edge1 - edge0), 0, 1);
		return t * t * (3 - 2 * t);
	}
};

template<typename SquareMatrix>
NBL_PARTIAL_REQ_TOP(matrix_traits<SquareMatrix>::Square)
struct determinant_helper<SquareMatrix NBL_PARTIAL_REQ_BOT(matrix_traits<SquareMatrix>::Square) >
{
	static typename matrix_traits<SquareMatrix>::scalar_type __call(NBL_CONST_REF_ARG(SquareMatrix) mat)
	{
		return glm::determinant(reinterpret_cast<typename SquareMatrix::Base const&>(mat));
	}
};

template<typename T>
requires concepts::FloatingPointVectorial<T>
struct faceForward_helper<T>
{
	using return_t = T;
	static inline return_t __call(const T N, const T I, const T Nref)
	{
		if (dot_helper<T>::__call(Nref, I) < 0.0)
			return N;
		else
			return -N;
	}
};

template<typename T>
requires concepts::FloatingPointVector<T>
struct reflect_helper<T>
{
	using return_t = T;
	static inline return_t __call(const T I, const T N)
	{
		return I - T(2.0 * dot_helper<T>::__call(N, I)) * N;
	}
};

template<typename T, typename U>
requires concepts::FloatingPointVector<T> && concepts::FloatingPointScalar<U>
struct refract_helper<T, U>
{
	using return_t = T;
	static inline return_t __call(const T I, const T N, const U eta)
	{
		U k = 1.0 - eta * eta * (1.0 - dot_helper<T>::__call(N, I) * dot_helper<T>::__call(N, I));
		if (k < 0.0)
			return T(0.0);
		
		return eta * I - (eta * dot_helper<T>::__call(N, I) + std::sqrt(k)) * N;
	}
};

#endif // C++ only specializations

// C++ and HLSL specializations

template<typename T, uint16_t Bits>
NBL_PARTIAL_REQ_TOP(concepts::UnsignedIntegralScalar<T> && (Bits <= sizeof(T) * 8))
struct bitReverseAs_helper<T, Bits NBL_PARTIAL_REQ_BOT(concepts::UnsignedIntegralScalar<T> && (Bits <= sizeof(T) * 8)) >
{
	static T __call(NBL_CONST_REF_ARG(T) val)
	{
		return bitReverse_helper<T>::__call(val) >> promote<T, scalar_type_t<T> >(scalar_type_t <T>(sizeof(T) * 8 - Bits));
	}

	static T __call(NBL_CONST_REF_ARG(T) val, uint16_t bits)
	{
		return bitReverse_helper<T>::__call(val) >> promote<T, scalar_type_t<T> >(scalar_type_t <T>(sizeof(T) * 8 - bits));
	}
};

template<typename Vectorial>
NBL_PARTIAL_REQ_TOP(concepts::Vectorial<Vectorial>)
struct dot_helper<Vectorial NBL_PARTIAL_REQ_BOT(concepts::Vectorial<Vectorial>) >
{
	using scalar_type = typename vector_traits<Vectorial>::scalar_type;

	static inline scalar_type __call(NBL_CONST_REF_ARG(Vectorial) lhs, NBL_CONST_REF_ARG(Vectorial) rhs)
	{
		static const uint32_t ArrayDim = vector_traits<Vectorial>::Dimension;
		static array_get<Vectorial, scalar_type> getter;

		scalar_type retval = getter(lhs, 0) * getter(rhs, 0);
		for (uint32_t i = 1; i < ArrayDim; ++i)
			retval = retval + getter(lhs, i) * getter(rhs, i);

		return retval;
	}
};
template<typename FloatingPointLikeVectorial>
NBL_PARTIAL_REQ_TOP(concepts::FloatingPointLikeVectorial<FloatingPointLikeVectorial> && (vector_traits<FloatingPointLikeVectorial>::Dimension == 3))
struct cross_helper<FloatingPointLikeVectorial NBL_PARTIAL_REQ_BOT(concepts::FloatingPointLikeVectorial<FloatingPointLikeVectorial> && (vector_traits<FloatingPointLikeVectorial>::Dimension == 3)) >
{
	static FloatingPointLikeVectorial __call(NBL_CONST_REF_ARG(FloatingPointLikeVectorial) lhs, NBL_CONST_REF_ARG(FloatingPointLikeVectorial) rhs)
	{
#ifdef __HLSL_VERSION
		return spirv::cross(lhs, rhs);
#else
		using traits = hlsl::vector_traits<FloatingPointLikeVectorial>;
		array_get<FloatingPointLikeVectorial, typename traits::scalar_type> getter;
		array_set<FloatingPointLikeVectorial, typename traits::scalar_type> setter;

		FloatingPointLikeVectorial output;
		setter(output, 0, getter(lhs, 1) * getter(rhs, 2) - getter(rhs, 1) * getter(lhs, 2));
		setter(output, 1, getter(lhs, 2) * getter(rhs, 0) - getter(rhs, 2) * getter(lhs, 0));
		setter(output, 2, getter(lhs, 0) * getter(rhs, 1) - getter(rhs, 0) * getter(lhs, 1));

		return output;
#endif
	}
};

#ifdef __HLSL_VERSION
// SPIR-V already defines specializations for builtin vector types
#define VECTOR_SPECIALIZATION_CONCEPT concepts::Vectorial<T> && !is_vector_v<T>
#else
#define VECTOR_SPECIALIZATION_CONCEPT concepts::Vectorial<T>
#endif

template<typename T>
NBL_PARTIAL_REQ_TOP(VECTOR_SPECIALIZATION_CONCEPT)
struct clamp_helper<T NBL_PARTIAL_REQ_BOT(VECTOR_SPECIALIZATION_CONCEPT) >
{
	using return_t = T;
	static return_t __call(NBL_CONST_REF_ARG(T) val, NBL_CONST_REF_ARG(T) min, NBL_CONST_REF_ARG(T) max)
	{
		using traits = hlsl::vector_traits<T>;
		array_get<T, typename traits::scalar_type> getter;
		array_set<return_t, typename traits::scalar_type> setter;

		return_t output;
		for (uint32_t i = 0; i < traits::Dimension; ++i)
			setter(output, i, clamp_helper<typename traits::scalar_type>::__call(getter(val, i), getter(min, i), getter(max, i)));

		return output;
	}
};

template<typename T>
NBL_PARTIAL_REQ_TOP(VECTOR_SPECIALIZATION_CONCEPT)
struct min_helper<T NBL_PARTIAL_REQ_BOT(VECTOR_SPECIALIZATION_CONCEPT) >
{
	static T __call(NBL_CONST_REF_ARG(T) a, NBL_CONST_REF_ARG(T) b)
	{
		using traits = hlsl::vector_traits<T>;
		array_get<T, typename traits::scalar_type> getter;
		array_set<T, typename traits::scalar_type> setter;

		T output;
		for (uint32_t i = 0; i < traits::Dimension; ++i)
			setter(output, i, min_helper<typename traits::scalar_type>::__call(getter(a, i), getter(b, i)));

		return output;
	}
};
template<typename T>
NBL_PARTIAL_REQ_TOP(VECTOR_SPECIALIZATION_CONCEPT)
struct max_helper<T NBL_PARTIAL_REQ_BOT(VECTOR_SPECIALIZATION_CONCEPT) >
{
	static T __call(NBL_CONST_REF_ARG(T) a, NBL_CONST_REF_ARG(T) b)
	{
		using traits = hlsl::vector_traits<T>;
		array_get<T, typename traits::scalar_type> getter;
		array_set<T, typename traits::scalar_type> setter;

		T output;
		for (uint32_t i = 0; i < traits::Dimension; ++i)
			setter(output, i, max_helper<typename traits::scalar_type>::__call(getter(a, i), getter(b, i)));

		return output;
	}
};

template<typename LhsT, typename RhsT>
NBL_PARTIAL_REQ_TOP(concepts::Matrix<LhsT> && concepts::Vector<RhsT> && (matrix_traits<LhsT>::ColumnCount == vector_traits<RhsT>::Dimension))
struct mul_helper<LhsT, RhsT NBL_PARTIAL_REQ_BOT(concepts::Matricial<LhsT> && concepts::Vectorial<RhsT> && (matrix_traits<LhsT>::column_count == vector_traits<RhsT>::Dimension)) >
{
	using lhs_traits = matrix_traits<LhsT>;
	using rhs_traits = vector_traits<RhsT>;
	using return_t = vector<typename lhs_traits::scalar_type, lhs_traits::RowCount>;
	static inline return_t __call(LhsT lhs, RhsT rhs)
	{
		return mul(lhs, rhs);
	}
};

template<typename LhsT, typename RhsT>
NBL_PARTIAL_REQ_TOP(concepts::Matrix<LhsT> && concepts::Matrix<RhsT> && (matrix_traits<LhsT>::ColumnCount == matrix_traits<RhsT>::RowCount))
struct mul_helper<LhsT, RhsT NBL_PARTIAL_REQ_BOT(concepts::Matrix<LhsT> && concepts::Matrix<RhsT> && (matrix_traits<LhsT>::ColumnCount == matrix_traits<RhsT>::RowCount)) >
{
	using lhs_traits = matrix_traits<LhsT>;
	using rhs_traits = matrix_traits<RhsT>;
	using return_t = matrix<typename lhs_traits::scalar_type, lhs_traits::RowCount, rhs_traits::ColumnCount>;
	static inline return_t __call(LhsT lhs, RhsT rhs)
	{
		return mul(lhs, rhs);
	}
};

#define AUTO_SPECIALIZE_HELPER_FOR_VECTOR(HELPER_NAME, REQUIREMENT, RETURN_TYPE)\
template<typename T>\
NBL_PARTIAL_REQ_TOP(REQUIREMENT)\
struct HELPER_NAME<T NBL_PARTIAL_REQ_BOT(REQUIREMENT) >\
{\
	using return_t = RETURN_TYPE;\
	static return_t __call(NBL_CONST_REF_ARG(T) vec)\
	{\
		using traits = hlsl::vector_traits<T>;\
		using return_t_traits = hlsl::vector_traits<return_t>;\
		array_get<T, typename traits::scalar_type> getter;\
		array_set<return_t, typename return_t_traits::scalar_type> setter;\
\
		return_t output;\
		for (uint32_t i = 0; i < traits::Dimension; ++i)\
			setter(output, i, HELPER_NAME<typename traits::scalar_type>::__call(getter(vec, i)));\
\
		return output;\
	}\
};

AUTO_SPECIALIZE_HELPER_FOR_VECTOR(rsqrt_helper, concepts::FloatingPointVectorial<T> && VECTOR_SPECIALIZATION_CONCEPT, T)
AUTO_SPECIALIZE_HELPER_FOR_VECTOR(bitReverse_helper, VECTOR_SPECIALIZATION_CONCEPT, T)
AUTO_SPECIALIZE_HELPER_FOR_VECTOR(frac_helper, VECTOR_SPECIALIZATION_CONCEPT,T)
AUTO_SPECIALIZE_HELPER_FOR_VECTOR(sign_helper, VECTOR_SPECIALIZATION_CONCEPT, T)
AUTO_SPECIALIZE_HELPER_FOR_VECTOR(degrees_helper, VECTOR_SPECIALIZATION_CONCEPT, T)
AUTO_SPECIALIZE_HELPER_FOR_VECTOR(radians_helper, VECTOR_SPECIALIZATION_CONCEPT, T)
#define INT32_VECTOR_TYPE vector<int32_t, hlsl::vector_traits<T>::Dimension>
AUTO_SPECIALIZE_HELPER_FOR_VECTOR(bitCount_helper, VECTOR_SPECIALIZATION_CONCEPT, INT32_VECTOR_TYPE)
AUTO_SPECIALIZE_HELPER_FOR_VECTOR(find_msb_helper, VECTOR_SPECIALIZATION_CONCEPT, INT32_VECTOR_TYPE)
AUTO_SPECIALIZE_HELPER_FOR_VECTOR(find_lsb_helper, VECTOR_SPECIALIZATION_CONCEPT, INT32_VECTOR_TYPE)
#undef INT32_VECTOR_TYPE
#undef AUTO_SPECIALIZE_HELPER_FOR_VECTOR

template<typename BooleanVector>
NBL_PARTIAL_REQ_TOP(concepts::Vectorial<BooleanVector> && is_same_v<typename vector_traits<BooleanVector>::scalar_type, bool>)
struct all_helper<BooleanVector NBL_PARTIAL_REQ_BOT(concepts::Vectorial<BooleanVector> && is_same_v<typename vector_traits<BooleanVector>::scalar_type, bool>) >
{
	static bool __call(NBL_CONST_REF_ARG(BooleanVector) x)
	{
		using traits = hlsl::vector_traits<BooleanVector>;
		array_get<BooleanVector, typename traits::scalar_type> getter;
		array_set<BooleanVector, typename traits::scalar_type> setter;

		bool output = true;
		for (uint32_t i = 0; i < traits::Dimension; ++i)
			output = output && getter(x, i);

		return output;
	}
};

template<typename BooleanVector>
NBL_PARTIAL_REQ_TOP(concepts::Vectorial<BooleanVector> && is_same_v<typename vector_traits<BooleanVector>::scalar_type, bool>)
struct any_helper<BooleanVector NBL_PARTIAL_REQ_BOT(concepts::Vectorial<BooleanVector> && is_same_v<typename vector_traits<BooleanVector>::scalar_type, bool>) >
{
	static bool __call(NBL_CONST_REF_ARG(BooleanVector) x)
	{
		using traits = hlsl::vector_traits<BooleanVector>;
		array_get<BooleanVector, typename traits::scalar_type> getter;
		array_set<BooleanVector, typename traits::scalar_type> setter;

		bool output = false;
		for (uint32_t i = 0; i < traits::Dimension; ++i)
			output = output || getter(x, i);

		return output;
	}
};

template<typename T>
NBL_PARTIAL_REQ_TOP(VECTOR_SPECIALIZATION_CONCEPT)
struct step_helper<T NBL_PARTIAL_REQ_BOT(VECTOR_SPECIALIZATION_CONCEPT) >
{
	using return_t = T;
	static return_t __call(NBL_CONST_REF_ARG(T) edge, NBL_CONST_REF_ARG(T) x)
	{
		using traits = hlsl::vector_traits<T>;
		array_get<T, typename traits::scalar_type> getter;
		array_set<return_t, typename traits::scalar_type> setter;
		
		return_t output;
		for (uint32_t i = 0; i < traits::Dimension; ++i)
			setter(output, i, step_helper<typename traits::scalar_type>::__call(getter(edge, i), getter(x, i)));
		
		return output;
	}
};

template<typename T>
NBL_PARTIAL_REQ_TOP(VECTOR_SPECIALIZATION_CONCEPT)
struct smoothStep_helper<T NBL_PARTIAL_REQ_BOT(VECTOR_SPECIALIZATION_CONCEPT) >
{
	using return_t = T;
	static return_t __call(NBL_CONST_REF_ARG(T) edge0, NBL_CONST_REF_ARG(T) edge1, NBL_CONST_REF_ARG(T) x)
	{
		using traits = hlsl::vector_traits<T>;
		array_get<T, typename traits::scalar_type> getter;
		array_set<return_t, typename traits::scalar_type> setter;

		return_t output;
		for (uint32_t i = 0; i < traits::Dimension; ++i)
			setter(output, i, smoothStep_helper<typename traits::scalar_type>::__call(getter(edge0, i), getter(edge1, i), getter(x, i)));

		return output;
	}
};

}
}
}

#endif