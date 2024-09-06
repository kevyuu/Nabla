#ifndef _NBL_BUILTIN_HLSL_EMULATED_FLOAT64_T_INCLUDED_
#define _NBL_BUILTIN_HLSL_EMULATED_FLOAT64_T_INCLUDED_

#include <nbl/builtin/hlsl/emulated/emulated_float64_t_impl.hlsl>

namespace nbl
{
namespace hlsl
{
    template<bool FastMath = true, bool FlushDenormToZero = true>
    struct emulated_float64_t
    {
        using storage_t = uint64_t;
        using this_t = emulated_float64_t<FastMath, FlushDenormToZero>;

        storage_t data;

        // constructors
        /*static emulated_float64_t create(uint16_t val)
        {
            return emulated_float64_t(bit_cast<uint64_t>(float64_t(val)));
        }*/

        NBL_CONSTEXPR_STATIC_INLINE this_t create(this_t val)
        {
            return val;
        }

        NBL_CONSTEXPR_STATIC_INLINE this_t create(int32_t val)
        {
            return bit_cast<this_t>(impl::castToUint64WithFloat64BitPattern(int64_t(val)));
        }

        NBL_CONSTEXPR_STATIC_INLINE this_t create(int64_t val)
        {
            return bit_cast<this_t>(impl::castToUint64WithFloat64BitPattern(val));
        }

        NBL_CONSTEXPR_STATIC_INLINE this_t create(uint32_t val)
        {
            return bit_cast<this_t>(impl::castToUint64WithFloat64BitPattern(uint64_t(val)));
        }

        NBL_CONSTEXPR_STATIC_INLINE this_t create(uint64_t val)
        {
            return bit_cast<this_t>(impl::castToUint64WithFloat64BitPattern(val));
        }

        NBL_CONSTEXPR_STATIC_INLINE this_t create(float32_t val)
        {
            this_t output;
            output.data = impl::castFloat32ToStorageType<FlushDenormToZero>(val);
            return output;
        }

        NBL_CONSTEXPR_STATIC_INLINE this_t create(float64_t val)
        {
#ifdef __HLSL_VERSION
            emulated_float64_t retval;
            uint32_t lo, hi;
            asuint(val, lo, hi);
            retval.data = (uint64_t(hi) << 32) | lo;
            return retval;
#else
            return bit_cast<this_t>(reinterpret_cast<uint64_t&>(val));
#endif
        }
        
        // TODO: unresolved external symbol imath_half_to_float_table
        /*static emulated_float64_t create(float16_t val)
        {
            return emulated_float64_t(bit_cast<uint64_t>(float64_t(val)));
        }*/

        // arithmetic operators
        this_t operator+(const emulated_float64_t rhs) NBL_CONST_MEMBER_FUNC
        {
            if (FlushDenormToZero)
            {
                if (!FastMath && (tgmath::isnan(data) || tgmath::isnan(rhs.data)))
                    return bit_cast<this_t>(ieee754::traits<float64_t>::quietNaN);

                if (!FastMath && impl::areBothInfinity(data, rhs.data))
                {
                    uint64_t lhsSign = data & ieee754::traits<float64_t>::signMask;
                    uint64_t rhsSign = rhs.data & ieee754::traits<float64_t>::signMask;

                    if(lhsSign == rhsSign)
                        return bit_cast<this_t>(ieee754::traits<float64_t>::inf | lhsSign);
                    else if(lhsSign || rhsSign)
                        return bit_cast<this_t>(ieee754::traits<float64_t>::quietNaN | ieee754::traits<float64_t>::signMask);
                }

                if (!FastMath && tgmath::isinf(data))
                    return bit_cast<this_t>(data);

                if (!FastMath && tgmath::isinf(rhs.data))
                    return bit_cast<this_t>(rhs.data);

                const int lhsBiasedExp = ieee754::extractBiasedExponent(data);
                const int rhsBiasedExp = ieee754::extractBiasedExponent(rhs.data);

                uint64_t lhsData = impl::flushDenormToZero(lhsBiasedExp, data);
                uint64_t rhsData = impl::flushDenormToZero(rhsBiasedExp, rhs.data);

                uint64_t lhsSign = ieee754::extractSignPreserveBitPattern(lhsData);
                uint64_t rhsSign = ieee754::extractSignPreserveBitPattern(rhsData);

                if (!FastMath && impl::areBothZero(lhsData, rhsData))
                {
                    if (lhsSign == rhsSign)
                        return bit_cast<this_t>(lhsSign);
                    else
                        return bit_cast<this_t>(0ull);
                }

                if (!FastMath && tgmath::isinf(lhsData))
                    return bit_cast<this_t>(ieee754::traits<float64_t>::inf | ieee754::extractSignPreserveBitPattern(max(lhsData, rhsData)));

                uint64_t lhsNormMantissa = ieee754::extractNormalizeMantissa(lhsData);
                uint64_t rhsNormMantissa = ieee754::extractNormalizeMantissa(rhsData);

                const int expDiff = lhsBiasedExp - rhsBiasedExp;

                const int exp = max(lhsBiasedExp, rhsBiasedExp) - ieee754::traits<float64_t>::exponentBias;
                const uint32_t shiftAmount = abs(expDiff);

                if (expDiff < 0)
                {
                    // so lhsNormMantissa always holds mantissa of number with greater exponent
                    swap<uint64_t>(lhsNormMantissa, rhsNormMantissa);
                    swap<uint64_t>(lhsSign, rhsSign);
                }

                rhsNormMantissa >>= shiftAmount;

                uint64_t resultMantissa;
                if (lhsSign != rhsSign)
                {
                    int64_t mantissaDiff = lhsNormMantissa - rhsNormMantissa;
                    if (mantissaDiff < 0)
                    {
                        swap<uint64_t>(lhsNormMantissa, rhsNormMantissa);
                        swap<uint64_t>(lhsSign, rhsSign);
                    }

                    lhsNormMantissa <<= 10;
                    rhsNormMantissa <<= 10;
                    resultMantissa = uint64_t(int64_t(lhsNormMantissa) - int64_t(rhsNormMantissa));
                    resultMantissa >>= 10;
                }
                else
                {
                    resultMantissa = lhsNormMantissa + rhsNormMantissa;
                }

                uint64_t resultBiasedExp = uint64_t(exp) + ieee754::traits<float64_t>::exponentBias;

                if (resultMantissa & 1ull << 53)
                {
                    ++resultBiasedExp;
                    resultMantissa >>= 1;
                }

                while (resultMantissa < (1ull << 52))
                {
                    --resultBiasedExp;
                    resultMantissa <<= 1;
                }

                resultMantissa &= ieee754::traits<float64_t>::mantissaMask;
                uint64_t output = impl::assembleFloat64(lhsSign, resultBiasedExp << ieee754::traits<float64_t>::mantissaBitCnt, resultMantissa);
                return bit_cast<this_t>(output);
            }

            // not implemented
            if (!FlushDenormToZero)
                return bit_cast<this_t>(0xdeadbeefbadcaffeull);
        }

        emulated_float64_t operator+(float rhs)
        {
            return bit_cast<this_t>(data) + create(rhs);
        }

        emulated_float64_t operator-(emulated_float64_t rhs) NBL_CONST_MEMBER_FUNC
        {
            emulated_float64_t lhs = bit_cast<this_t>(data);
            emulated_float64_t rhsFlipped = rhs.flipSign();
            
            return lhs + rhsFlipped;
        }

        emulated_float64_t operator-(float rhs) NBL_CONST_MEMBER_FUNC
        {
            return bit_cast<this_t>(data) - create(rhs);
        }

        emulated_float64_t operator*(emulated_float64_t rhs) NBL_CONST_MEMBER_FUNC
        {
            if(FlushDenormToZero)
            {
                emulated_float64_t retval = this_t::create(0ull);

                int lhsBiasedExp = ieee754::extractBiasedExponent(data);
                int rhsBiasedExp = ieee754::extractBiasedExponent(rhs.data);

                uint64_t lhsData = impl::flushDenormToZero(lhsBiasedExp, data);
                uint64_t rhsData = impl::flushDenormToZero(rhsBiasedExp, rhs.data);

                uint64_t lhsSign = lhsData & ieee754::traits<float64_t>::signMask;
                uint64_t rhsSign = rhsData & ieee754::traits<float64_t>::signMask;

                uint64_t lhsMantissa = ieee754::extractMantissa(lhsData);
                uint64_t rhsMantissa = ieee754::extractMantissa(rhsData);

                int exp = int(lhsBiasedExp + rhsBiasedExp) - ieee754::traits<float64_t>::exponentBias;
                uint64_t sign = (lhsData ^ rhsData) & ieee754::traits<float64_t>::signMask;

                if (!FastMath && (tgmath::isnan(lhsData) || tgmath::isnan(rhsData)))
                    return bit_cast<this_t>(ieee754::traits<float64_t>::quietNaN | sign);
                if (!FastMath && (tgmath::isinf(lhsData) || tgmath::isinf(rhsData)))
                    return bit_cast<this_t>(ieee754::traits<float64_t>::inf | sign);
                if (!FastMath && impl::areBothZero(lhsData, rhsData))
                    return bit_cast<this_t>(sign);

                const uint64_t hi_l = (lhsMantissa >> 21) | (1ull << 31);
                const uint64_t lo_l = lhsMantissa & ((1ull << 21) - 1);
                const uint64_t hi_r = (rhsMantissa >> 21) | (1ull << 31);
                const uint64_t lo_r = rhsMantissa & ((1ull << 21) - 1);

                //const uint64_t RoundToNearest = (1ull << 31) - 1;
                uint64_t newPseudoMantissa = ((hi_l * hi_r) >> 10) + ((hi_l * lo_r + lo_l * hi_r/* + RoundToNearest*/) >> 31);

                if (newPseudoMantissa & (0x1ull << 53))
                {
                    newPseudoMantissa >>= 1;
                    ++exp;
                }
                newPseudoMantissa &= (ieee754::traits<float64_t>::mantissaMask);

                return bit_cast<this_t>(impl::assembleFloat64(sign, uint64_t(exp) << ieee754::traits<float64_t>::mantissaBitCnt, newPseudoMantissa));
            }
            else
            {
                //static_assert(false, "not implemented yet");
                return bit_cast<this_t>(0xdeadbeefbadcaffeull);
            }
        }

        emulated_float64_t operator*(float rhs)
        {
            return _static_cast<this_t>(data) * create(rhs);
        }

        emulated_float64_t operator/(const emulated_float64_t rhs) NBL_CONST_MEMBER_FUNC
        {
            if (FlushDenormToZero)
            {
                if (!FastMath && (tgmath::isnan<uint64_t>(data) || tgmath::isnan<uint64_t>(rhs.data)))
                    return bit_cast<this_t>(ieee754::traits<float64_t>::quietNaN);

                const uint64_t sign = (data ^ rhs.data) & ieee754::traits<float64_t>::signMask;

                if (!FastMath && impl::isZero(rhs.data))
                    return bit_cast<this_t>(ieee754::traits<float64_t>::quietNaN | sign);
                if (!FastMath && impl::areBothInfinity(data, rhs.data))
                    return bit_cast<this_t>(ieee754::traits<float64_t>::quietNaN | ieee754::traits<float64_t>::signMask);
                if (!FastMath && tgmath::isinf(data))
                    return bit_cast<this_t>(ieee754::traits<float64_t>::inf | sign);
                if (!FastMath && tgmath::isinf(rhs.data))
                    return bit_cast<this_t>(0ull | sign);
                if (!FastMath && impl::isZero(rhs.data))
                    return bit_cast<this_t>(ieee754::traits<float64_t>::quietNaN | sign);

                int lhsBiasedExp = ieee754::extractBiasedExponent(data);
                int rhsBiasedExp = ieee754::extractBiasedExponent(rhs.data);

                uint64_t lhsData = impl::flushDenormToZero(lhsBiasedExp, data);
                uint64_t rhsData = impl::flushDenormToZero(rhsBiasedExp, rhs.data);

                const uint64_t lhsRealMantissa = (ieee754::extractMantissa(lhsData) | (1ull << ieee754::traits<float64_t>::mantissaBitCnt));
                const uint64_t rhsRealMantissa = ieee754::extractMantissa(rhsData) | (1ull << ieee754::traits<float64_t>::mantissaBitCnt);

                int exp = lhsBiasedExp - rhsBiasedExp + int(ieee754::traits<float64_t>::exponentBias);

                uint64_t2 lhsMantissaShifted = impl::shiftMantissaLeftBy53(lhsRealMantissa);
                uint64_t mantissa = impl::divmod128by64(lhsMantissaShifted.x, lhsMantissaShifted.y, rhsRealMantissa);

                while (mantissa < (1ull << 52))
                {
                    mantissa <<= 1;
                    exp--;
                }

                mantissa &= ieee754::traits<float64_t>::mantissaMask;

                return bit_cast<this_t>(impl::assembleFloat64(sign, uint64_t(exp) << ieee754::traits<float64_t>::mantissaBitCnt, mantissa));
            }
            else
            {
                //static_assert(false, "not implemented yet");
                return bit_cast<this_t>(0xdeadbeefbadcaffeull);
            }
        }

        // relational operators
        // TODO: should `FlushDenormToZero` affect relational operators?
        bool operator==(this_t rhs) NBL_CONST_MEMBER_FUNC
        {
            if (!FastMath && (tgmath::isnan<uint64_t>(data) || tgmath::isnan<uint64_t>(rhs.data)))
                return false;
            if (!FastMath && impl::areBothZero(data, rhs.data))
                return true;

            const emulated_float64_t xored = bit_cast<this_t>(data ^ rhs.data);
            if ((xored.data & 0x7FFFFFFFFFFFFFFFull) == 0ull)
                return true;

            return !(xored.data);
        }
        bool operator!=(emulated_float64_t rhs) NBL_CONST_MEMBER_FUNC
        {
            if (!FastMath && (tgmath::isnan<uint64_t>(data) || tgmath::isnan<uint64_t>(rhs.data)))
                return false;

            return !(bit_cast<this_t>(data) == rhs);
        }
        bool operator<(emulated_float64_t rhs) NBL_CONST_MEMBER_FUNC
        {
            if (!FastMath && (tgmath::isnan<uint64_t>(data) || tgmath::isnan<uint64_t>(rhs.data)))
                return false;
            if (!FastMath && impl::areBothInfinity(data, rhs.data))
                return false;
            if (!FastMath && impl::areBothZero(data, rhs.data))
                return false;

            const uint64_t lhsSign = ieee754::extractSign(data);
            const uint64_t rhsSign = ieee754::extractSign(rhs.data);

            // flip bits of negative numbers and flip signs of all numbers
            uint64_t lhsFlipped = data ^ ((0x7FFFFFFFFFFFFFFFull * lhsSign) | ieee754::traits<float64_t>::signMask);
            uint64_t rhsFlipped = rhs.data ^ ((0x7FFFFFFFFFFFFFFFull * rhsSign) | ieee754::traits<float64_t>::signMask);

            return lhsFlipped < rhsFlipped;
        }
        bool operator>(emulated_float64_t rhs) NBL_CONST_MEMBER_FUNC
        {
            if (!FastMath && (tgmath::isnan<uint64_t>(data) || tgmath::isnan<uint64_t>(rhs.data)))
                return false;
            if (!FastMath && impl::areBothInfinity(data, rhs.data))
                return false;
            if (!FastMath && impl::areBothZero(data, rhs.data))
                return false;

            const uint64_t lhsSign = ieee754::extractSign(data);
            const uint64_t rhsSign = ieee754::extractSign(rhs.data);

            // flip bits of negative numbers and flip signs of all numbers
            uint64_t lhsFlipped = data ^ ((0x7FFFFFFFFFFFFFFFull * lhsSign) | ieee754::traits<float64_t>::signMask);
            uint64_t rhsFlipped = rhs.data ^ ((0x7FFFFFFFFFFFFFFFull * rhsSign) | ieee754::traits<float64_t>::signMask);

            return lhsFlipped > rhsFlipped;
        }
        bool operator<=(emulated_float64_t rhs) NBL_CONST_MEMBER_FUNC 
        { 
            if (!FastMath && (tgmath::isnan<uint64_t>(data) || tgmath::isnan<uint64_t>(rhs.data)))
                return false;

            return !(bit_cast<this_t>(data) > bit_cast<this_t>(rhs.data));
        }
        bool operator>=(emulated_float64_t rhs)
        {
            if (!FastMath && (tgmath::isnan<uint64_t>(data) || tgmath::isnan<uint64_t>(rhs.data)))
                return false;

            return !(bit_cast<this_t>(data) < bit_cast<this_t>(rhs.data));
        }

        //logical operators
        bool operator&&(emulated_float64_t rhs) NBL_CONST_MEMBER_FUNC { return bool(data) && bool(rhs.data); }
        bool operator||(emulated_float64_t rhs) NBL_CONST_MEMBER_FUNC { return bool(data) || bool(rhs.data); }
        bool operator!() NBL_CONST_MEMBER_FUNC { return !bool(data); }

        emulated_float64_t flipSign()
        {
            return bit_cast<this_t>(data ^ ieee754::traits<float64_t>::signMask);
        }

        NBL_CONSTEXPR_STATIC_INLINE bool supportsFastMath()
        {
            return FastMath;
        }

        enum E_ROUNDING_MODE
        {
            FLOAT_ROUND_NEAREST_EVEN,
            FLOAT_ROUND_TO_ZERO,
            FLOAT_ROUND_DOWN,
            FLOAT_ROUND_UP
        };

        static const E_ROUNDING_MODE RoundingMode = E_ROUNDING_MODE::FLOAT_ROUND_TO_ZERO;
    };

#define IMPLEMENT_IEEE754_FUNC_SPEC_FOR_EMULATED_F64_TYPE(...) \
template<>\
struct traits_base<__VA_ARGS__ >\
{\
    NBL_CONSTEXPR_STATIC_INLINE int16_t exponentBitCnt = 11;\
    NBL_CONSTEXPR_STATIC_INLINE int16_t mantissaBitCnt = 52;\
};\
template<>\
inline uint32_t extractBiasedExponent(__VA_ARGS__ x)\
{\
    return extractBiasedExponent<uint64_t>(x.data);\
}\
\
template<>\
inline int extractExponent(__VA_ARGS__ x)\
{\
    return extractExponent(x.data);\
}\
\
template<>\
NBL_CONSTEXPR_INLINE_FUNC __VA_ARGS__ replaceBiasedExponent(__VA_ARGS__ x, typename unsigned_integer_of_size<sizeof(__VA_ARGS__)>::type biasedExp)\
{\
    return __VA_ARGS__(replaceBiasedExponent(x.data, biasedExp));\
}\
\
template <>\
NBL_CONSTEXPR_INLINE_FUNC __VA_ARGS__ fastMulExp2(__VA_ARGS__ x, int n)\
{\
    return __VA_ARGS__(replaceBiasedExponent(x.data, extractBiasedExponent(x) + uint32_t(n)));\
}\
\
template <>\
NBL_CONSTEXPR_INLINE_FUNC unsigned_integer_of_size<sizeof(__VA_ARGS__)>::type extractMantissa(__VA_ARGS__ x)\
{\
    return extractMantissa(x.data);\
}\
\
template <>\
NBL_CONSTEXPR_INLINE_FUNC uint64_t extractNormalizeMantissa(__VA_ARGS__ x)\
{\
    return extractNormalizeMantissa(x.data);\
}\
\

#define DEFINE_BIT_CAST_SPEC(...)\
template<>\
NBL_CONSTEXPR_FUNC __VA_ARGS__ bit_cast<__VA_ARGS__, uint64_t>(NBL_CONST_REF_ARG(uint64_t) val)\
{\
__VA_ARGS__ output;\
output.data = val;\
\
return output;\
}\
\
template<>\
NBL_CONSTEXPR_FUNC __VA_ARGS__ bit_cast<__VA_ARGS__, float64_t>(NBL_CONST_REF_ARG(float64_t) val)\
{\
__VA_ARGS__ output;\
output.data = bit_cast<uint64_t>(val);\
\
return output;\
}\
\

namespace impl
{

template<typename To, bool FastMath, bool FlushDenormToZero>
struct static_cast_helper<To,emulated_float64_t<FastMath,FlushDenormToZero>,void>
{
    static_assert(is_scalar<To>::value);

    using From = emulated_float64_t<FastMath,FlushDenormToZero>;

    static inline To cast(From v)
    {
        using ToAsFloat = typename float_of_size<sizeof(To)>::type;
        using ToAsUint = typename unsigned_integer_of_size<sizeof(To)>::type;


        if (is_same_v<To, float64_t>)
            return To(bit_cast<float64_t>(v.data));

        if (is_floating_point<To>::value)
        {

            const int exponent = ieee754::extractExponent(v.data);
            if (!From::supportsFastMath())
            {
                if (exponent > ieee754::traits<ToAsFloat>::exponentMax)
                    return bit_cast<To>(ieee754::traits<ToAsFloat>::inf);
                if (exponent < ieee754::traits<ToAsFloat>::exponentMin)
                    return bit_cast<To>(-ieee754::traits<ToAsFloat>::inf);
                if (tgmath::isnan(v.data))
                    return bit_cast<To>(ieee754::traits<ToAsFloat>::quietNaN);
            }


            const uint32_t toBitSize = sizeof(To) * 8;
            const ToAsUint sign = ToAsUint(ieee754::extractSign(v.data) << (toBitSize - 1));
            const ToAsUint biasedExponent = ToAsUint(exponent + ieee754::traits<ToAsFloat>::exponentBias) << ieee754::traits<ToAsFloat>::mantissaBitCnt;
            const ToAsUint mantissa = ToAsUint(v.data >> (ieee754::traits<float64_t>::mantissaBitCnt - ieee754::traits<ToAsFloat>::mantissaBitCnt)) & ieee754::traits<ToAsFloat>::mantissaMask;

            return bit_cast<ToAsFloat>(sign | biasedExponent | mantissa);
        }

        // NOTE: casting from negative float to unsigned int is an UB, function will return abs value in this case
        if (is_integral<To>::value)
        {
            const int exponent = ieee754::extractExponent(v.data);
            if (exponent < 0)
                return 0;

            uint64_t unsignedOutput = ieee754::extractMantissa(v.data) & 1ull << ieee754::traits<float64_t>::mantissaBitCnt;
            const int shiftAmount = exponent - int(ieee754::traits<float64_t>::mantissaBitCnt);

            if (shiftAmount < 0)
                unsignedOutput <<= -shiftAmount;
            else
                unsignedOutput >>= shiftAmount;

            if (is_signed<To>::value)
            {
                int64_t signedOutput64 = unsignedOutput & ((1ull << 63) - 1);
                To signedOutput = To(signedOutput64);
                if (ieee754::extractSignPreserveBitPattern(v.data) != 0)
                    signedOutput = -signedOutput;

                return signedOutput;
            }

            return To(unsignedOutput);
        }

        // assert(false);
        return To(0xdeadbeefbadcaffeull);
    }
};

template<bool FastMath, bool FlushDenormToZero>
struct static_cast_helper<emulated_float64_t<FastMath, FlushDenormToZero>, float32_t, void>
{
    using To = emulated_float64_t<FastMath, FlushDenormToZero>;

    static inline To cast(float32_t v)
    {
        return To::create(v);
    }
};

template<bool FastMath, bool FlushDenormToZero>
struct static_cast_helper<emulated_float64_t<FastMath, FlushDenormToZero>, float64_t, void>
{
    using To = emulated_float64_t<FastMath, FlushDenormToZero>;

    static inline To cast(float64_t v)
    {
        return To::create(v);
    }
};

template<bool FastMath, bool FlushDenormToZero>
struct static_cast_helper<emulated_float64_t<FastMath, FlushDenormToZero>, uint32_t, void>
{
    using To = emulated_float64_t<FastMath, FlushDenormToZero>;

    static inline To cast(uint32_t v)
    {
        return To::create(v);
    }
};

template<bool FastMath, bool FlushDenormToZero>
struct static_cast_helper<emulated_float64_t<FastMath, FlushDenormToZero>, uint64_t, void>
{
    using To = emulated_float64_t<FastMath, FlushDenormToZero>;

    static inline To cast(uint64_t v)
    {
        return To::create(v);
    }
};

template<bool FastMath, bool FlushDenormToZero>
struct static_cast_helper<emulated_float64_t<FastMath, FlushDenormToZero>, emulated_float64_t<FastMath, FlushDenormToZero>, void>
{
    static inline emulated_float64_t<FastMath, FlushDenormToZero> cast(emulated_float64_t<FastMath, FlushDenormToZero> v)
    {
        return v;
    }
};

}

DEFINE_BIT_CAST_SPEC(emulated_float64_t<true, true>);
DEFINE_BIT_CAST_SPEC(emulated_float64_t<false, false>);
DEFINE_BIT_CAST_SPEC(emulated_float64_t<true, false>);
DEFINE_BIT_CAST_SPEC(emulated_float64_t<false, true>);

//template<bool FastMath, bool FlushDenormToZero>
//struct is_floating_point<emulated_float64_t<FastMath, FlushDenormToZero> > : bool_constant<true> {};

namespace ieee754
{
IMPLEMENT_IEEE754_FUNC_SPEC_FOR_EMULATED_F64_TYPE(emulated_float64_t<true, true>);
IMPLEMENT_IEEE754_FUNC_SPEC_FOR_EMULATED_F64_TYPE(emulated_float64_t<false, false>);
IMPLEMENT_IEEE754_FUNC_SPEC_FOR_EMULATED_F64_TYPE(emulated_float64_t<true, false>);
IMPLEMENT_IEEE754_FUNC_SPEC_FOR_EMULATED_F64_TYPE(emulated_float64_t<false, true>);
}

}
}

#undef IMPLEMENT_IEEE754_FUNC_SPEC_FOR_EMULATED_F64_TYPE
#undef DEFINE_BIT_CAST_SPEC

#endif