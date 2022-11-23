// Copyright (C) 2022 - DevSH Graphics Programming Sp. z O.O.
// This file is part of the "Nabla Engine".
// For conditions of distribution and use, see copyright notice in nabla.h
#ifndef _NBL_BUILTIN_HLSL_BINOPS_INCLUDED_
#define _NBL_BUILTIN_HLSL_BINOPS_INCLUDED_

namespace nbl
{
namespace hlsl
{
namespace binops
{
template<typename T>
struct bitwise_and
{
    T operator()(const T lhs, const T rhs)
    {
        return lhs&rhs;
    }
	
	T identity()
	{
		return ~0;
	}
};

template<typename T>
struct bitwise_or
{
    T operator()(const T lhs, const T rhs)
    {
        return lhs|rhs;
    }
	
	T identity()
	{
		return 0;
	}
};

template<typename T>
struct bitwise_xor
{
    T operator()(const T lhs, const T rhs)
    {
        return lhs^rhs;
    }
	
	T identity()
	{
		return 0;
	}
};

template<typename T>
struct add
{
    T operator()(const T lhs, const T rhs)
    {
        return lhs+rhs;
    }
	
	T identity()
	{
		return 0;
	}
};

template<typename T>
struct mul
{
    T operator()(const T lhs, const T rhs)
    {
        return lhs*rhs;
    }
	
	T identity()
	{
		return 1;
	}
};

template<typename T, class Comparator>
struct min
{
    T operator()(const T lhs, const T rhs, inout Comparator comp)
    {
        return comp(lhs, rhs) ? lhs : rhs;
    }
};

template<typename T>
struct min
{
    T operator()(const T lhs, const T rhs)
    {
		comparator_lt_t comp;
        return bitwise_min(lhs, rhs, comp);
    }
	
	T identity()
	{
		return ~0;
	}
};

template<typename T, class Comparator>
struct max
{
    T operator()(const T lhs, const T rhs, inout Comparator comp)
    {
        return comp(lhs, rhs) ? lhs : rhs;
    }
};

template<typename T>
struct max
{
    T operator()(const T lhs, const T rhs)
    {
		comparator_gt_t comp;
        return bitwise_max(lhs, rhs, comp);
    }
	
	T identity()
	{
		return 0; // REVIEW: This assumes T = unsigned but what if we got T = signed ?
	}
};

template<typename T>
struct comparator_lt_t
{
    bool operator()(const T lhs, const T rhs)
    {
        return lhs<rhs;
    }
};

template<typename T>
struct comparator_gt_t
{
    bool operator()(const T lhs, const T rhs)
    {
        return lhs>rhs;
    }
};

template<typename T>
struct comparator_lte_t
{
    bool operator()(const T lhs, const T rhs)
    {
        return lhs<=rhs;
    }
};

template<typename T>
struct comparator_gte_t
{
    bool operator()(const T lhs, const T rhs)
    {
        return lhs>=rhs;
    }
};

}
}
}
#endif