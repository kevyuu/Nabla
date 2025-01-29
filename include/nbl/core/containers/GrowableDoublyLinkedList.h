// Copyright (C) 2018-2020 - DevSH Graphics Programming Sp. z O.O.
// This file is part of the "Nabla Engine".
// For conditions of distribution and use, see copyright notice in nabla.h

#ifndef __NBL_CORE_GROWABLE_DOUBLY_LINKED_LIST_H_INCLUDED__
#define __NBL_CORE_GROWABLE_DOUBLY_LINKED_LIST_H_INCLUDED__


#include "nbl/core/containers/lists/DoublyLinkedListBase.h"

namespace nbl
{
namespace core
{

template<typename Value>
class GrowableDoublyLinkedList : public DoublyLinkedListBase<Value>
{
public:
	using base_t = DoublyLinkedListBase<Value>;
	using value_t = Value;
	using node_t = typename base_t::node_t;
	using address_allocator_t = typename base_t::address_allocator_t;
	using disposal_func_t = typename base_t::disposal_func_t;

	//Constructor, capacity determines the amount of allocated space
	GrowableDoublyLinkedList(const uint32_t capacity, disposal_func_t&& dispose_f = disposal_func_t()) : base_t(capacity, std::move(dispose_f))
	{}

	GrowableDoublyLinkedList() = default;

	GrowableDoublyLinkedList(const GrowableDoublyLinkedList& other) = delete;

	GrowableDoublyLinkedList& operator=(const GrowableDoublyLinkedList& other) = delete;

	GrowableDoublyLinkedList& operator=(GrowableDoublyLinkedList&& other)
	{
		base_t::operator=(other);
	}

	~GrowableDoublyLinkedList() = default;

	/**
	* @brief Resizes the list by extending its capacity so it can hold more elements. Returns a bool indicating if capacity was indeed increased.
	*
	* @param [in] newCapacity New number of elements to hold. MUST be greater than current list capacity.
	*/
	inline bool grow(uint32_t newCapacity)
	{
		// Must at least make list grow
		if (newCapacity <= this->m_cap)
			return false;
		// Same as code found in ContiguousMemoryLinkedListBase to create aligned space
		const auto firstPart = core::alignUp(address_allocator_t::reserved_size(1u, newCapacity, 1u), alignof(node_t));
		void* newReservedSpace = _NBL_ALIGNED_MALLOC(firstPart + newCapacity * sizeof(node_t), alignof(node_t));

		// Malloc failed, not possible to grow
		if (!newReservedSpace)
			return false;

		node_t* newArray = reinterpret_cast<node_t*>(reinterpret_cast<uint8_t*>(newReservedSpace) + firstPart);

		// Copy memory over to new buffer, then free old one
		memcpy(reinterpret_cast<void*>(newArray), reinterpret_cast<void*>(this->m_array), this->m_cap * sizeof(node_t));
		_NBL_ALIGNED_FREE(this->m_reservedSpace);

		// Finally, create new address allocator from old one
		this->m_addressAllocator = std::unique_ptr<address_allocator_t>(new address_allocator_t(newCapacity, std::move(*(this->m_addressAllocator)), newReservedSpace));
		this->m_cap = newCapacity;
		this->m_array = newArray;
		this->m_reservedSpace = newReservedSpace;

		return true;
	}
};


}
}


#endif
