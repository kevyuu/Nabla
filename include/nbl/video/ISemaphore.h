#ifndef _NBL_VIDEO_I_SEMAPHORE_H_INCLUDED_
#define _NBL_VIDEO_I_SEMAPHORE_H_INCLUDED_


#include "nbl/core/IReferenceCounted.h"

#include <chrono>

#include "nbl/video/decl/IBackendObject.h"


namespace nbl::video
{

class ISemaphore : public IBackendObject
{
    public:
        // basically a pool function
        virtual uint64_t getCounterValue() const = 0;

        //! Basically the counter can only monotonically increase with time (ergo the "timeline"):
        // 1. `value` must have a value greater than the current value of the semaphore (what you'd get from `getCounterValue()`)
        // 2. `value` must be less than the value of any pending semaphore signal operations (this is actually more complicated)
        // Current pending signal operations can complete in any order, unless there's an execution dependency between them,
        // this will change the current value of the semaphore. Consider a semaphore with current value of 2 and pending signals of 3,4,5;
        // without any execution dependencies, you can only signal a value higher than 2 but less than 3 which is impossible.
        virtual void signal(const uint64_t value) = 0;

        // We don't provide waits as part of the semaphore (cause you can await multiple at once with ILogicalDevice),
        // but don't want to pollute ILogicalDevice with lots of enums and structs
        struct SWaitInfo
        {
            const ISemaphore* semaphore = nullptr;
            uint64_t value = 0;
        };
        enum class WAIT_RESULT : uint8_t
        {
            TIMEOUT,
            SUCCESS,
            DEVICE_LOST,
            _ERROR
        };

        // Utility structs
        class future_base_t
        {
            public:
                inline void set(const SWaitInfo& wait)
                {
                    m_semaphore = core::smart_refctd_ptr<const ISemaphore>(wait.semaphore);
                    m_waitValue = wait.value;
                }
                inline bool ready() const
                {
                    return m_semaphore->getCounterValue()>=m_waitValue;
                }
                WAIT_RESULT wait() const;

            protected:
                // the base class is not directly usable
                inline future_base_t() = default;
                // derived won't be copyable
                future_base_t(const future_base_t&) = delete;
                inline future_base_t(future_base_t&&) = default;
                future_base_t& operator=(const future_base_t&) = delete;
                inline future_base_t& operator=(future_base_t&&) = default;

                // smartpointer cause lifetime needs to be maintained
                core::smart_refctd_ptr<const ISemaphore> m_semaphore;
                uint64_t m_waitValue;
        };
        // Similar to `ISystem::future_t` but NOT thread-safe and the condition we wait on is signalled by the Device
        template<typename T>
        class future_t : private core::StorageTrivializer<T>, public future_base_t
        {
                using storage_t = core::StorageTrivializer<T>;
                using this_t = future_t<T>;

            public:
                template<typename... Args>
                inline future_t(Args&&... args)
                {
                    storage_t::construct(std::forward(args)...);
                }
                inline future_t(this_t&& other) : future_base_t(std::move<future_base_t>(other))
                {
                    if constexpr (!std::is_void_v<T>)
                        storage_t::construct(std::move(*other.getStorage()));
                }
                inline ~future_t()
                {
                    const bool success = wait();
                    assert(success);
                    storage_t::destruct();
                }

                inline this_t& operator=(this_t&& rhs)
                {
                    future_base_t::operator=(std::move<future_base_t>(rhs));
                    if constexpr (!std::is_void_v<T>)
                        *storage_t::getStorage() = *rhs.getStorage();
                    return *this;
                }

                inline const T* get() const
                {
                    if (ready())
                        return storage_t::getStorage();
                    return nullptr;
                }

                template<std::copyable U> requires std::is_same_v<T,U>
                inline T copy() const
                {
                    const bool success = wait();
                    assert(success);
                    return *get();
                }

                template<std::movable U> requires std::is_same_v<T,U>
                inline void move_into(U& dst)
                {
                    dst = std::move(*get());
                }
        };

        // Vulkan: const VkSemaphore*
        virtual const void* getNativeHandle() const = 0;

    protected:
        inline ISemaphore(core::smart_refctd_ptr<const ILogicalDevice>&& dev) : IBackendObject(std::move(dev)) {}
        virtual ~ISemaphore() = default;
};

}
#endif