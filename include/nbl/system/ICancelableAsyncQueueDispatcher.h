#ifndef __NBL_I_CANCELABLE_ASYNC_QUEUE_DISPATCHER_H_INCLUDED__
#define __NBL_I_CANCELABLE_ASYNC_QUEUE_DISPATCHER_H_INCLUDED__

#include "nbl/system/IAsyncQueueDispatcher.h"
#include "nbl/system/SReadWriteSpinLock.h"

namespace nbl {
namespace system
{

namespace impl
{
    class ICancelableAsyncQueueDispatcherBase
    {
    public:
        class future_base_t;

        struct request_base_t : impl::IAsyncQueueDispatcherBase::request_base_t
        {
            //! Atomically cancels this request
            bool set_cancel();
            bool query_cancel() const
            {
                return future == nullptr;
            }

        private:
            friend future_base_t;
            friend ICancelableAsyncQueueDispatcherBase;

            //! See ICancelableAsyncQueueDispatcher::associate_request_with_future() docs
            void associate_future_object(future_base_t* _future)
            {
                future = _future;
            }

            future_base_t* future = nullptr;
        };

        class future_base_t
        {
            friend request_base_t;

        protected:
            std::atomic_bool valid_flag = false;
            std::atomic<request_base_t*> request = nullptr;

            // future_t is non-copyable and non-movable
            future_base_t(const future_base_t&) = delete;

        public:
            future_base_t() = default;
            ~future_base_t()
            {
                request_base_t* req = request.exchange(nullptr);
                if (req)
                    req->set_cancel();
            }

            bool ready() const { return !request.load() || request.load()->ready; }
            bool valid() const { return valid_flag; }

            void wait()
            {
                if (!ready())
                    request.load()->wait();
            }
        };

    protected:
        void request_associate_future_object(request_base_t& req, future_base_t* future)
        {
            req.associate_future_object(future);
        }
    };
}

template <typename CRTP, typename RequestType, uint32_t BufferSize = 256u, typename InternalStateType = void>
class ICancelableAsyncQueueDispatcher : public IAsyncQueueDispatcher<CRTP, RequestType, BufferSize, InternalStateType>, public impl::ICancelableAsyncQueueDispatcherBase
{
    using this_async_queue_t = ICancelableAsyncQueueDispatcher<CRTP, RequestType, BufferSize, InternalStateType>;
    using base_t = IAsyncQueueDispatcher<CRTP, RequestType, BufferSize, InternalStateType>;
    friend base_t;

    template <typename T>
    class future_storage_t
    {
    public:
        alignas(T) uint8_t storage[sizeof(T)];

        T* getStorage() { return reinterpret_cast<T*>(storage); }
    };

public:
    using request_base_t = impl::ICancelableAsyncQueueDispatcherBase::request_base_t;

    static_assert(std::is_base_of_v<request_base_t, RequestType>, "Request type must derive from request_base_t!");

    template <typename T>
    class future_t : private future_storage_t<T>, public impl::ICancelableAsyncQueueDispatcherBase::future_base_t
    {
        friend this_async_queue_t;

        template <typename... Args>
        void notify(Args&&... args)
        {
            new (future_storage_t<T>::getStorage()) T(std::forward<Args>(args)...);
            valid_flag = true;
        }

        //! See ICancelableAsyncQueueDispatcher::associate_request_with_future() docs
        void associate_request(RequestType* req)
        {
            request = req;
        }

    public:
        using value_type = T;

        future_t() = default;

        ~future_t()
        {
            if (valid_flag)
                future_storage_t<T>::getStorage()->~T();
        }

        T& get()
        {
            future_base_t::wait();
            assert(valid_flag);
            T* ptr = future_storage_t<T>::getStorage();
            return (ptr[0]);
        }
    };

    using base_t::base_t;

protected:
    bool process_request_predicate(const RequestType& req)
    {
        return !req.query_cancel();
    }

    //! Must be called from within process_request()
    //! User is responsible for providing a value into the associated future object
    template <typename T, typename... Args>
    void notify_future(RequestType& req, Args&&... args)
    {
        auto& req_base = static_cast<request_base_t&>(req);
        auto* future = static_cast<future_t<T>*>(req.future);
        future->notify(std::forward(args)...);
    }

    //! Must be called from within request_impl()
    //! User is responsible for associating future object with a request
    //! Request is automatically cancelled if it is not associated with any future object
    //! More than one request associated with the same future object is undefined behaviour
    template <typename T>
    void associate_request_with_future(RequestType& req, future_t<T>& future)
    {
        assert(!future.valid());
        future_base_t* future_ptr = static_cast<future_base_t*>(&future);
        impl::ICancelableAsyncQueueDispatcherBase::request_associate_future_object(static_cast<request_base_t&>(req), future_ptr);
        future.associate_request(&req);
    }
};

inline bool impl::ICancelableAsyncQueueDispatcherBase::request_base_t::set_cancel()
{
    auto lk = lock();
    if (ready.load())
        return false;
    if (future)
        future->request = nullptr;
    future = nullptr;
    return true;
}

}}

#endif
