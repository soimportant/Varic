#pragma once

#include <vector>
#include <map>
#include <memory>
#include <concepts>
#include <exception>
#include <algorithm>
#include "worker.hpp"

namespace biovoltron {

namespace detail {
  namespace policy {
  /**
   * @brief RoundRobin policy that used to choose approriate worker when a new
   * task has been submit into threadpool
   */
  class RoundRobin {
  private:
    using IdType = std::size_t;
    std::atomic<IdType> id_counter { 0 };

    auto get_next_id(const std::size_t max_id) -> IdType {
      auto old_id = id_counter.load(std::memory_order_acquire);
      auto new_id = old_id + 1;
      if (new_id >= max_id) {
        new_id = 0;
      }
      while (!id_counter.compare_exchange_weak(old_id,
                                               new_id));
      // ? survey compare_exchange_strong and compare_exchange_weak

      return new_id;
    }

  public:
    auto operator()(const std::vector<Worker>& workers) {
      return get_next_id(workers.size());
    }
  };

  class Random {
  private:
    using IdType = std::size_t;

    auto get_next_id(const std::size_t max_id) -> IdType {
      thread_local auto rng = std::mt19937_64(std::random_device{}());
      return rng() % max_id;
    }
  public:
    auto operator()(const std::vector<Worker>& workers) {
      return get_next_id(workers.size());
    }
  };

  /**
   * ! QueueSize policy failed the test, it will make a worker in waiting status
   * ! even the queue is empty
   */
  // class QueueSize {
  // public:
  //   auto operator()(const std::vector<Worker>& workers) {
  //     auto iter = std::ranges::min_element(workers, [](auto& a, auto& b) {
  //       return a.size() < b.size();
  //     });
  //     return std::ranges::distance(workers.begin(), iter);
  //   }
  // };
  }
}
  

template<class P>
concept WorkerSelectable = std::convertible_to<
  std::invoke_result_t<P, std::vector<Worker>&>, std::size_t
>;

template<WorkerSelectable Policy = detail::policy::RoundRobin>
class ThreadPool {
public:
  /* disable copy */
  ThreadPool(const ThreadPool&) = delete;
  ThreadPool& operator =(const ThreadPool&) = delete;

  /* disable move */
  ThreadPool(ThreadPool&&) = delete;
  ThreadPool& operator =(ThreadPool&&) = delete;

  ThreadPool(const std::size_t& size) : pool_size(size), workers(size) {}

protected:
  /**
   * @brief the number of threads in threadpool
   */
  std::size_t pool_size;
  /**
   * @brief the workers -> the threads
   */
  std::vector<Worker> workers;

public:
  /**
   * @brief submit a callable object and its arguments to a specific worker 
   * inside the threadpool.
   * 
   * @tparam F the type of callable object
   * @tparam ArgType the type of arguments of the callable object
   * @param id the worker id
   * @param func the callable object
   * @param args the arguments of `func`
   * @return a pair contains the worker id which and std::future of the task
   * @note If you want to pass the argument as a reference, use std::ref to
   * wrap your argument.
   */
  template<class F, class... ArgType>
  auto submit(const std::size_t& id, F&& func, ArgType&&... args) 
    requires std::invocable<F, ArgType...> {
    auto res = workers[id].push(std::forward<F>(func),
                                std::forward<ArgType>(args)...);
    return std::make_pair(id, std::move(res));
  }

  /**
   * @brief submit a callable object and its argument to the threadpool, the
   * worker is been determined by `Policy`
   * 
   * @tparam F the type of callable object
   * @tparam ArgType the type of arguments of the callable object
   * @param id the worker id
   * @param func the callable object
   * @param args the arguments of `func`
   * @return a `std::pair` contains the worker id and std::future of the task
   * @note If you want to pass the argument as a reference, use `std::ref` to
   * wrap your argument.
   */
  template<class F, class... ArgType>
  auto submit(F&& func, ArgType&&... args) 
    requires std::invocable<F, ArgType...> {
    static Policy p;
    const std::size_t worker_id = p(workers);
    return submit(worker_id, std::forward<F>(func),
                             std::forward<ArgType>(args)...);
  }

  /**
   * The request stop function need to redesign, and deal with the following 
   * scenario appropriately.
   * 1. call request stop on same id twice
   *   - std::stop_token maybe helpful
   * 2. If user submit a job, we cannot submit the job to the stopped worker
   *   - need to handle a data structure to store stopped id
   *   - or throw a exception when `Policy` selected a stopped id
   */
  // /**
  //  * @brief make a stop request for a specific worker
  //  * @param id the worker id
  //  * @return `true` if this invocation made a stop request, `false` otherwise
  //  */
  // auto request_stop(int id) const noexcept {
  //   return workers[id].request_stop();
  // }

  /**
   * @brief the number of threads in threadpool
   * @return std::size_t 
   */
  auto size() const noexcept {
    return pool_size;
  }
};


/**
 * @brief An api for making thread pool
 * @tparam Policy the class that used to choose the appropriate worker to run 
 * the task, default is biovoltron::detail::Policy::RoundRobin.
 * @param threads how many threads in the threadpool, default is 
 * `std::thread::hardware_concurrency`
 * @return a thread pool with `threads` workers
 */
template<WorkerSelectable Policy = detail::policy::RoundRobin>
auto make_threadpool(
  const std::size_t threads = std::thread::hardware_concurrency()) {
  return ThreadPool<Policy>(threads);
}

} // namespace biovoltron
