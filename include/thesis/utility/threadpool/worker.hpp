#pragma once

#include <functional>
#include <future>
#include <chrono>
#include <mutex>
#include <memory>
#include <thread>
#include <concepts> 
#include <queue>
#include <cassert>

using namespace std::literals;

namespace biovoltron {

class Worker {

protected:
  using Thread = std::jthread;
  template<class T>
  using Queue = std::queue<T>;
  using TaskType = std::packaged_task<void(void)>;
  using TaskQueue = Queue<TaskType>;

protected:
  /**
   * @brief whether the worker is running
   */
  std::atomic_bool running { false };
  /**
   * @brief whether the worker is waiting
   */
  std::atomic_bool waiting { false };
  /**
   * @brief the pointer of underlying thread
   */
  std::unique_ptr<Thread> worker { nullptr };
  /**
   * @brief the task queue of this worker
   */
  TaskQueue task_queue;
  /**
   * @brief the lock for task_queue
   */
  std::mutex queue_mutex;
  /**
   * @brief the conditional variable to wake up the worker when a new job is
   * pushed into empty task_queue
   */
  std::condition_variable_any run_cv;

  /**
   * @brief the maximum waiting time for `run_cv`
   */
  std::chrono::milliseconds stop_duration = 100ms;

public:
  /* disable copy */
  Worker(const Worker&) = delete;
  Worker& operator=(const Worker&) = delete;

  /* disable move */
  Worker(Worker&&) = delete;
  Worker& operator=(Worker&&) = delete;

  Worker() {
    /**
     * load a underlying thread
     */
    worker.reset(new Thread([this](std::stop_token stop_token) {
      auto m = std::mutex {};
      auto lock = std::unique_lock { m };
      while (!stop_token.stop_requested()) {
        if (!run_single_task()) {
          running = false;
          /**
           * there's no job to do
           * 1. Notify the join conditional variable
           * 2. wait until a new job
           */
          run_cv.wait_for(lock, stop_token, stop_duration, [this]() {
            return !task_queue.empty();
          });
          running = true;
        }
      }
    }));
  }

protected:
  /**
   * @brief consume a task from the queue and run it
   * @return `true` if worker consume a job and run it successively, `false`
   * otherwise
   */
  auto run_single_task() -> bool {
    auto task = TaskType();
    if (auto lock = std::scoped_lock { queue_mutex }; !task_queue.empty()) {
      task = std::move(task_queue.front());
      task_queue.pop();
    } else {
      return false;
    }
    task();
    return true;
  }

public:

  /**
   * @brief push a callable object and its arguments to worker queue
   * @tparam F the type of callable object
   * @tparam Args... the type of arguments of callable objects
   * @param func callable object
   * @param args the arguments of `func`
   * @return the std::future of the job.
   */
  template<class F, class... Args>
  auto push(F&& func, Args&&... args) {

    /* pack the task */
    using ret_type = std::invoke_result_t<F, Args...>;
    auto task = std::packaged_task<ret_type()>(
      std::bind(std::forward<F>(func), std::forward<Args>(args)...)
    );
    auto res = task.get_future();

    /**
     * packaged_task magic!
     * see the move constructor of packaged_task at
     * https://en.cppreference.com/w/cpp/thread/packaged_task/packaged_task
     *
     * here we use `pakcaged_task<void(void)>` as a wrapper to wrap
     * packaged_task in any type, then push it into task queue.
     */
    auto wrapper = std::packaged_task<void(void)>(std::move(task));

    std::scoped_lock guard { queue_mutex };
    task_queue.emplace(std::move(wrapper));
    /**
     * notify the waiting worker to start working
     */
    run_cv.notify_one();
    return res;
  }

  /**
   * @brief whether the worker is running or not
   * @return `true` if the worker is currently running, `false` otherwise.
   */
  auto is_running() const noexcept -> bool {
    return running;
  }

  /**
   * @brief the number of task inside the queue, which can be used to measure
   * the loading of worker
   * @return the size of worker queue
   */
  auto size() const noexcept {
    return task_queue.size();
  }

  auto get_id() const noexcept {
    return worker->get_id();
  }

  auto join() {
    worker->request_stop();
    run_cv.notify_one();
    worker->join();
  }

  /**
   * @brief made a request to the worker for stop working
   * @return `true` if this invocation made a stop request or it doesn't running,
   * `false` otherwise
   */
  auto request_stop() const noexcept {
    if (!running) {
      return true;
    }
    return worker->request_stop();
  }

  auto stop_requested() const noexcept {
    return worker->get_stop_token().stop_requested();
  }

  ~Worker() {
    worker->request_stop();
    run_cv.notify_one();
  }
};

} // namespace biovoltron




// template<bool CommonQueue>
// class Worker {
// protected:
//   using Thread = std::jthread;
//   template<class T>
//   using Queue = std::queue<T>;
//   using TaskType = std::packaged_task<void(void)>;
//   using TaskQueue = Queue<TaskType>;
//   template<class T>
//   using PtrType = std::conditional_t<CommonQueue,
//                                      std::shared_ptr<T>,
//                                      std::unique_ptr<T>>;

// protected:
//   bool running;
//   std::unique_ptr<Thread> worker;
//   // std::unique_ptr<TaskQueue> task_queue;
//   // std::mutex queue_mutex;

//   PtrType<TaskQueue> task_queue;
//   PtrType<std::mutex> queue_mutex;
//   PtrType<std::condition_variable> run_cv;
// public:
//   /* disable copy */
//   Worker(const Worker&) = delete;
//   Worker& operator=(const Worker&) = delete;

//   /* disable move */
//   Worker(Worker&&) = delete;
//   Worker& operator=(Worker&&) = delete;

//   Worker() : worker(nullptr) {
//     if constexpr (CommonQueue) {
//       static auto common_cv = std::make_shared<std::condition_variable>();
//       run_cv = std::shared_ptr(common_cv);

//       static auto queue = std::make_shared<TaskQueue>();
//       task_queue = std::shared_ptr(queue);

//       static auto mutex = std::make_shared<std::mutex>();
//       queue_mutex = std::shared_ptr(mutex);
//       load();
//     } else {
//       run_cv = std::make_unique<std::condition_variable>();
//       task_queue = std::make_unique<TaskQueue>();
//       queue_mutex = std::make_unique<std::mutex>();
//     }
//   }

// protected:
//   /**
//    * @brief consume a task from the queue and run it
//    * @return `true` if worker consume a job and run it successively, `false`
//    * otherwise
//    */
//   auto run_single_task() {
//     queue_mutex->lock();
//     if (!task_queue->empty()) {
//       TaskType job = std::move(task_queue->front());
//       task_queue->pop();
//       queue_mutex->unlock();
//       job();
//       return true;
//     }
//     queue_mutex->unlock();
//     return false;
//   }

// public:
//   /**
//    * @brief push a callable object and its arguments to worker queue
//    * @tparam F the type of callable object
//    * @tparam Args... the type of arguments of callable objects
//    * @param func callable object
//    * @param args the arguments of `func`
//    * @return the std::future of the job.
//    */
//   template<class F, class... Args>
//   auto push(F&& func, Args&&... args)
//     requires std::invocable<F, Args...> {

//     /* pack the task */
//     using ret_type = std::invoke_result_t<F, Args...>;
//     auto task = std::packaged_task<ret_type()>(
//       std::bind(std::forward<F>(func), std::forward<Args...>(args...))
//     );
//     auto res = task.get_future();

//     /**
//      * packaged_task magic!
//      * see the move constructor of packaged_task at
//      * https://en.cppreference.com/w/cpp/thread/packaged_task/packaged_task
//      *
//      * here we use `pakcaged_task<void(void)>` as a wrapper to wrap
//      * packaged_task in any type, then push it into task queue.
//      */
//     auto wrapper = std::packaged_task<void(void)>(std::move(task));

//     queue_mutex->lock();
//     task_queue->emplace(std::move(wrapper));
//     queue_mutex->unlock();

//     /**
//      * notify the waiting worker to start working
//      */
//     run_cv->notify_one();
//     if constexpr (!CommonQueue) {
//       if (worker == nullptr) {
//         load();
//       }
//     }
//     return res;
//   }

  // /**
  //  * @brief load a underlying thread
  //  * @return None
  //  */
  // auto load( /* std::unique_lock<std::mutex>& */) -> void {
  //   worker.reset(new Thread([this](std::stop_token stop_token) {
  //     while (!stop_token.stop_requested()) {
  //       if (!run_single_task()) {
  //         running = false;
  //         /**
  //          * there's no job to do
  //          * 1. Notify the join conditional variable
  //          * 2. wait until a new job
  //          */
  //         // join_cv.notify_one();

  //         // queue_mutex->unlock();

  //         std::mutex m;
  //         std::unique_lock<std::mutex> lock(m);
  //         run_cv->wait(lock);
  //         running = true;
  //       }
  //     }
  //   }));
  // }

  // /**
  //  * @brief whether the worker is running or not
  //  * @return `true` if this thread is currently running, `false` otherwise.
  //  */
  // auto is_running() {
  //   return running;
  // }

  // /**
  //  * @brief made a request to the worker for stop working
  //  * @return `true` if this invocation made a stop request or it doesn't running,
  //  * `false` otherwise
  //  */
  // auto request_stop() noexcept {
  //   if (!running) {
  //     return true;
  //   }
  //   return worker->request_stop();
  // }

//   ~Worker() {
//     worker->request_stop();
//     run_cv->notify_all();

//     // segmentation fault in destructor
//     if constexpr (CommonQueue) {

//     }
//   }
// };

// }