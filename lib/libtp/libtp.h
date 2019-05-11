/*********************************************************
*
*  Copyright (C) 2014 by Vitaliy Vitsentiy
*  Modifications copyright (C) 2019 by Minesh Patel (SAFARI Group at ETH Zurich)
*
*  Licensed under the Apache License, Version 2.0 (the "License");
*  you may not use this file except in compliance with the License.
*  You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
*  Unless required by applicable law or agreed to in writing, software
*  distributed under the License is distributed on an "AS IS" BASIS,
*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*  See the License for the specific language governing permissions and
*  limitations under the License.
*
*********************************************************/
#ifndef LIBTP_H
#define LIBTP_H

/*
 * This code is inspired by the CTPL library found at:
 *      https://github.com/vit-vit/CTPL
 * 
 * We license this file using the same Apache 2.0 license as the CTPL library 
 * in order to mark this heritage.
 */
#include <functional>
#include <thread>
#include <atomic>
#include <vector>
#include <memory>
#include <exception>
#include <future>
#include <mutex>
#include <queue>

class fn_obj
{
public:
    int p;
    std::function<void(int id)> *f;

    fn_obj(std::function<void(int id)> *fn_ptr, int priority) : p(priority), f(fn_ptr) {}
    ~fn_obj() {}

    bool operator<(const fn_obj& rhs) const
    {
        return p < rhs.p;
    }
};

class thread_pool 
{
public:

    thread_pool(int num_threads) 
    { 
        this->n_threads = num_threads;
        this->n_threads_running = 0;
        this->n_jobs_completed = 0;
        this->running = false;
        this->terminate = false;

        // launch N threads, which will wait for work
        auto thread_func = 
            [this](int tid)
            {
                std::unique_lock<std::mutex> q_lk(q_mx);
                // printf("Worker %d spawned\n", tid);

                while(1)
                {   
                    // wait for a new queue element to be available
                    while(q.empty() || running == false) 
                    {
                        q_cv.wait(q_lk);
                        if(terminate)
                            return;
                    }
                    fn_obj fobj = q.top();
                    q.pop();
                    n_threads_running++;
                    q_lk.unlock();

                    // perform dequeued work
                    // printf("Worker %d allocated work\n", tid);
                    (*fobj.f)(tid);
                    delete fobj.f;
                    // printf("Worker %d finished work\n", tid);

                    q_lk.lock();
                    n_threads_running--;
                    n_jobs_completed++;
                    thread_cv.notify_one();                 
                }
                // printf("Worker %d ABOUT TO EXIT\n", tid);
            };

        // launch n_threads
        for (int i = 0; i < num_threads; ++i)
            threads.push_back(std::thread(thread_func, i));
    }

    // the destructor waits for all jobs to finish
    ~thread_pool() 
    {
        this->wait();

        std::unique_lock<std::mutex> lk(q_mx);
        this->terminate = true;
        q_cv.notify_all();
        lk.unlock();

        for (int i = 0; i < this->n_threads; ++i)
            threads[i].join(); // cleanly exit
    }

    // get the number of running threads in the pool
    int get_n_jobs_completed(void) 
    { 
        std::unique_lock<std::mutex> lk(q_mx);
        size_t ret = n_jobs_completed;
        lk.unlock();

        return (int)ret;
    }

    int get_n_jobs_outstanding(void) 
    { 
        std::unique_lock<std::mutex> lk(q_mx);
        size_t ret = q.size() + n_threads_running;
        lk.unlock();

        return (int)ret;
    }

    // wait for all computing threads to finish and stop all threads
    // may be called asynchronously to not pause the calling thread while waiting
    // if isWait == true, all the functions in the queue are run, otherwise the queue is cleared without running the functions
    void wait(bool pause = false) 
    {
        // wait on condvar for all thrads to be completed
        std::unique_lock<std::mutex> lk(q_mx);
        if(pause)
        {
            running = false;
            while(this->n_threads_running) 
                thread_cv.wait(lk);
        }
        else
        {
            while(this->n_threads_running || !q.empty()) 
                thread_cv.wait(lk);         
        }
    }

    void start(void)
    {
        std::unique_lock<std::mutex> lk(q_mx);
        running = true;
        q_cv.notify_all();
    }

    // returns std::future<>
    // this function is modeled from CTPL (https://github.com/vit-vit/CTPL)
    template<typename F, typename... Args>
    auto add(F &&f, int priority, Args &&... rest) -> std::future<decltype(f(0, rest...))> 
    {
        auto packed_fn = std::make_shared<std::packaged_task<decltype(f(0, rest...))(int)>>
            (std::bind(std::forward<F>(f), std::placeholders::_1, std::forward<Args>(rest)...));
        auto fn_wrapper = new std::function<void(int id)>
        (
            [packed_fn](int id) 
            {
                //printf("Thread %d about to launch: %p!\n", id, packed_fn.get());
                (*packed_fn)(id);
            }
        );
        fn_obj fobj(fn_wrapper, priority);

        std::unique_lock<std::mutex> lk(this->q_mx);
        q.push(fobj);
        q_cv.notify_one(); // tell the pool that there's work available!
        return packed_fn->get_future();
    }

    void reset_stats(void)
    {
        n_jobs_completed = 0;
    }

private:
    thread_pool(const thread_pool &);
    thread_pool(thread_pool &&);
    thread_pool & operator=(const thread_pool &);
    thread_pool & operator=(thread_pool &&);

    std::priority_queue<fn_obj> q;
    std::mutex q_mx;
    std::condition_variable q_cv;
    std::condition_variable thread_cv;
    std::vector< std::thread > threads;

    int n_threads;
    int n_threads_running;
    int n_jobs_completed;
    bool running;
    bool terminate;
};


#endif /* LIBTP_H */
