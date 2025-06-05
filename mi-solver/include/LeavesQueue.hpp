#ifndef _MI_QPIPMPC_LEAVESQUEUE_HPP_
#define _MI_QPIPMPC_LEAVESQUEUE_HPP_

#include "Node.hpp"
#include "QP_consts.hpp"
#include "DataVolatile.hpp"
#include <vector>
#include <algorithm>
#include <deque>
#include <cmath>
#include <utility>
#include <mutex>
#include <memory>

namespace MI_QPIPMPC
{

class DataVolatile; // forward declaration

class LeavesQueue
{
    public:
        void push(std::unique_ptr<Node> n);
        void push_top(std::unique_ptr<Node> n);
        void pop();
        std::unique_ptr<Node> pop_top(DataVolatile &data_volatile, int thread_number);
        void prune(double upper_glob);
        void clear();
        std::pair<double,bool> get_lower_glob() const;
        int size() const;
        bool empty() const;

    protected: 
        std::deque<std::unique_ptr<Node>> deque;

    private:
        mutable std::mutex mtx;
};

} // namespace MI_QPIPMPC

#endif