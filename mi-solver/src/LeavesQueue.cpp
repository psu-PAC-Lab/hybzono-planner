#include "LeavesQueue.hpp"

using namespace MI_QPIPMPC;

void LeavesQueue::push(std::unique_ptr<Node> n)
{
    std::lock_guard<std::mutex> lock(mtx);

    if (deque.empty())
    {
        deque.push_back(std::move(n));
    }
    else if (*n < *(deque.front()))
    {
        deque.push_front(std::move(n));
    }
    else if (!(*n < *(deque.back())))
    {
        deque.push_back(std::move(n));
    }
    else 
    {
        // check whether n is closer to front or back of queue
        if (abs(*n - *(deque.front())) < abs(*n - *(deque.back())))
        {
            // search for insertion point from front
            for (auto it=deque.begin(); it!=deque.end(); ++it)
            {
                if (*n < **it)
                {
                    deque.insert(it, std::move(n));
                    return;
                }
            }
        }
        else
        {
            // search for insertion point from back
            for (auto it=deque.rbegin(); it!=deque.rend(); ++it)
            {
                if (!(*n < **it))
                {
                    deque.insert(it.base(), std::move(n));
                    return;
                }
            }
        }
    }
}

void LeavesQueue::push_top(std::unique_ptr<Node> n)
{
    std::lock_guard<std::mutex> lock(mtx);
    deque.push_back(std::move(n));
}

void LeavesQueue::pop()
{
    std::lock_guard<std::mutex> lock(mtx);
    deque.pop_back();
}

std::unique_ptr<Node> LeavesQueue::pop_top(DataVolatile &data_volatile, int thread_number)
{
    std::lock_guard<std::mutex> lock(mtx);

    // protect against empty queue
    if (deque.empty())
    {
        data_volatile.reset_lower_glob_threads(thread_number);
        return std::unique_ptr<Node> (new Node);
    }
    else
    {
        // get node
        std::unique_ptr<Node> n = std::move(deque.back());
        deque.pop_back();

        // set lower_glob vector in data_volatile
        data_volatile.set_lower_glob_threads(thread_number, n->lower);

        // return
        return std::move(n);
    }
}

void LeavesQueue::prune(double upper_glob)
{
    std::lock_guard<std::mutex> lock(mtx);

    // remove leaves with lower bound above upper_glob
    while (!deque.empty())
    {
        if (deque.front()->lower > upper_glob)
            deque.pop_front();
        else
            break;
    }
}

void LeavesQueue::clear()
{
    std::lock_guard<std::mutex> lock(mtx);
    deque.clear();
}

std::pair<double, bool> LeavesQueue::get_lower_glob() const
{
    std::lock_guard<std::mutex> lock(mtx);

    // protect against empty queue
    if (deque.empty())
        return std::make_pair(QP::inf, false); // invalid, empty case
    
    // init
    auto it = deque.rbegin();
    double lower = QP::inf;

    // find lowest lower bound
    while (it != deque.rend() && (**it).is_priority())
    {
        lower = std::min((**it).lower, lower);
        ++it;
    }

    // return
    if (it == deque.rend())
        return std::make_pair(lower, true); // valid, all priority cases
    else
        return std::make_pair(std::min((**it).lower, lower), true); // valid, nominal case
}

int LeavesQueue::size() const
{
    std::lock_guard<std::mutex> lock(mtx);
    return deque.size();
}

bool LeavesQueue::empty() const
{
    std::lock_guard<std::mutex> lock(mtx);
    return deque.empty();
}