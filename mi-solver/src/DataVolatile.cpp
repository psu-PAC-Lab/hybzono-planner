#include "DataVolatile.hpp"

using namespace MI_QPIPMPC;

// constructor
DataVolatile::DataVolatile()
{
    std::lock_guard<std::mutex> lock(mtx);
    first_run = true;
    iter_num = 0;
    qp_solve_time = 0.0;
    qp_iter = 0;
    qp_iter_avg = 0;
    lower_glob = -QP::inf;
    upper_glob = QP::inf;
    x = Eigen::VectorXd(0);
    x0 = Eigen::VectorXd(0);
    status = MIQP_NO_SOL;
    setup_time = 0.0;
    solve_time = 0.0;
    run_time = 0.0;
}

DataVolatile::DataVolatile(const Data& data)
{
    std::lock_guard<std::mutex> lock(mtx);
    first_run = true;
    iter_num = 0;
    qp_solve_time = 0.0;
    qp_iter = 0;
    qp_iter_avg = 0;
    lower_glob = -QP::inf;
    upper_glob = QP::inf;
    x = Eigen::VectorXd::Zero(data.n);
    x0 = Eigen::VectorXd::Zero(data.n);
    status = MIQP_NO_SOL;
    setup_time = 0.0;
    solve_time = 0.0;
    run_time = 0.0;
}

// move constructor
DataVolatile::DataVolatile(DataVolatile &&other)
{
    std::lock_guard<std::mutex> lock(mtx);
    first_run = other.first_run;
    iter_num = other.iter_num;
    qp_solve_time = other.qp_solve_time;
    qp_iter = other.qp_iter;
    qp_iter_avg = other.qp_iter_avg;
    lower_glob = other.lower_glob;
    upper_glob = other.upper_glob;
    x = std::move(other.x);
    x0 = std::move(other.x0);
    status = other.status;
    setup_time = other.setup_time;
    solve_time = other.solve_time;
    run_time = other.run_time;
    region_vec = std::move(other.region_vec);
}

DataVolatile& DataVolatile::operator=(const DataVolatile &other)
{
    std::lock_guard<std::mutex> lock(mtx);
    first_run = other.first_run;
    iter_num = other.iter_num;
    qp_solve_time = other.qp_solve_time;
    qp_iter = other.qp_iter;
    qp_iter_avg = other.qp_iter_avg;
    lower_glob = other.lower_glob;
    upper_glob = other.upper_glob;
    x = other.x;
    x0 = other.x0;
    status = other.status;
    setup_time = other.setup_time;
    solve_time = other.solve_time;
    run_time = other.run_time;
    region_vec = other.region_vec;

    return *this;
}

// get methods
int DataVolatile::get_iter_num() const
{
    std::lock_guard<std::mutex> lock(mtx);
    return iter_num;
}

double DataVolatile::get_lower_glob() const
{
    std::lock_guard<std::mutex> lock(mtx);

    return lower_glob;
}

double DataVolatile::get_upper_glob() const
{
    std::lock_guard<std::mutex> lock(mtx);
    return upper_glob;
}

MIQP_Status DataVolatile::get_status() const
{
    std::lock_guard<std::mutex> lock(mtx);
    return status;
}

int DataVolatile::get_qp_iter() const
{
    std::lock_guard<std::mutex> lock(mtx);
    return qp_iter;
}

double DataVolatile::get_run_time() const
{
    std::lock_guard<std::mutex> lock(mtx);
    return run_time;
}

Eigen::VectorXd DataVolatile::get_x() const
{
    std::lock_guard<std::mutex> lock(mtx);
    return x;
}

double DataVolatile::get_qp_solve_time() const
{
    std::lock_guard<std::mutex> lock(mtx);
    return qp_solve_time;
}

int DataVolatile::get_qp_iter_avg() const
{
    std::lock_guard<std::mutex> lock(mtx);
    return qp_iter_avg;
}

std::vector<int> DataVolatile::get_region_vec() const
{
    std::lock_guard<std::mutex> lock(mtx);
    return region_vec;
}

// set methods
void DataVolatile::set_setup_time(double setup_time)
{
    // lock guard
    std::lock_guard<std::mutex> lock(mtx);
    this->setup_time = setup_time;

}

void DataVolatile::set_status(MIQP_Status status)
{
    // lock guard
    std::lock_guard<std::mutex> lock(mtx);
    this->status = status;
}

void DataVolatile::set_lower_glob(std::pair<double, bool> lower_glob_valid)
{
    // lock guard
    std::lock_guard<std::mutex> lock(mtx);

    // get lower bound from threads
    double lower_thread = *std::min_element(lower_glob_threads.begin(), lower_glob_threads.end());

    // get lower bound from leaves
    if (lower_glob_valid.second) // valid lower bound
        this->lower_glob = std::min(lower_glob_valid.first, lower_thread);
    else
        this->lower_glob = lower_thread;
}

void DataVolatile::set_upper_glob(double upper_glob)
{
    // lock guard
    std::lock_guard<std::mutex> lock(mtx);
    this->upper_glob = upper_glob;
}

void DataVolatile::set_solve_time(double solve_time)
{
    // lock guard
    std::lock_guard<std::mutex> lock(mtx);
    this->solve_time = solve_time;
}

void DataVolatile::update_optimum(const std::unique_ptr<Node> &leaf)
{
    std::lock_guard<std::mutex> lock(mtx);
    if (leaf->obj < upper_glob)
    {
        x = leaf->x;
        upper_glob = leaf->obj;

        this->region_vec.clear();
        for (auto it = leaf->region_vec.begin(); it != leaf->region_vec.end(); it++)
            this->region_vec.push_back(*(it->second.begin()));
    }
}

void DataVolatile::init_lower_glob_threads(int num_threads)
{
    std::lock_guard<std::mutex> lock(mtx);
    lower_glob_threads = std::vector<double>(num_threads, QP::inf);
}

void DataVolatile::set_lower_glob_threads(int thread_number, double lower_glob)
{
    std::lock_guard<std::mutex> lock(mtx);
    lower_glob_threads[thread_number] = lower_glob;
}

void DataVolatile::reset_lower_glob_threads(int thread_number)
{
    std::lock_guard<std::mutex> lock(mtx);
    lower_glob_threads[thread_number] = QP::inf;
}

void DataVolatile::increment_iter_num()
{
    std::lock_guard<std::mutex> lock(mtx);
    iter_num++;
}

void DataVolatile::increment_qp_iter(int n)
{
    std::lock_guard<std::mutex> lock(mtx);
    qp_iter += n;
}

void DataVolatile::increment_qp_solve_time(double t)
{
    std::lock_guard<std::mutex> lock(mtx);
    qp_solve_time += t;
}

void DataVolatile::calc_qp_iter_avg()
{
    std::lock_guard<std::mutex> lock(mtx);
    if (iter_num != 0)
        qp_iter_avg = qp_iter / iter_num;
    else
        qp_iter_avg = 0;
}

void DataVolatile::calc_run_time()
{
    std::lock_guard<std::mutex> lock(mtx);

    // calculate run time
    if (first_run)
    {
        run_time = setup_time + solve_time;
        first_run = false;
    }
    else
    {
        run_time = solve_time;
    }
}
    
// reset
void DataVolatile::reset(const Data& data)
{
    std::lock_guard<std::mutex> lock(mtx);
    first_run = true;
    iter_num = 0;
    qp_solve_time = 0.0;
    qp_iter = 0;
    qp_iter_avg = 0;
    lower_glob = -QP::inf;
    upper_glob = QP::inf;
    x = Eigen::VectorXd::Zero(data.n);
    x0 = Eigen::VectorXd::Zero(data.n);
    status = MIQP_NO_SOL;
    setup_time = 0.0;
    solve_time = 0.0;
    run_time = 0.0;
    region_vec.clear();
}