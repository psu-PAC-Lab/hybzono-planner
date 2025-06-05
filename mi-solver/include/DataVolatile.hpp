#ifndef _MI_QPIPMPC_DATAVOLATILE_HPP_
#define _MI_QPIPMPC_DATAVOLATILE_HPP_

#include "Eigen/Dense"
#include "Eigen/Sparse"
#include <map>
#include <vector>
#include <mutex>
#include <memory>
#include "Node.hpp"
#include "Data.hpp"
#include "LeavesQueue.hpp"
#include "QP_consts.hpp"

namespace MI_QPIPMPC
{

class DataVolatile
{
    public:
        DataVolatile();
        DataVolatile(const Data& data);
        DataVolatile(DataVolatile &&other); // move constructor
        DataVolatile& operator=(const DataVolatile &other);
        
        // get methods
        int get_iter_num() const;
        double get_lower_glob() const;
        double get_upper_glob() const;  
        MIQP_Status get_status() const;
        int get_qp_iter() const;
        double get_run_time() const;
        Eigen::VectorXd get_x() const;
        double get_qp_solve_time() const;
        int get_qp_iter_avg() const; 
        std::vector<int> get_region_vec() const;

        // set methods
        void set_status(MIQP_Status status);
        void set_lower_glob(std::pair<double, bool> lower_glob_valid);
        void set_upper_glob(double upper_glob);
        void set_setup_time(double setup_time);
        void set_solve_time(double solve_time);
        void update_optimum(const std::unique_ptr<Node> &leaf);
        void increment_qp_iter(int n);
        void increment_iter_num();
        void increment_qp_solve_time(double t);
        void init_lower_glob_threads(int num_threads);
        void set_lower_glob_threads(int thread_number, double lower_glob);
        void reset_lower_glob_threads(int thread_number);

        // calculate
        void calc_qp_iter_avg();
        void calc_run_time();

        // reset statistics
        void reset(const Data& data);

    protected:
        bool first_run;
        int iter_num;
        double qp_solve_time;
        int qp_iter, qp_iter_avg;
        double lower_glob, upper_glob;
        Eigen::VectorXd x, x0;
        MIQP_Status status;
        double setup_time, solve_time, run_time;
        std::vector<int> region_vec;
        std::vector<double> lower_glob_threads;

    private:
        mutable std::mutex mtx;
};


} // end namespace MI_QPIPMPC


#endif