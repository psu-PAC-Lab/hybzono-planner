#ifndef _MI_QPIPMPC_MI_MPC_HPP_
#define _MI_QPIPMPC_MI_MPC_HPP_

#include "Node.hpp"
#include "Data.hpp"
#include "LeavesQueue.hpp"
#include "DataVolatile.hpp"
#include "VerboseOutputs.hpp"
#include "QP_consts.hpp"
#include <chrono>
#include <utility>
#include <string>
#include <thread>
#include <future>
#include <numeric>
#include <sstream>
#include <memory>
#include <cmath>
#include <exception>
#include <map>

// appear to need mex.h included when compliling for Matlab from Windows
// do not need this when using Linux
#if defined __has_include
#  if __has_include ("mex.h")
#    include "mex.h"
#  endif
#endif

namespace MI_QPIPMPC
{

struct Results
{
    Eigen::VectorXd x;
    double upper_glob;
    double lower_glob;
    double run_time;
    MIQP_Status status;
    double qp_solve_time;
    int qp_iter_avg;
    int iter_num;
    std::vector<int> region_vec;

    // constructors
    Results() = default;
    Results(const Eigen::Ref<const Eigen::VectorXd> x, double upper_glob, double lower_glob,
        double run_time, MIQP_Status status, double qp_solve_time, int qp_iter_avg,
        int iter_num, std::vector<int> region_vec) : 
        x(x), upper_glob(upper_glob), lower_glob(lower_glob), run_time(run_time), status(status),
        qp_solve_time(qp_solve_time), qp_iter_avg(qp_iter_avg), iter_num(iter_num),
        region_vec(region_vec) {}
};

class MI_MPC
{
    public:
        // constructor
        MI_MPC() = default;
        MI_MPC(MI_MPC &&other); // move constructor

        // copy assignment
        MI_MPC& operator=(const MI_MPC &other);

        // setup
        void setup(const QPIPMPC_Data &qpipmpc_data, const MI_QPIPMPC_Settings &settings, 
            const QP_IP_MPC::QP_settings &qp_settings, 
            const std::map<int, std::map<int, std::vector<int>>> &R_region2region,
            const std::map<int, std::vector<int>> &R_x02region, int region_0);

        void setup(const QPIPMPC_Data &qpipmpc_data, const MI_QPIPMPC_Settings &settings, 
            const QP_IP_MPC::QP_settings &qp_settings, 
            const std::map<int, std::map<int, std::vector<int>>> &R_region2region,
            const std::map<int, std::vector<int>> &R_x02region, int region_0,
            const std::map<int, Eigen::MatrixXd> &A_Hrep,
            const std::map<int, Eigen::VectorXd> &b_Hrep);

        void setup(const QPIPMPC_Data &qpipmpc_data, const MI_QPIPMPC_Settings &settings, 
            const QP_IP_MPC::QP_settings &qp_settings, 
            const std::map<int, std::map<int, std::vector<int>>> &R_region2region,
            const std::map<int, std::vector<int>> &R_x02region, int region_0,
            const Eigen::SparseMatrix<double> &C_term,
            const Eigen::SparseMatrix<double> &D_term,
            const Eigen::SparseMatrix<double> &G_term,
            const Eigen::SparseMatrix<double> &P_term,
            const Eigen::Ref<const Eigen::VectorXd> crhs_term,
            const Eigen::Ref<const Eigen::VectorXd> w_term,
            const Eigen::Ref<const Eigen::VectorXd> q_term,
            const std::vector<int> &idx_term,
            const std::vector<int> &idx_binvar_term);

        void setup(const QPIPMPC_Data &qpipmpc_data, const MI_QPIPMPC_Settings &settings, 
            const QP_IP_MPC::QP_settings &qp_settings, 
            const std::map<int, std::map<int, std::vector<int>>> &R_region2region,
            const std::map<int, std::vector<int>> &R_x02region, int region_0,
            const std::map<int, Eigen::MatrixXd> &A_Hrep,
            const std::map<int, Eigen::VectorXd> &b_Hrep,
            const Eigen::SparseMatrix<double> &C_term,
            const Eigen::SparseMatrix<double> &D_term,
            const Eigen::SparseMatrix<double> &G_term,
            const Eigen::SparseMatrix<double> &P_term,
            const Eigen::Ref<const Eigen::VectorXd> crhs_term,
            const Eigen::Ref<const Eigen::VectorXd> w_term,
            const Eigen::Ref<const Eigen::VectorXd> q_term,
            const std::vector<int> &idx_term,
            const std::vector<int> &idx_binvar_term);

        // warm start
        void set_warm_start_regions(const std::vector<std::vector<int>> &region_vec_warm_start, int term_region_warm_start);

        // solve
        Results solve();

        // update methods
        void update_P_i_nom(const Eigen::SparseMatrix<double> &P_i_nom);
        void update_q_i_nom(const Eigen::Ref<const Eigen::VectorXd> q_i_nom);
        void update_P_i_vec(const std::map<int, Eigen::SparseMatrix<double>> &P_i_vec);
        void update_q_i_vec(const std::map<int, Eigen::VectorXd> &q_i_vec);
        void update_C_i_nom(const Eigen::SparseMatrix<double> &C_i_nom);
        void update_D_i_nom(const Eigen::SparseMatrix<double> &D_i_nom);
        void update_crhs_i_nom(const Eigen::Ref<const Eigen::VectorXd> crhs_i_nom);
        void update_C_i_vec(const std::map<int, Eigen::SparseMatrix<double>> &C_i_vec);
        void update_D_i_vec(const std::map<int, Eigen::SparseMatrix<double>> &D_i_vec);
        void update_crhs_i_vec(const std::map<int, Eigen::VectorXd> &crhs_i_vec);
        void update_G_i_nom(const Eigen::SparseMatrix<double> &G_i_nom);
        void update_w_i_nom(const Eigen::Ref<const Eigen::VectorXd> w_i_nom);
        void update_G_i_vec(const std::map<int, Eigen::SparseMatrix<double>> &G_i_vec);
        void update_w_i_vec(const std::map<int, Eigen::VectorXd> &w_i_vec);
        void update_b(double b);
        void update_idx_i_nom(const std::vector<int> &idx_i_nom);
        void update_idx_i_vec(const std::map<int, std::vector<int>> &idx_i_vec);
        void update_idx_state(const std::map<int, std::vector<int>> &idx_state);
        void update_idx_input(const std::map<int, std::vector<int>> &idx_input);
        void update_idx_binvar(const std::map<int, std::vector<int>> &idx_binvar);
        void update_idx_y(const std::map<int, std::vector<int>> &idx_y);
        void update_idx_eq(const std::map<int, std::vector<int>> &idx_eq);
        void update_idx_ineq(const std::map<int, std::vector<int>> &idx_ineq);
        void update_reachability(const std::map<int, std::vector<int>> &R_x02region,
            const std::map<int, std::map<int, std::vector<int>>> &R_region2region,
            int region_0);
        void update_Hrep(const std::map<int, Eigen::MatrixXd> &A_Hrep,
            const std::map<int, Eigen::VectorXd> &b_Hrep);
        void update_terminal_constraint(const Eigen::SparseMatrix<double> &C_term,
            const Eigen::SparseMatrix<double> &D_term,
            const Eigen::SparseMatrix<double> &G_term,
            const Eigen::SparseMatrix<double> &P_term,
            const Eigen::Ref<const Eigen::VectorXd> crhs_term,
            const Eigen::Ref<const Eigen::VectorXd> w_term,
            const Eigen::Ref<const Eigen::VectorXd> q_term,
            const std::vector<int> &idx_term,
            const std::vector<int> &idx_binvar_term);

        // get logging string when verbose
        std::string get_log() const;
        void set_output_stream_fcn(void (*funcptr)(const std::string &str));

    protected:

        // solve
        void main_solve_loop(int thread_number);

        // flags
        bool update_pending = true;

        // internal update method
        void update();

        // verbose output
        VerboseOutputs verbose_output;
        
        // problem data
        Data data;
        DataVolatile data_volatile;

        // leaves queue
        LeavesQueue leaves;

        // branch and bound
        void branch_and_bound(const std::unique_ptr<Node> &leaf);
        void branch_at_timestep(const std::unique_ptr<Node> &leaf, int kb);
        void branch_at_terminal_constraint(const std::unique_ptr<Node> &leaf);
        void branch_most_fractional_Hrep(const std::unique_ptr<Node> &leaf);

        // workspace methods
        void update_return_status();
        bool is_converged() const; 
        bool can_continue(int n_dead_threads, int n_threads) const; 

};

} // end namespace MI_QPIPMPC

#endif