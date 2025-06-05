#ifndef _MI_QPIPMPC_DATA_HPP_
#define _MI_QPIPMPC_DATA_HPP_

#include "QP_primal_dual_mpc.hpp"
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include <map>
#include <vector>
#include <exception>
#include <limits>

namespace MI_QPIPMPC
{

struct QPIPMPC_Data
{
    Eigen::SparseMatrix<double> P_i_nom;
    Eigen::VectorXd q_i_nom;
    double b=0;
    std::map<int, Eigen::SparseMatrix<double>> P_i_vec;
    std::map<int, Eigen::VectorXd> q_i_vec;
    Eigen::SparseMatrix<double> C_i_nom;
    Eigen::SparseMatrix<double> D_i_nom;
    Eigen::VectorXd crhs_i_nom;
    std::map<int, Eigen::SparseMatrix<double>> C_i_vec;
    std::map<int, Eigen::SparseMatrix<double>> D_i_vec;
    std::map<int, Eigen::VectorXd> crhs_i_vec;
    Eigen::SparseMatrix<double> G_i_nom;
    Eigen::VectorXd w_i_nom;
    std::map<int, Eigen::SparseMatrix<double>> G_i_vec;
    std::map<int, Eigen::VectorXd> w_i_vec;
    std::vector<int> idx_i_nom;
    std::map<int, std::vector<int>> idx_i_vec;
    std::map<int, std::vector<int>> idx_state;
    std::map<int, std::vector<int>> idx_input;
    std::map<int, std::vector<int>> idx_binvar;
    std::map<int, std::vector<int>> idx_y;
    std::map<int, std::vector<int>> idx_eq;
    std::map<int, std::vector<int>> idx_ineq;
    Eigen::MatrixXd Pc;
    int n_horizon;
}; // end struct QPIPMPC_Data

struct MI_QPIPMPC_Settings
{
    double eps_feas = 1e-3;
    int max_iter_bb = 1e4;    
    bool verbose = true;
    double conv_rel = 1e-2; // objective function convergence, relative to magnitude
    double conv_abs = 1e-1; // objective function convergence, absolute
    double T_max = 0; // max execution time
    unsigned int n_threads = 1; // number of threads
    double alpha_most_frac = 0.2; // when using most fracitional branching, create an integer feasible node if most fractional var < this
    double max_cost = std::numeric_limits<double>::max(); // max cost for the node
}; // end struct MI_QPIPMPC_Settings

enum MIQP_Status
{
    MIQP_NO_SOL,
    MIQP_SOLVED,
    MIQP_INFEASIBLE,
    MIQP_MAX_ITER_FEASIBLE,
    MIQP_MAX_ITER_INFEASIBLE,
    MIQP_MAX_COST
}; // end enum MI_Status

enum QP_Status
{
    QP_SOLVED,
    QP_NO_SOL,
    QP_INFEASIBLE,
    QP_SUBOPTIMAL
}; // end enum QP_Status

class Data
{
    public:
        
        // fields
        QPIPMPC_Data qpipmpc_data;
        QP_IP_MPC::QP_settings qp_settings;
        std::map<int, std::vector<int>> R_x02region;
        std::map<int, std::map<int, std::vector<int>>> R_region2region;
        int region_0;
        int n_horizon;
        int n_regions;
        int n;
        std::map<int, Eigen::MatrixXd> A_Hrep;
        std::map<int, Eigen::VectorXd> b_Hrep;
        MI_QPIPMPC_Settings settings;
        QP_IP_MPC::QP_primal_dual_mpc solver;
        bool Hrep_defined = false;

        bool terminal_constraint = false;
        Eigen::SparseMatrix<double> C_term, D_term, G_term, P_term;
        Eigen::VectorXd crhs_term, w_term, q_term;
        std::vector<int> idx_term;
        std::vector<int> idx_binvar_term;
        int n_term_regions = 0;

        bool warm_start = false;
        std::vector<std::vector<int>> region_vec_warm_start;
        int term_region_warm_start;

        std::map<int, std::vector<int>> idx_pos_vars;

        // constructors
        Data(const QPIPMPC_Data &qpipmpc_data, 
            const std::map<int, std::vector<int>> &R_x02region, 
            const std::map<int, std::map<int, std::vector<int>>> &R_region2region,
            int region_0, const MI_QPIPMPC_Settings &settings,
            const QP_IP_MPC::QP_settings &qp_settings,
            const std::map<int, Eigen::MatrixXd> &A_Hrep,
            const std::map<int, Eigen::VectorXd> &b_Hrep);
        Data(const QPIPMPC_Data &qpipmpc_data, 
            const std::map<int, std::vector<int>> &R_x02region, 
            const std::map<int, std::map<int, std::vector<int>>> &R_region2region,
            int region_0, const MI_QPIPMPC_Settings &settings,
            const QP_IP_MPC::QP_settings &qp_settings);
        Data(Data &&other); // move constructor
        Data(); // default constructor
        Data(const Data &other); // copy constructor

        // copy assignment
        Data& operator=(const Data &other);
        
        // setup methods
        void set_qp_settings(QP_IP_MPC::QP_settings &qp_settings);
        void set_Hrep(const std::map<int, Eigen::MatrixXd> &A_Hrep,
            const std::map<int, Eigen::VectorXd> &b_Hrep);
        void set_terminal_constraint(const Eigen::SparseMatrix<double> &C_term,
            const Eigen::SparseMatrix<double> &D_term,
            const Eigen::SparseMatrix<double> &G_term,
            const Eigen::SparseMatrix<double> &P_term,
            const Eigen::Ref<const Eigen::VectorXd> crhs_term,
            const Eigen::Ref<const Eigen::VectorXd> w_term,
            const Eigen::Ref<const Eigen::VectorXd> q_term,
            const std::vector<int> &idx_term,
            const std::vector<int> &idx_binvar_term);

        // get problem dimension
        void get_problem_dimension();
        
        // update/set methods
        void set_P_i_nom(const Eigen::SparseMatrix<double> &P_i_nom);
        void set_q_i_nom(const Eigen::Ref<const Eigen::VectorXd> q_i_nom);
        void set_b(double b);
        void set_P_i_vec(const std::map<int, Eigen::SparseMatrix<double>> &P_i_vec);
        void set_q_i_vec(const std::map<int, Eigen::VectorXd> &q_i_vec);
        void set_C_i_nom(const Eigen::SparseMatrix<double> &C_i_nom);
        void set_D_i_nom(const Eigen::SparseMatrix<double> &D_i_nom);
        void set_crhs_i_nom(const Eigen::Ref<const Eigen::VectorXd> crhs_i_nom);
        void set_C_i_vec(const std::map<int, Eigen::SparseMatrix<double>> &C_i_vec);
        void set_D_i_vec(const std::map<int, Eigen::SparseMatrix<double>> &D_i_vec);
        void set_crhs_i_vec(const std::map<int, Eigen::VectorXd> crhs_i_vec);
        void set_G_i_nom(const Eigen::SparseMatrix<double> &G_i_nom);
        void set_w_i_nom(const Eigen::Ref<const Eigen::VectorXd> w_i_nom);
        void set_G_i_vec(const std::map<int, Eigen::SparseMatrix<double>> &G_i_vec);
        void set_w_i_vec(const std::map<int, Eigen::VectorXd> &w_i_vec);
        void set_idx_i_nom(const std::vector<int> &idx_i_nom);
        void set_idx_i_vec(const std::map<int, std::vector<int>> &idx_i_vec);
        void set_idx_state(const std::map<int, std::vector<int>> &idx_state);
        void set_idx_input(const std::map<int, std::vector<int>> &idx_input);
        void set_idx_binvar(const std::map<int, std::vector<int>> &idx_binvar);
        void set_idx_y(const std::map<int, std::vector<int>> &idx_y);
        void set_idx_eq(const std::map<int, std::vector<int>> &idx_eq);
        void set_idx_ineq(const std::map<int, std::vector<int>> &idx_ineq);

        void update_reachability(const std::map<int, std::vector<int>> &R_x02region,
            const std::map<int, std::map<int, std::vector<int>>> &R_region2region,
            int region_0);

        void set_warm_start_regions(const std::vector<std::vector<int>> &region_vec_warm_start, int term_region_warm_start);

    protected:
        
        // terminal constraints
        void apply_terminal_constraints_to_problem_data();

    private:

        // check for positive variables
        void get_pos_vars();

}; // end class Data

} // end namespace MI_QPIPMPC

#endif