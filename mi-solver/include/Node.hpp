#ifndef _MI_QPIPMPC_NODE_HPP_
#define _MI_QPIPMPC_NODE_HPP_

#include "Data.hpp"
#include "QP_primal_dual_mpc.hpp"
#include "QP_consts.hpp"
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include <utility>
#include <map>
#include <vector>
#include <numeric>
#include <exception>
#include <algorithm>
#include <memory>

namespace MI_QPIPMPC
{

class Node
{
    public:
        
        // fields
        std::map<int, Eigen::SparseMatrix<double>> P_i_vec, C_i_vec, D_i_vec, G_i_vec;
        std::map<int, Eigen::VectorXd> q_i_vec, crhs_i_vec, w_i_vec;
        std::vector<std::pair<int,int>> idx_node2problem;

        QP_IP_MPC::QP_primal_dual_mpc solver;
        int depth;
        double obj;
        double lower;
        double frac_idx;
        int num_iter;
        double qp_solve_time;
        Eigen::VectorXd x;
        QP_Status status;
        QP_IP_MPC::QP_results qp_results;            
        double frac = 0;
        std::map<int, std::vector<int>> region_vec;
        std::map<int, std::vector<int>> region_vec_no_violation;
        bool cons_violated = true;

        std::vector<int> term_region_vec;
        int term_constraint_depth;
        bool terminal_constraint = false;
        double term_frac = 0;

        // constructor
        Node(const Data& data, const QP_IP_MPC::QP_primal_dual_mpc &solver,
            double lower, const std::map<int, std::vector<int>> &region_vec);
        Node() = default;

        Node(const Data& data, const QP_IP_MPC::QP_primal_dual_mpc &solver,
            double lower, const std::map<int, std::vector<int>> &region_vec,
            const std::vector<int> &term_region_vec);

        // sorting methods - returns true if this worse that n
        bool operator < (const Node &n) const;
        double operator - (const Node &n) const;

        // set methods
        void set_region_vec(const std::map<int, std::vector<int>> &region_vec);
        void set_frac(double frac);

        // get methods
        void get_depth(const Data& data);

        // reset solver
        void reset_solver(const Data& data);

        // solve
        void solve(const Data& data);

        // check for H-rep constraint violations at given time step
        std::pair<bool, int> check_Hrep_constraint_violation(int k_timestep, const Data& data);

        // get regions without constraint violation
        void get_regions_no_constraint_violation(const Data& data);
        int k_first_cons_violation;

        // get fractionality vector
        Eigen::VectorXd get_frac_vec(int k, const std::vector<int> &regions, const Data& data) const;

        // reachability logic
        bool init_reachability_constraints(const Data& data);
        void apply_reachability_constraints(const Data& data);
        void make_reachability_constraints_consistent(const Data& data, std::map<int, std::vector<int>>& region_vec);
        bool update_reachability_constraints(const Data& data, int kb, 
            const std::vector<int> &regions_kb);
        bool check_reachability_valid(const Data& data);

        // terminal constraint logic
        void update_terminal_constraint(const std::vector<int> &term_region_vec);
        void apply_terminal_constraint(const Data& data); 
        void get_terminal_constraint_depth();
        Eigen::VectorXd get_term_frac_vec(const std::vector<int> &regions, const Data& data) const;
        void set_term_frac(double term_frac);

        // get approximate cost of a region_vec selection
        double get_approx_cost(const Data& data, const std::map<int, std::vector<int>> &region_vec) const;

        // give priority
        void give_priority();
        void remove_priority();
        bool is_priority() const;

        // validity check
        bool is_valid() const;

    protected:
        bool priority = false;

        // add back zeroed variables to solution
        void add_zeroed_variables(const Data& data, Eigen::VectorXd &x); // cannot be Eigen::Ref b/c changes size

    private:
        bool valid = false;

        // utilities
        void remove_cols(Eigen::SparseMatrix<double> * M, const std::vector<int> &cols);
        void remove_rows(Eigen::SparseMatrix<double> * M, const std::vector<int> &rows);
        void remove_rows_cols(Eigen::SparseMatrix<double> * M, const std::vector<int> &rows, const std::vector<int> &cols);
        void remove_elements(Eigen::VectorXd * v, const std::vector<int> &elems);
        std::vector<int> get_zero_rows(const Eigen::SparseMatrix<double> * M);
        std::vector<int> get_indices(int offset, int range);
        void hor_cat_row_major(const Eigen::SparseMatrix<double> * A_ptr, const Eigen::SparseMatrix<double> * B_ptr, std::shared_ptr<Eigen::SparseMatrix<double, Eigen::RowMajor>> &AB_ptr);
};

} // end namespace MI_QPIPMPC

#endif