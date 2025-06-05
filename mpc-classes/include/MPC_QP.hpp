#ifndef _MPC_QP_HPP_
#define _MPC_QP_HPP_

#include "QP_primal_dual_mpc.hpp"
#include "QP_consts.hpp"
#include "ZonoCpp.hpp"
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <map>
#include <vector>
#include <utility>
#include <stdexcept>
#include <string>
#include <algorithm>
#include <sstream>

namespace HybZonoMPC
{

struct MPC_Settings
{
    bool warm_start = false;
    bool u1_control = false;
    int n_horizon = 15;
    bool ref_dep_term_constraint = false;
};

class MPC_QP
{
    public:

        // constructor
        MPC_QP() = default;

        // set methods
        void set_dyn_matrices(const Eigen::SparseMatrix<double> &A_dyn, const Eigen::SparseMatrix<double> &B_dyn);
        void set_dyn_matrices(const std::vector<Eigen::SparseMatrix<double>> &A_dyn_vec, const std::vector<Eigen::SparseMatrix<double>> &B_dyn_vec);
        void set_stage_cost(const Eigen::SparseMatrix<double> &Q_cost, const Eigen::SparseMatrix<double> &R_cost,
            const Eigen::SparseMatrix<double> &N_cost=Eigen::SparseMatrix<double>(), const Eigen::Ref<const Eigen::VectorXd> q_x_cost=Eigen::VectorXd(),
            const Eigen::Ref<const Eigen::VectorXd> q_u_cost=Eigen::VectorXd());
        void set_terminal_cost(const Eigen::SparseMatrix<double> &P_cost, const Eigen::Ref<const Eigen::VectorXd> q_xN_cost=Eigen::VectorXd());
        void set_state_inequality_constraints(const Eigen::SparseMatrix<double> &Ax_ineq, const Eigen::VectorXd &bx_ineq);
        void set_input_inequality_constraints(const Eigen::SparseMatrix<double> &Au_ineq, const Eigen::VectorXd &bu_ineq);
        void set_terminal_state_inequality_constraints(const Eigen::SparseMatrix<double> &Ax_term_ineq, 
            const Eigen::Ref<const Eigen::VectorXd> bx_term_ineq);
        void set_state_inequality_constraints(const ZonoCpp::ConZono &Zc_x, const Eigen::SparseMatrix<double> &Pp_x, 
            const Eigen::Ref<const Eigen::VectorXd> x_min_hard, const Eigen::Ref<const Eigen::VectorXd> x_max_hard);
        void set_input_inequality_constraints(const ZonoCpp::ConZono &Zc_u, const Eigen::SparseMatrix<double> &Pp_u,
            const Eigen::Ref<const Eigen::VectorXd> u_min_hard, const Eigen::Ref<const Eigen::VectorXd> u_max_hard);
        void set_terminal_state_inequality_constraints(const ZonoCpp::ConZono &Zc_x_term, const Eigen::SparseMatrix<double> &Pp_x_term,
            const Eigen::Ref<const Eigen::VectorXd> x_term_min_hard, const Eigen::Ref<const Eigen::VectorXd> x_term_max_hard);
        void set_state_inequality_constraints_slack_vars_cost(const Eigen::SparseMatrix<double> &Qx_constraint_cost, 
            Eigen::Ref<const Eigen::VectorXd> sigma_max_x);
        void set_input_inequality_constraints_slack_vars_cost(const Eigen::SparseMatrix<double> &Qu_constraint_cost,
            Eigen::Ref<const Eigen::VectorXd> sigma_max_u);
        void set_terminal_state_inequality_constraints_slack_vars_cost(const Eigen::SparseMatrix<double> &Qx_term_constraint_cost,
            Eigen::Ref<const Eigen::VectorXd> sigma_max_x_term);

        void set_mpc_settings(const MPC_Settings &mpc_settings);
        void set_qp_settings(const QP_IP_MPC::QP_settings &qp_settings);

        // build controller
        void build_controller();

        // control methods
        std::pair<Eigen::VectorXd, bool> control(const Eigen::Ref<const Eigen::VectorXd> x, 
            const std::vector<Eigen::VectorXd> &x_ref=std::vector<Eigen::VectorXd>(), 
            const Eigen::Ref<const Eigen::VectorXd> u=Eigen::VectorXd(),
            const std::vector<Eigen::SparseMatrix<double>> &A_dyn_vec=std::vector<Eigen::SparseMatrix<double>>(), 
            const std::vector<Eigen::SparseMatrix<double>> &B_dyn_vec=std::vector<Eigen::SparseMatrix<double>>());
        
        // get methods
        std::pair<std::vector<Eigen::VectorXd>, std::vector<Eigen::VectorXd>> get_trajectory();
        double get_objective();


    protected:

        // problem size
        int n_state, n_input;
        int n_y_tot, n_ineq_tot, n_eq_tot;
        std::map<int, int> n_y_vec, n_ineq_vec, n_eq_vec;
        std::map<int, std::vector<int>> idx_state, idx_input, idx_y, idx_ineq, idx_eq, idx_dyn_vec;
        std::vector<int> idx_dyn_nom;
        int start_idx_term_cons;

        // flags
        bool LTV = false;
        bool soft_state_ineq = false;
        bool soft_input_ineq = false;
        bool soft_term_state_ineq = false;
        bool terminal_cost_specified = false;
        bool state_ineq_specified = false;
        bool input_ineq_specified = false;
        bool terminal_state_ineq_specified = false;
        bool conzono_constraints = false;
        bool Hrep_constraints = false;
        bool ref_dep_term_constraint = false;

        // initial condition
        Eigen::VectorXd x0, u0;

        // reference
        std::vector<Eigen::VectorXd> x_ref;

        // dynamics matrices
        Eigen::SparseMatrix<double> A_dyn, B_dyn;
        std::vector<Eigen::SparseMatrix<double>> A_dyn_vec, B_dyn_vec;
        Eigen::SparseMatrix<double> Q_cost, R_cost, P_cost, N_cost;
        Eigen::VectorXd q_x_cost, q_u_cost, q_xN_cost;

        // softened constraint costs
        // J_k = (1/2)*(s_k^T) Q (s_k)
        // where Ax*x_k - sx_k <= bx
        Eigen::SparseMatrix<double> Qx_ineq_cost, Qx_term_ineq_cost, Qu_ineq_cost;

        // state constraints
        Eigen::SparseMatrix<double> Ax_ineq, Ax_term_ineq;
        Eigen::VectorXd bx_ineq, bx_term_ineq;
        ZonoCpp::ConZono Zc_x, Zc_x_term;
        Eigen::SparseMatrix<double> Pp_x, Pp_x_term;
        Eigen::VectorXd x_min_hard, x_max_hard, x_term_min_hard, x_term_max_hard;

        // input constraints
        Eigen::SparseMatrix<double> Au_ineq;
        Eigen::VectorXd bu_ineq;
        ZonoCpp::ConZono Zc_u;
        Eigen::SparseMatrix<double> Pp_u;
        Eigen::VectorXd u_min_hard, u_max_hard;

        // slack variable tbounds
        Eigen::VectorXd sigma_max_x, sigma_max_u, sigma_max_x_term;

        // solver
        QP_IP_MPC::QP_primal_dual_mpc solver_qp;
        Eigen::VectorXd solution;
        double objective;

        // solver problem definition
        Eigen::SparseMatrix<double> P_i_nom;
        Eigen::VectorXd q_i_nom;
        double b;
        std::map<int, Eigen::SparseMatrix<double>> P_i_vec;
        std::map<int, Eigen::VectorXd> q_i_vec;
        Eigen::SparseMatrix<double> C_i_nom, D_i_nom;
        Eigen::VectorXd crhs_i_nom;
        std::map<int, Eigen::SparseMatrix<double>> C_i_vec, D_i_vec;
        std::map<int, Eigen::VectorXd> crhs_i_vec;
        Eigen::SparseMatrix<double> G_i_nom;
        Eigen::VectorXd w_i_nom;
        std::map<int, Eigen::SparseMatrix<double>> G_i_vec;
        std::map<int, Eigen::VectorXd> w_i_vec;

        // settings
        QP_IP_MPC::QP_settings qp_settings;
        int n_horizon;  
        bool u1_control = false;
        bool warm_start = false;

        // solution vars
        std::vector<Eigen::VectorXd> x_vec, u_vec;

        // check problem validity
        std::pair<bool, std::string> check_problem_validity();

        // define optimization problem
        void make_inequality_constraints();
        void make_equality_constraints();
        void make_cost_function();
        void make_IC_constraint();
        void make_cost_gradient_from_reference();
        void make_const_cost_from_reference();
        void update_equality_constraints_LTV();
        void update_ref_dep_term_constraint();

        // get problem dimensions
        void get_problem_dimensions();

        // configure solver
        void configure_solver();

        // solve
        bool solve_optimization_problem();

        // get solution trajectory
        void get_state_input_trajectories();

        // utilities
        void get_triplets_offset(const Eigen::SparseMatrix<double> &mat, std::vector<Eigen::Triplet<double>> &triplets, 
            int i_offset, int j_offset);
        void get_triplets_offset(const Eigen::Ref<const Eigen::VectorXd> vec, std::vector<Eigen::Triplet<double>> &triplets, 
            int i_offset, int j_offset);
        void get_triplets_offset_identity(int n_dim, std::vector<Eigen::Triplet<double>> &triplets, int i_offset, int j_offset);
        void get_triplets_offset_neg_identity(int n_dim, std::vector<Eigen::Triplet<double>> &triplets, int i_offset, int j_offset);
        std::vector<int> get_indices(int offset, int range);

        // update sparse matrix
        template <typename T>
        void update_sparse_matrix(Eigen::SparseMatrix<T> &A, const Eigen::SparseMatrix<T> &A_sub_update, int i_offset, int j_offset)
        {
            // note: can probably do this more efficiently for cases where sparsity pattern can be guaranteed not to be changing

            // force all values in A over appropriate range to be zero
            for (int k=0; k<A.outerSize(); ++k)
            {
                for (typename Eigen::SparseMatrix<T>::InnerIterator it(A, k); it; ++it)
                {
                    if (it.row() >= i_offset && it.row() < i_offset + A_sub_update.rows() && 
                        it.col() >= j_offset && it.col() < j_offset + A_sub_update.cols())
                    {
                        it.valueRef() = 0;
                    }
                }
            }

            for (int k = 0; k < A_sub_update.outerSize(); ++k)
            {
                for (typename Eigen::SparseMatrix<T>::InnerIterator it(A_sub_update, k); it; ++it)
                {
                    try
                    {
                        // check if element exists in A
                        A.coeffRef(it.row() + i_offset, it.col() + j_offset) = it.value();
                    }
                    catch (std::exception &e)
                    {
                        // insert element in A
                        A.insert(it.row() + i_offset, it.col() + j_offset) = it.value();
                    }
                }
            }

            // Eigen::MatrixXd A_full = Eigen::MatrixXd(A);
            // A_full.block(i_offset, j_offset, A_sub_update.rows(), A_sub_update.cols()) = Eigen::MatrixXd(A_sub_update);
            // A = Eigen::SparseMatrix<T>(A_full.sparseView());
        }

        template <typename T>
        Eigen::SparseMatrix<T> select_cols(const Eigen::SparseMatrix<T> &A, const std::vector<int> &idx)
        {
            // enforce that idx is sorted
            if (!std::is_sorted(idx.begin(), idx.end()))
                throw std::invalid_argument("Indices must be sorted.");

            // enforce that A is compressed
            if (!A.isCompressed())
                throw std::invalid_argument("Matrix must be compressed.");

            // init output matrix
            Eigen::SparseMatrix<T> A_new (A.rows(), idx.size());
            std::vector<Eigen::Triplet<T>> tripvec_A_new;
            tripvec_A_new.reserve(A.nonZeros());

            // iterate through elements of matrix
            for (int k = 0; k < A.outerSize(); ++k)
            {
                for (typename Eigen::SparseMatrix<T>::InnerIterator it(A, k); it; ++it)
                {
                    // check if column index is in idx
                    auto it_col = std::find(idx.begin(), idx.end(), it.col());
                    if (it_col != idx.end())
                    {
                        // get column index in new matrix
                        int col = std::distance(idx.begin(), it_col);
                        tripvec_A_new.push_back(Eigen::Triplet<T>(it.row(), col, it.value()));
                    }
                }
            }

            // build matrix and return
            A_new.setFromTriplets(tripvec_A_new.begin(), tripvec_A_new.end());
            return A_new;
        }

        template <typename T>
        Eigen::SparseMatrix<T> select_rows_and_cols(const Eigen::SparseMatrix<T> &A, const std::vector<int> &idx_row, const std::vector<int> &idx_col)
        {
            // enforce that idx is sorted
            if (!std::is_sorted(idx_col.begin(), idx_col.end()))
                throw std::invalid_argument("Indices must be sorted.");

            if (!std::is_sorted(idx_row.begin(), idx_row.end()))
                throw std::invalid_argument("Indices must be sorted.");

            // enforce that A is compressed
            if (!A.isCompressed())
                throw std::invalid_argument("Matrix must be compressed.");

            // init output matrix
            Eigen::SparseMatrix<T> A_new (idx_row.size(), idx_col.size());
            std::vector<Eigen::Triplet<T>> tripvec_A_new;
            tripvec_A_new.reserve(A.nonZeros());

            // iterate through elements of matrix
            for (int k = 0; k < A.outerSize(); ++k)
            {
                for (typename Eigen::SparseMatrix<T>::InnerIterator it(A, k); it; ++it)
                {
                    // check if column and row index is in idx
                    auto it_col = std::find(idx_col.begin(), idx_col.end(), it.col());
                    auto it_row = std::find(idx_row.begin(), idx_row.end(), it.row());

                    if (it_col != idx_col.end() && it_row != idx_row.end())
                    {
                        // get column index in new matrix
                        int col = std::distance(idx_col.begin(), it_col);

                        // get row index in new matrix
                        int row = std::distance(idx_row.begin(), it_row);

                        tripvec_A_new.push_back(Eigen::Triplet<T>(row, col, it.value()));
                    }
                }
            }

            // build matrix and return
            A_new.setFromTriplets(tripvec_A_new.begin(), tripvec_A_new.end());
            return A_new;
        }
};

} // namespace HybZonoMPC

#endif