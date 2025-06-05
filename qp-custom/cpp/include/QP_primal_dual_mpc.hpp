#ifndef _QP_PRIMAL_DUAL_MPC_HPP_
#define _QP_PRIMAL_DUAL_MPC_HPP_

#include "Eigen/Sparse"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <cmath>
#include <map>
#include <chrono>
#include <utility>
#include "QP_consts.hpp"

#include "diag_flag.hpp"

namespace QP_IP_MPC
{

using QP::inf;
using QP::eps;

struct QP_settings
{
    double mu_term = 1e-3;
    double mu_feas = 1e0;
    double mu_max = 1e12;
    double mu_init = 1e7;
    unsigned int iter_max = 100;
    double gamma = 0.999;
    double t_ls = 0.9;
    double T_max = 0;
};

struct QP_results
{
    Eigen::VectorXd y;
    Eigen::VectorXd nu;
    Eigen::VectorXd lambda;
    Eigen::VectorXd s;
    double objective;
    bool converged;
    unsigned int num_iter;
    double sol_time;
    bool feas;
};

class QP_primal_dual_mpc
{
    private:
        // utilities
        void inline getTripletsForMatrix(const Eigen::Ref<const Eigen::MatrixXd> mat, std::vector<Eigen::Triplet<double>> &tripvec,
            int rowOffset, int colOffset);
        void inline getTripletsForMatrix(const Eigen::SparseMatrix<double> * mat_ptr, std::vector<Eigen::Triplet<double>> &tripvec,
            int rowOffset, int colOffset);

        bool checkIfBox(const Eigen::SparseMatrix<double> * M_ptr);
        bool checkIfDiag(const Eigen::SparseMatrix<double> * M_ptr);
        
    protected:
        // cost
        Eigen::SparseMatrix<double> P_i_nom;
        Eigen::VectorXd q_i_nom;
        std::map<int, Eigen::SparseMatrix<double>> P_i_vec;
        std::map<int, Eigen::VectorXd> q_i_vec;
        Eigen::SparseMatrix<double> P;
        Eigen::VectorXd q;
        double b;

        // equality constraints
        Eigen::SparseMatrix<double> C_i_nom;
        Eigen::SparseMatrix<double> D_i_nom;
        Eigen::VectorXd crhs_i_nom;
        std::map<int, Eigen::SparseMatrix<double>> C_i_vec;
        std::map<int, Eigen::SparseMatrix<double>> D_i_vec;
        std::map<int, Eigen::VectorXd> crhs_i_vec;
        Eigen::SparseMatrix<double> C;
        Eigen::VectorXd crhs;
        Eigen::MatrixXd C_i_nom_T_dense, D_i_nom_T_dense;
        std::map<int, Eigen::MatrixXd> C_i_vec_T_dense, D_i_vec_T_dense;
        Eigen::SparseMatrix<double> C_i_nom_T, D_i_nom_T;
        std::map<int, Eigen::SparseMatrix<double>> C_i_vec_T, D_i_vec_T;
        Eigen::SparseMatrix<double> C_T;

        // inequality constraints
        Eigen::SparseMatrix<double> G_i_nom;
        Eigen::VectorXd w_i_nom;
        std::map<int, Eigen::SparseMatrix<double>> G_i_vec;
        std::map<int, Eigen::VectorXd> w_i_vec;
        Eigen::SparseMatrix<double> G;
        Eigen::VectorXd w;
        Eigen::SparseMatrix<double> G_i_nom_T;
        std::map<int, Eigen::SparseMatrix<double>> G_i_vec_T;
        Eigen::SparseMatrix<double> G_T;

        // box constraint tracking
        bool isBox_i_nom;
        std::map<int, bool> isBox_i_vec;

        // diagonal cost tracking
        bool isDiag_i_nom;
        std::map<int, bool> isDiag_i_vec;
        std::map<int, bool> isDiag_phi_i_vec;

        // problem matrices
        std::map<int, Eigen::MatrixXd> CPC;
        std::map<int, Eigen::MatrixXd> DPD;
        std::map<int, Eigen::MatrixXd> DPC;
        std::map<int, Eigen::MatrixXd> L_i_i_vec;
        std::map<int, Eigen::MatrixXd> L_ip1_i_vec;
        std::map<int, Eigen::MatrixXd> L_phi_i_vec;
        std::map<int, Eigen::DiagonalMatrix<double, -1>> phi_i_inv_vec;
        Eigen::SparseMatrix<double> L_phi;
        std::map<int, Eigen::DiagonalMatrix<double, -1>> S_inv_Lambda_vec;
        Eigen::DiagonalMatrix<double, -1> S_inv, Lambda, S_inv_Lambda, S;
        
        // Cholesky decomposition solver
        Eigen::LLT<Eigen::MatrixXd> LLT_solver;

        // right hand side terms
        Eigen::VectorXd beta;
        Eigen::VectorXd r_C, r_E, r_I, r_S, r_d;
        Eigen::VectorXd nu_centering;
        Eigen::VectorXd phixb;

        // primal-dual vars
        Eigen::VectorXd y, nu, s, lambda;
        Eigen::VectorXd y0, nu0, s0, lambda0;

        // step direction
        Eigen::VectorXd del_y, del_nu, del_s, del_lambda;
        Eigen::DiagonalMatrix<double, -1> Del_s;

        // horizon length
        int n_horizon;

        // problem dimensions
        int n_y_tot, n_nu_tot, n_lambda_s_tot;
        std::map<int, int> n_y_vec;
        std::map<int, int> n_nu_vec;
        std::map<int, int> n_lambda_s_vec;

        // settings
        QP_settings settings;

        // flags
        bool P_updated, q_updated, C_updated, crhs_updated, G_updated, w_updated, b_updated;

        // internal methods
        void factorizeSystemMatrix();
        void generateRHS();
        void updateRHS();
        void computeFullStepDirectionFrom_del_nu();
        void solve_PHIx_eq_b(const Eigen::Ref<const Eigen::VectorXd> b, Eigen::Ref<Eigen::VectorXd> x);
        void solveCholeskyLinSys(const Eigen::SparseMatrix<double> * L_lowerTriangular_ptr, const Eigen::Ref<const Eigen::VectorXd> b, Eigen::Ref<Eigen::VectorXd> x);
        void solveCholeskyLinSys(const Eigen::Ref<const Eigen::MatrixXd> L_lowerTriangular, const Eigen::Ref<const Eigen::VectorXd> b, Eigen::Ref<Eigen::VectorXd> x);
        std::pair<double, bool> lineSearch(const Eigen::Ref<const Eigen::VectorXd> del_s, const Eigen::Ref<const Eigen::VectorXd> del_lambda);
        void computeProblemDimensions();
        double objective(const Eigen::Ref<const Eigen::VectorXd> x);
        double computeMu(const Eigen::Ref<const Eigen::VectorXd> s, const Eigen::Ref<const Eigen::VectorXd> lambda);
        bool make_S_inv_Lambda();
        void build_P();
        void build_q();
        void build_C();
        void build_crhs();
        void build_G();
        void build_w();
        void initializeWorkingMatrices();
        void solveYdelnu_eq_beta();

    public:
        // constructor
        QP_primal_dual_mpc();

        // destructor
        ~QP_primal_dual_mpc();

        // setup
        void setup(const Eigen::SparseMatrix<double> &P_i_nom, const Eigen::VectorXd &q_i_nom, 
            const Eigen::SparseMatrix<double> &C_i_nom, const Eigen::SparseMatrix<double> &D_i_nom,
            const Eigen::VectorXd &crhs_i_nom, const Eigen::SparseMatrix<double> &G_i_nom,
            const Eigen::VectorXd &w_i_nom, int n_horizon, double b=0);    
        void set_settings(const QP_settings &settings);

        // solve
        QP_results solve();

        // set methods
        void set_P_vec(const std::map<int, Eigen::SparseMatrix<double>> &P_vec);
        void set_q_vec(const std::map<int, Eigen::VectorXd> &q_vec);
        void set_C_vec(const std::map<int, Eigen::SparseMatrix<double>> &C_vec);
        void set_D_vec(const std::map<int, Eigen::SparseMatrix<double>> &D_vec);
        void set_crhs_vec(const std::map<int, Eigen::VectorXd> &crhs_vec);
        void set_G_vec(const std::map<int, Eigen::SparseMatrix<double>> &G_vec);
        void set_w_vec(const std::map<int, Eigen::VectorXd> &w_vec);
        void set_b(double b);
};



} // end namespace QP_IP_MPC

#endif