#ifndef _QP_PRIMAL_DUAL_HPP_
#define _QP_PRIMAL_DUAL_HPP_

#include "Eigen/Dense"
#include "Eigen/Sparse"
#include <vector>
#include <string>
#include <cmath>
#include <chrono>
#include <utility>
#include "QP_consts.hpp"

namespace QP
{

struct QP_settings
{
    double mu_term = 1e-3;
    double mu_feas = 1e0;
    double mu_max = 1e12;
    double mu_init = 1e7;
    double eps_feas = 1e-6;
    unsigned int iter_max = 100;
    double gamma = 0.999;
    double t_ls = 0.9;
    bool preprocessing_enable = true;
    double T_max = 0;
};

struct QP_results
{
    Eigen::VectorXd x;
    Eigen::VectorXd v;
    Eigen::VectorXd u;
    Eigen::VectorXd s;
    double objective;
    bool converged;
    unsigned int num_iter;
    double sol_time;
    bool feas;
};

class QP_primal_dual
{
    private:

        // utilities
        void inline getTripletsForMatrix(const Eigen::Ref<const Eigen::MatrixXd> mat, std::vector<Eigen::Triplet<double>> &tripvec,
            int rowOffset, int colOffset);
        void inline getTripletsForMatrix(const Eigen::SparseMatrix<double> * mat_ptr, std::vector<Eigen::Triplet<double>> &tripvec,
            int rowOffset, int colOffset);
        void inline getTripletsForMatrixDiagonal(const Eigen::Ref<const Eigen::VectorXd> d, std::vector<Eigen::Triplet<double>> &tripvec,
            int rowOffset, int colOffset);

    protected:

        // problem matrices / vectors
        Eigen::SparseMatrix<double> P;
        Eigen::VectorXd q;
        Eigen::SparseMatrix<double> A, A_T;
        Eigen::VectorXd b;
        Eigen::SparseMatrix<double> G, G_T;
        Eigen::VectorXd w;

        // linear system
        Eigen::SparseMatrix<double> M, M0, dM;
        std::vector<Eigen::Triplet<double>> tripvec_dM;
        Eigen::VectorXd r_C, r_E, r_I, r_S;
        Eigen::VectorXd bm;
        Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> LU_solver;
        Eigen::ComputationInfo LU_status = Eigen::ComputationInfo::Success;

        // preprocessing
        Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> QR_solver;
        Eigen::SparseMatrix<double> P_eq;

        // initial vars
        Eigen::VectorXd x0, v0, u0, s0;

        // working vars
        Eigen::VectorXd x, v, u, s;
        Eigen::VectorXd nu;
        Eigen::DiagonalMatrix<double, -1> S, Del_S;

        // duality measure
        double mu;

        // problem dimensions
        unsigned int n;
        unsigned int m_ineq;
        unsigned int m_eq;

        // flags
        bool equalityConstrained;
        bool A_updated;
        bool b_updated;
   
        // settings
        QP_settings settings;

        // helper methods
        double objective(const Eigen::Ref<const Eigen::VectorXd> x);
        void generateSystemMatrix();
        void generateRHS();
        void updateRHS();
        std::pair<double, bool> lineSearch(const Eigen::Ref<const Eigen::VectorXd> del_s, const Eigen::Ref<const Eigen::VectorXd> del_u);
        void getLinDepPermuteAndChopMatrix(const Eigen::SparseMatrix<double> * mat_in, Eigen::SparseMatrix<double> * mat_out);
        void getValidEqualityConstraints();
        void makeValid_A();
        void makeValid_b();
        double computeMu(const Eigen::Ref<const Eigen::VectorXd> s, const Eigen::Ref<const Eigen::VectorXd> u);
        void computeProblemDimensions();
        void initializeWorkingMatrices();

    public:

        // constructor
        QP_primal_dual();
        QP_primal_dual(const Eigen::SparseMatrix<double> &P,
            const Eigen::VectorXd &q,
            const Eigen::SparseMatrix<double> &A,
            const Eigen::VectorXd &b,
            const Eigen::SparseMatrix<double> &G,
            const Eigen::VectorXd &w);

        // destructor
        ~QP_primal_dual();

        // setup
        void setup(const Eigen::SparseMatrix<double> &P,
            const Eigen::VectorXd &q,
            const Eigen::SparseMatrix<double> &A,
            const Eigen::VectorXd &b,
            const Eigen::SparseMatrix<double> &G,
            const Eigen::VectorXd &w);

        // set settings
        void set_settings(const QP_settings *settings);

        // set methods
        void set_P(const Eigen::SparseMatrix<double> &P);
        void set_q(const Eigen::VectorXd &q);
        void set_A(const Eigen::SparseMatrix<double> &A);
        void set_b(const Eigen::VectorXd &b);
        void set_G(const Eigen::SparseMatrix<double> &G);
        void set_w(const Eigen::VectorXd &w);
        
        // solve method
        QP_results solve();

};

} // end namespace

#endif