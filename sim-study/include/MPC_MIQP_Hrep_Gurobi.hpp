#ifndef _MPC_MIQP_HREP_GUROBI_HPP_
#define _MPC_MIQP_HREP_GUROBI_HPP_

#include "MPC_MIQP_Hrep.hpp"
#include "MPC_MIQP_Gurobi.hpp"
#include "gurobi_c++.h"
#include <sstream>

namespace HybZonoMPC
{

class MPC_MIQP_Hrep_Gurobi : public MPC_MIQP_Hrep
{

public:

    // constructor
    MPC_MIQP_Hrep_Gurobi() = default;

    // destructor
    ~MPC_MIQP_Hrep_Gurobi();

    // setup methods
    void set_MI_settings(const Gurobi_Settings &settings);

    // overwritten methods
    void build_controller();

    // unsupported inherited methods
    void set_terminal_hybzono_constraints(const ZonoCpp::HybZono &X_hz_term, const Eigen::Ref<const Eigen::VectorXd> term_region_cost_vec,
            const Eigen::SparseMatrix<double> &Pp_term, const Eigen::SparseMatrix<double> &Q_term_slack_cost);
    int get_terminal_region_selection();

    // control methods (LTI)
    std::pair<Eigen::VectorXd, bool>  control(const Eigen::Ref<const Eigen::VectorXd> x, const std::vector<Eigen::VectorXd> &x_ref);
    std::pair<Eigen::VectorXd, bool>  control(const Eigen::Ref<const Eigen::VectorXd> x);

    // control methods (LTV)
    std::pair<Eigen::VectorXd, bool>  control(const Eigen::Ref<const Eigen::VectorXd> x, const std::vector<Eigen::VectorXd> &x_ref, 
        const std::vector<Eigen::SparseMatrix<double>> &A_dyn_vec, const std::vector<Eigen::SparseMatrix<double>> &B_dyn_vec);
    std::pair<Eigen::VectorXd, bool> control(const Eigen::Ref<const Eigen::VectorXd> x, const std::vector<Eigen::SparseMatrix<double>> &A_dyn_vec,
        const std::vector<Eigen::SparseMatrix<double>> &B_dyn_vec);

protected:

    // problem matrices
    Eigen::SparseMatrix<double> P, A, G;
    Eigen::VectorXd q, b, w;

    // make optimizaton problem
    void make_full_system_matrices();
    void make_optim_problem();

    // overwritten methods
    bool solve_optimization_problem();

private:

    // Gurobi state / settings
    GRBEnv env = GRBEnv();
    GRBModel * model_ptr = nullptr;
    Gurobi_Settings gurobi_settings;
    GRBVar * x_gur = nullptr;

    void set_gurobi_settings();

    // utilities
    void quad_cost(GRBModel * model_ptr, const GRBVar *x,
        const Eigen::SparseMatrix<double> * P, const Eigen::Ref<const Eigen::VectorXd> q);
    void lin_cons(GRBModel * model_ptr, const GRBVar *x, 
        const Eigen::Ref<const Eigen::MatrixXd> A_eq, const Eigen::Ref<const Eigen::VectorXd> b_eq,
        const Eigen::Ref<const Eigen::MatrixXd> A_ineq, const Eigen::Ref<const Eigen::VectorXd> b_ineq);


}; // end class

} // end namespace HybZonoMPC


#endif // _MPC_MIQP_HREP_GUROBI_HPP_