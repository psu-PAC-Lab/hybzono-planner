#ifndef _MPC_MIQP_HREP_HPP_
#define _MPC_MIQP_HREP_HPP_

#include "MPC_MIQP.hpp"

namespace HybZonoMPC
{

class MPC_MIQP_Hrep : public MPC_MIQP
{

public:

    // constructor
    MPC_MIQP_Hrep();

    // setup methods
    void setup_Hrep_constraints(const std::map<int, Eigen::MatrixXd> &A_Hrep, const std::map<int, Eigen::VectorXd> &b_Hrep,
        const Eigen::SparseMatrix<double> &Pp, double big_M, const Eigen::SparseMatrix<double> &Q_hr);
    void setup_Hrep_constraints(const std::map<int, Eigen::MatrixXd> &A_Hrep, const std::map<int, Eigen::VectorXd> &b_Hrep,
        const Eigen::SparseMatrix<double> &Pp, double big_M);
    void setup_Hrep_constraints(const std::map<int, Eigen::MatrixXd> &A_Hrep, const std::map<int, Eigen::VectorXd> &b_Hrep,
        const Eigen::SparseMatrix<double> &Pp, double big_M, const Eigen::SparseMatrix<double> &Q_hr,
        const Eigen::Ref<const Eigen::VectorXd> region_cost_vec);
    void setup_Hrep_constraints(const std::map<int, Eigen::MatrixXd> &A_Hrep, const std::map<int, Eigen::VectorXd> &b_Hrep,
        const Eigen::SparseMatrix<double> &Pp, double big_M,
        const Eigen::Ref<const Eigen::VectorXd> region_cost_vec);

    void build_controller();

    // unsupported inherited methods
    void set_terminal_hybzono_constraints(const ZonoCpp::HybZono &X_hz_term, const Eigen::Ref<const Eigen::VectorXd> term_region_cost_vec,
            const Eigen::SparseMatrix<double> &Pp_term, const Eigen::SparseMatrix<double> &Q_term_slack_cost);
    int get_terminal_region_selection();

protected:

    // constraint violation cost
    Eigen::SparseMatrix<double> Q_hr;
    bool softened_hrep_constraints = false;
    double sigma_max_hr = 1e4;

    // big M
    double M;

    // number of H-rep constraints
    int n_cons;

    // flags
    bool Hrep_setup = false;

    // internal methods
    void get_problem_dimensions();
    void make_inequality_constraints();
    void make_equality_constraints();
    void make_cost_function();

}; // end class

} // end namespace

#endif // _MPC_MIQP_HREP_HPP_