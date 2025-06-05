#include "MPC_QP.hpp"
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "ZonoCpp.hpp"
#include <iostream>
#include <chrono>

using namespace ZonoCpp;

int main()
{
    // problem definition
    Eigen::Matrix<double, 2, 2> A_dyn;
    Eigen::Matrix<double, 2, 1> B_dyn;
    A_dyn << 1, 1, 0, 1;
    B_dyn << 0, 1;

    // cost matrices
    Eigen::Matrix<double, 2, 2> Q_cost, P_cost;
    Eigen::Matrix<double, 1, 1> R_cost;
    Q_cost << 1, 0, 0, 1;
    R_cost << 0.1;
    P_cost << 1, 0, 0, 1;

    // H-rep state and input constraints

    // // state constraint matrices: x1 in [-1, 1], x2 in [-1, 1]
    // Eigen::Matrix<double, 4, 2> Ax_ineq = Eigen::Matrix<double, 4, 2>::Zero();
    // Eigen::Vector<double, 4> bx_ineq;
    // Ax_ineq.block(0, 0, 2, 2) << Eigen::Matrix<double, 2, 2>::Identity();
    // Ax_ineq.block(2, 0, 2, 2) << -Eigen::Matrix<double, 2, 2>::Identity();
    // bx_ineq << 1, 1, 1, 1;

    // // terminal state constraint matrices, same as state constraint matrices here
    // Eigen::Matrix<double, 4, 2> Ax_term_ineq = Ax_ineq;
    // Eigen::Vector<double, 4> bx_term_ineq = bx_ineq;

    // // state constraint softening
    // Eigen::Matrix<double, 4, 4> Qx_ineq_cost = 1e6*Eigen::Matrix<double, 4, 4>::Identity();
    // Eigen::Matrix<double, 4, 4> Qx_term_ineq_cost = 1e6*Eigen::Matrix<double, 4, 4>::Identity();
    // Eigen::Matrix<double, 2, 2> Qu_ineq_cost = 1e6*Eigen::Matrix<double, 2, 2>::Identity();

    // Eigen::Vector<double, 4> sigma_max_x, sigma_max_x_term;
    // Eigen::Vector<double, 2> sigma_max_u;
    // sigma_max_x = 1e4*Eigen::Vector<double, 4>::Ones();
    // sigma_max_x_term = 1e4*Eigen::Vector<double, 4>::Ones();
    // sigma_max_u = 1e4*Eigen::Vector<double, 2>::Ones();

    // // input constraint matrices: u in [-1, 1]
    // Eigen::Matrix<double, 2, 1> Au_ineq;
    // Au_ineq << 1, -1;
    // Eigen::Vector<double, 2> bu_ineq;
    // bu_ineq << 1, 1;

    // conzono state and input constraints

    // state constraints: x1 in [-1, 1], x2 in [-1, 1]
    Eigen::Matrix<double, 2, 2> Gx = Eigen::Matrix<double, 2, 2>::Identity();
    Eigen::Vector<double, 2> cx = Eigen::Vector<double, 2>::Zero();
    Eigen::Matrix<double, 0, 2> Ax = Eigen::Matrix<double, 0, 2>::Zero();
    Eigen::Vector<double, 0> bx = Eigen::Vector<double, 0>::Zero();
    ConZono Zc_x (Gx.sparseView(), cx, Ax.sparseView(), bx);
    Eigen::Vector<double, 2> x_max_hard = 1e4*Eigen::Vector<double, 2>::Ones();
    Eigen::Vector<double, 2> x_min_hard = -1e4*Eigen::Vector<double, 2>::Ones();
    Eigen::Matrix<double, 2, 2> Pp_x = Eigen::Matrix<double, 2, 2>::Identity();

    // terminal state constraints
    ConZono Zc_x_term = Zc_x;
    Eigen::Vector<double, 2> x_term_max_hard = x_max_hard;
    Eigen::Vector<double, 2> x_term_min_hard = x_min_hard;
    Eigen::Matrix<double, 2, 2> Pp_x_term = Pp_x;

    // input constraints: u in [-1, 1]
    Eigen::Matrix<double, 1, 1> Gu = Eigen::Matrix<double, 1, 1>::Identity();
    Eigen::Vector<double, 1> cu = Eigen::Vector<double, 1>::Zero();
    Eigen::Matrix<double, 0, 1> Au = Eigen::Matrix<double, 0, 1>::Zero();
    Eigen::Vector<double, 0> bu = Eigen::Vector<double, 0>::Zero();
    ConZono Zc_u (Gu.sparseView(), cu, Au.sparseView(), bu);
    Eigen::Vector<double, 1> u_max_hard = 1e4*Eigen::Vector<double, 1>::Ones();
    Eigen::Vector<double, 1> u_min_hard = -1e4*Eigen::Vector<double, 1>::Ones();
    Eigen::Matrix<double, 1, 1> Pp_u = Eigen::Matrix<double, 1, 1>::Identity();

    // constraint softening
    Eigen::Matrix<double, 2, 2> Qx_ineq_cost = 1e6*Eigen::Matrix<double, 2, 2>::Identity();
    Eigen::Matrix<double, 2, 2> Qx_term_ineq_cost = Qx_ineq_cost;
    Eigen::Matrix<double, 1, 1> Qu_ineq_cost = 1e6*Eigen::Matrix<double, 1, 1>::Identity();
    Eigen::Vector<double, 2> sigma_max_x = 1e4*Eigen::Vector<double, 2>::Ones();
    Eigen::Vector<double, 2> sigma_max_x_term = sigma_max_x;
    Eigen::Vector<double, 1> sigma_max_u = 1e4*Eigen::Vector<double, 1>::Ones();

    // prediction horizon
    int n_horizon = 10;

    // loop time of controller
    double T_max = 0.1;

    // initial conditon
    Eigen::Vector<double, 2> x0;
    x0 << 0.9, 0;
    Eigen::Vector<double, 1> u0;
    u0 << 0;

    // create object
    HybZonoMPC::MPC_QP mpc;

    // define problem
    mpc.set_dyn_matrices(A_dyn.sparseView(), B_dyn.sparseView());
    mpc.set_stage_cost(Q_cost.sparseView(), R_cost.sparseView());
    mpc.set_terminal_cost(P_cost.sparseView());
    
    // // H-rep constraints
    // mpc.set_state_inequality_constraints(Ax_ineq.sparseView(), bx_ineq);
    // mpc.set_input_inequality_constraints(Au_ineq.sparseView(), bu_ineq);
    // mpc.set_terminal_state_inequality_constraints(Ax_term_ineq.sparseView(), bx_term_ineq);

    // conzono constraints
    mpc.set_state_inequality_constraints(Zc_x, Pp_x.sparseView(), x_min_hard, x_max_hard);
    mpc.set_input_inequality_constraints(Zc_u, Pp_u.sparseView(), u_min_hard, u_max_hard);
    mpc.set_terminal_state_inequality_constraints(Zc_x_term, Pp_x_term.sparseView(), x_term_min_hard, x_term_max_hard);

    mpc.set_state_inequality_constraints_slack_vars_cost(Qx_ineq_cost.sparseView(), sigma_max_x);
    mpc.set_input_inequality_constraints_slack_vars_cost(Qu_ineq_cost.sparseView(), sigma_max_u);
    mpc.set_terminal_state_inequality_constraints_slack_vars_cost(Qx_term_ineq_cost.sparseView(), sigma_max_x_term);
    
    
    HybZonoMPC::MPC_Settings mpc_settings;
    mpc_settings.warm_start = false;
    mpc_settings.n_horizon = n_horizon;
    mpc_settings.u1_control = true;
    mpc.set_mpc_settings(mpc_settings);

    // build controller
    mpc.build_controller();

    // simulate
    Eigen::VectorXd x = x0;
    Eigen::VectorXd u = u0;
    Eigen::VectorXd u1 = u0;

    // time execution
    double tot_time = 0;
    double max_time = 0;

    // run MPC controller in loop and print control inputs to screen
    const int n_cycles = 19;
    for (int i = 0; i < n_cycles; i++)
    {
        // update control
        u = u1;

        // start timer
        auto start_cycle = std::chrono::high_resolution_clock::now();

        // get control input
        auto [u1_sol, feasible] = mpc.control(x, std::vector<Eigen::VectorXd>(), u);
        u1 = u1_sol;

        // log execution time in second
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start_cycle);
        max_time = std::max(max_time, duration.count()/1e6);
        tot_time += duration.count()/1e6;

        // update state
        x = A_dyn*x + B_dyn*u;

        // print control input
        std::cout << "u = " << u.transpose() << std::endl;
    }

    // get total execution

    // print max execution time
    std::cout << "Max execution time: " << max_time << " sec" << std::endl;
    std::cout << "Average execution time: " << tot_time/n_cycles << " sec" << std::endl;

    return 0;
}