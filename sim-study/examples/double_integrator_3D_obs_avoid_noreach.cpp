// solvers
#include "MPC_MIQP.hpp"

// maps
#include "zono3D.hpp"

// other includes
#include "ZonoCpp.hpp"
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "vrep_calcs.hpp"
#include <iostream>
#include <chrono>
#include <map>
#include <vector>
#include <cmath>
#include <utility>
#include <fstream>
#include <string>
#include <memory>

using namespace HybZonoMPC;

// reference
std::vector<Eigen::VectorXd> get_reference(int n_horizon, const Eigen::Ref<const Eigen::VectorXd> xr)
{
    return std::vector<Eigen::VectorXd> (n_horizon, xr);
}

// map enum
enum Maps
{
    map_zono3D
};

int main()
{
    // comparison params
    const unsigned int n_trials = 5;
    const double T_max = 60;
    const unsigned int n_threads = 16;
    const bool verbose = false;
    const double conv_rel = 0.01;
    const double conv_abs = 0.1;

    // output solution times to file
    std::ofstream file, file_traj;
    file.open("../examples/data/solver_noreach_3D.txt");
    file_traj.open("../examples/data/solver_noreach_traj_3D.txt");

    // problem definition
    const double max_speed = 1;
    const double dT = 1;
    const std::array<int, 4> n_horizon_arr = {5, 10, 15, 20};
    
    // dynamics matrices
    Eigen::Matrix<double, 6, 6> A_dyn;
    Eigen::Matrix<double, 6, 3> B_dyn;
    A_dyn << 1, dT, 0, 0, 0, 0,
             0, 1, 0, 0, 0, 0,
             0, 0, 1, dT, 0, 0,
             0, 0, 0, 1, 0, 0,
             0, 0, 0, 0, 1, dT,
             0, 0, 0, 0, 0, 1;
    B_dyn << 0.5*dT*dT, 0, 0,
             dT, 0, 0,
             0, 0.5*dT*dT, 0,
             0, dT, 0,
             0, 0, 0.5*dT*dT,
             0, 0, dT;

    // cost matrices
    Eigen::Matrix<double, 6, 6> Q_cost, P_cost;
    Eigen::Matrix<double, 3, 3> R_cost;
    Q_cost = Eigen::Matrix<double, 6, 6>::Zero();
    P_cost = Eigen::Matrix<double, 6, 6>::Zero();
    R_cost = Eigen::Matrix<double, 3, 3>::Zero();
    Q_cost.diagonal() << 0.1, 0, 0.1, 0, 0.1, 0;
    R_cost.diagonal() << 10, 10, 10;
    P_cost.diagonal() << 10, 0, 10, 0, 10, 0;

    // conzono constraints
    
    // state constraints - only on velocity
    Eigen::Matrix<double, 3, 3> Gx = Eigen::Matrix<double, 3, 3>::Zero();
    Gx.diagonal() << max_speed, max_speed, max_speed;
    Eigen::Vector<double, 3> cx = Eigen::Vector<double, 3>::Zero();
    Eigen::Matrix<double, 0, 3> Ax = Eigen::Matrix<double, 0, 3>::Zero();
    Eigen::Vector<double, 0> bx = Eigen::Vector<double, 0>::Zero();
    ZonoCpp::ConZono Zc_x (Gx.sparseView(), cx, Ax.sparseView(), bx);
    Eigen::Vector<double, 6> x_max_hard = 1e4*Eigen::Vector<double, 6>::Ones();
    Eigen::Vector<double, 6> x_min_hard = -1e4*Eigen::Vector<double, 6>::Ones();
    Eigen::Matrix<double, 3, 6> Pp_x;
    Pp_x << 0, 1, 0, 0, 0, 0,
            0, 0, 0, 1, 0, 0,
            0, 0, 0, 0, 0, 1;
    
    // soft state constraints
    Eigen::Matrix<double, 3, 3> Qx_ineq_cost = 1e6*Eigen::Matrix<double, 3, 3>::Identity();
    Eigen::Vector<double, 3> sigma_max_x = 1e4*Eigen::Vector<double, 3>::Ones();

    // terminal state constraints - require at rest for persistent feasibility
    Eigen::Matrix<double, 3, 3> Gx_term = Eigen::Matrix<double, 3, 3>::Zero();
    Gx_term.diagonal() << 0, 0, 0;
    ZonoCpp::ConZono Zc_x_term (Gx_term.sparseView(), cx, Ax.sparseView(), bx);

    // soft terminal state constraints
    Eigen::Matrix<double, 3, 3> Qx_term_ineq_cost = 1e6*Eigen::Matrix<double, 3, 3>::Identity();
    Eigen::Vector<double, 3> sigma_max_x_term = 1e4*Eigen::Vector<double, 3>::Ones();

    // input constraints
    Eigen::Matrix<double, 3, 3> Gu;
    Gu.diagonal() << 1, 1, 1;
    Eigen::Vector<double, 3> cu = Eigen::Vector<double, 3>::Zero();
    Eigen::Matrix<double, 0, 3> Au = Eigen::Matrix<double, 0, 3>::Zero();
    Eigen::Vector<double, 0> bu = Eigen::Vector<double, 0>::Zero();
    ZonoCpp::ConZono Zc_u (Gu.sparseView(), cu, Au.sparseView(), bu);
    Eigen::Vector<double, 3> u_max_hard = 1e4*Eigen::Vector<double, 3>::Ones();
    Eigen::Vector<double, 3> u_min_hard = -1e4*Eigen::Vector<double, 3>::Ones();
    Eigen::Matrix<double, 3, 3> Pp_u = Eigen::Matrix<double, 3, 3>::Identity();

    // mpc settings
    MPC_Settings mpc_settings, mpc_settings_warmstart;
    mpc_settings.warm_start = false;
    mpc_settings_warmstart.warm_start = true;

    // mi settings
    MI_QPIPMPC::MI_QPIPMPC_Settings mi_settings_ours;
    mi_settings_ours.n_threads = n_threads;
    mi_settings_ours.T_max = T_max; 
    mi_settings_ours.verbose = verbose;
    mi_settings_ours.conv_abs = conv_abs;
    mi_settings_ours.conv_rel = conv_rel;
    mi_settings_ours.max_iter_bb = 50000;

    // declare mpc objects
    std::unique_ptr<MPC_MIQP> mpc_ours_hz;

    // vector of maps
    std::vector<Maps> map_vec = {map_zono3D};
   
    std::map<Maps, std::string> map_name_map = {
        {map_zono3D, "zono3D"}
    };  

    // declare params that vary with map
    ZonoCpp::HybZono X_hz;
    std::map<int, Eigen::MatrixXd> A_Hrep;
    std::map<int, Eigen::VectorXd> b_Hrep;
    Eigen::VectorXd pos_0, pos_ref;
    double T_sim;
    Eigen::VectorXd x0, xr;
    Eigen::VectorXd cost_vec;
    Eigen::MatrixXd d_mat;
    std::vector<Eigen::MatrixXd> V_obs, V_free;
    Eigen::MatrixXd V_lim;

    // loop through time steps
    for (const int & n_horizon : n_horizon_arr)
    {
        // print horizon to file / console
        file << "n_horizon: " << n_horizon << std::endl;
        std::cout << std::endl << "n_horizon: " << n_horizon << std::endl << std::endl;

        // update settings struct
        mpc_settings.n_horizon = n_horizon;
        mpc_settings_warmstart.n_horizon = n_horizon;

        // loop through maps
        for (auto const & map : map_vec)
        {
            // print map name to file / console
            file << map_name_map[map] << std::endl;
            std::cout << std::endl << map_name_map[map] << std::endl << std::endl;

            if (n_horizon == 15)
                file_traj << map_name_map[map] << std::endl;

            switch (map)
            {
                case map_zono3D:
                {
                    zono3D(X_hz, A_Hrep, b_Hrep, pos_0, pos_ref, T_sim, d_mat, V_obs, V_free, V_lim);
                    break;
                }
            }

            // get initial condition / reference / number of sim cycles
            int n_cycles = (int) (T_sim / dT);
            x0 = Eigen::VectorXd::Zero(6);
            x0(0) = pos_0(0);
            x0(2) = pos_0(1);
            x0(4) = pos_0(2);
            xr = Eigen::VectorXd::Zero(6);
            xr(0) = pos_ref(0);
            xr(2) = pos_ref(1);
            xr(4) = pos_ref(2);

            // H-rep constraint params
            int n_cons = 0; // init
            for (auto const &[key, val] : A_Hrep)
            {
                n_cons += val.rows();
            }
            Eigen::MatrixXd Q_hr = 1e6*Eigen::MatrixXd::Identity(n_cons, n_cons); // slack cost
            double big_M = 1e3;

            // hybzono constraint params
            Eigen::Matrix<double, 3, 3> Q_hz = 1e6*Eigen::Matrix<double, 3, 3>::Identity(); // slack cost
            
            // constraint mapping
            Eigen::Matrix<double, 3, 6> Pp;
            Pp << 1, 0, 0, 0, 0, 0,
                0, 0, 1, 0, 0, 0,
                0, 0, 0, 0, 1, 0;

            // loop through trials
            for (int i=0; i < n_trials; i++)
            {
                // create object
                mpc_ours_hz = std::make_unique<MPC_MIQP>();

                // dynamics and cost
                mpc_ours_hz->set_dyn_matrices(A_dyn.sparseView(), B_dyn.sparseView());
                mpc_ours_hz->set_stage_cost(Q_cost.sparseView(), R_cost.sparseView());
                mpc_ours_hz->set_terminal_cost(P_cost.sparseView());

                // state/input constraints
                mpc_ours_hz->set_state_inequality_constraints(Zc_x, Pp_x.sparseView(), x_min_hard, x_max_hard);
                mpc_ours_hz->set_terminal_state_inequality_constraints(Zc_x_term, Pp_x.sparseView(), x_min_hard, x_max_hard);
                mpc_ours_hz->set_input_inequality_constraints(Zc_u, Pp_u.sparseView(), u_min_hard, u_max_hard);
                mpc_ours_hz->set_state_inequality_constraints_slack_vars_cost(Qx_ineq_cost.sparseView(), sigma_max_x);
                mpc_ours_hz->set_terminal_state_inequality_constraints_slack_vars_cost(Qx_term_ineq_cost.sparseView(), sigma_max_x_term);

                // distance limit
                mpc_ours_hz->set_distance_limit(QP::inf); // DISABLE REACHABILITY CONSTRAINTS
                if (map == map_zono3D)
                    mpc_ours_hz->set_dist_between_regions(d_mat);

                // verbosity
                mpc_ours_hz->set_output_stream_fcn([](const std::string &str) { std::cout << str << std::endl; });

                // mpc settings
                mpc_ours_hz->set_mpc_settings(mpc_settings);
                
                // obs avoidance
                mpc_ours_hz->setup_hybzono_Hrep_constraints(X_hz, A_Hrep, b_Hrep, Pp.sparseView(), Q_hz.sparseView());   
                
                mpc_ours_hz->set_MI_settings(mi_settings_ours);

                // build controller
                mpc_ours_hz->build_controller();
          

                // sim setup

                // print trial number to console
                std::cout << std::endl << "Trial " << i << std::endl;

                // reference
                std::vector<Eigen::VectorXd> x_ref = get_reference(n_horizon, xr);

                // simulate
                Eigen::VectorXd x = x0;
                Eigen::VectorXd u;

                // declare controller output
                std::pair<Eigen::VectorXd, bool> control_output;
                bool feasible;

                // declare logging variables
                double tot_time = 0, max_time = 0, tot_time_warmstart = 0, max_time_warmstart = 0;

                MI_QPIPMPC::Results results;
                int mi_iter;
                double qp_sol_time;
                int tot_mi_iter = 0, tot_mi_iter_warmstart = 0;
                int max_mi_iter = 0, max_mi_iter_warmstart = 0;
                double tot_qp_sol_time = 0, tot_qp_sol_time_warmstart = 0;

                // log traj to file
                if (i==0 && n_horizon==15)
                    file_traj << x(0) << " " << x(2) << " " << x(4) << std::endl;

                // run MPC controller in loop and print control inputs to screen
                int n_cycles_completed = 0;
                for (int k = 0; k < n_cycles; k++)
                {
                    // start timer
                    std::chrono::time_point<std::chrono::high_resolution_clock> start_cycle = std::chrono::high_resolution_clock::now();
                    std::chrono::time_point<std::chrono::high_resolution_clock> end_cycle;

                    // get control input
                    control_output = mpc_ours_hz->control(x, x_ref);
                    end_cycle = std::chrono::high_resolution_clock::now();
                    results = mpc_ours_hz->get_miqp_results();

                    u = std::get<0>(control_output);
                    feasible = std::get<1>(control_output);
                    
                    // performance logging computations
                    mi_iter = results.iter_num;
                    qp_sol_time = results.qp_solve_time;

                    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_cycle - start_cycle);
                    double sol_time = duration.count()/1e6;
                    max_time = std::max(max_time, sol_time);
                    tot_time += sol_time;
                    if (k != 0)
                    {
                        max_time_warmstart = std::max(max_time_warmstart, sol_time);
                        tot_time_warmstart += sol_time;
                    }

                    max_mi_iter = std::max(max_mi_iter, mi_iter);
                    tot_mi_iter += mi_iter;
                    tot_qp_sol_time += qp_sol_time;
                    if (k != 0)
                    {
                        max_mi_iter_warmstart = std::max(max_mi_iter_warmstart, mi_iter);
                        tot_mi_iter_warmstart += mi_iter;
                        tot_qp_sol_time_warmstart += qp_sol_time;
                    }

                    // exception handling
                    if (!feasible || sol_time > T_max)
                        break;

                    // update state
                    x = A_dyn*x + B_dyn*u;

                    // print control input
                    std::cout << "u = " << u.transpose() << std::endl;

                    // log traj to file
                    if (i==0 && n_horizon==15)
                        file_traj << x(0) << " " << x(2) << " " << x(4) << std::endl;

                    // increment
                    n_cycles_completed++;
                }

                // logging
                double avg_time = tot_time / n_cycles_completed;
                double avg_mi_iter = ((double) tot_mi_iter) / n_cycles_completed;
                double avg_qp_sol_time = (tot_mi_iter == 0) ? 0 : tot_qp_sol_time / ((double) tot_mi_iter);

                double avg_time_warmstart = tot_time_warmstart / (n_cycles_completed-1);
                double avg_mi_iter_warmstart = ((double) tot_mi_iter_warmstart) / (n_cycles_completed-1);
                double avg_qp_sol_time_warmstart = (tot_mi_iter_warmstart == 0) ? 0 : tot_qp_sol_time_warmstart / ((double) tot_mi_iter_warmstart);

                file << avg_time << " " << max_time << " " << 
                    avg_mi_iter << " " << max_mi_iter << " " << 
                    avg_qp_sol_time << " " <<
                    avg_time_warmstart << " " << max_time_warmstart << " " <<
                    avg_mi_iter_warmstart << " " << max_mi_iter_warmstart << " " <<
                    avg_qp_sol_time_warmstart << std::endl;

                // free memory
                mpc_ours_hz.reset();
            }
        }

    } // end for n_horizon

    // close files
    file.close();
    file_traj.close();

    return 0;
}