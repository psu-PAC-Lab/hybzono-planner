#ifndef _MPC_MIQP_HPP_
#define _MPC_MIQP_HPP_

#include "MPC_QP.hpp"
#include "MI_MPC.hpp"
#include "QP_primal_dual.hpp"
#include "ZonoCpp.hpp"
#include "Eigen/Sparse"
#include "Eigen/Dense"
#include <vector>
#include <map>
#include <iostream>
#include <string>

namespace HybZonoMPC
{

class MPC_MIQP : public MPC_QP
{
    public:
        
        // constructor
        MPC_MIQP() = default;

        // setup methods
        void setup_hybzono_Hrep_constraints(const ZonoCpp::HybZono &X_hz, 
            const std::map<int, Eigen::MatrixXd> &A_Hrep, const std::map<int, Eigen::VectorXd> &b_Hrep,
            const Eigen::SparseMatrix<double> &Pc,
            const Eigen::SparseMatrix<double> &Q_hz=Eigen::SparseMatrix<double>(), double slack_min_hz=-1e4, double slack_max_hz=1e4, 
            const Eigen::Ref<const Eigen::VectorXd> region_cost_vec=Eigen::VectorXd(),
            const Eigen::SparseMatrix<double> &Pp=Eigen::SparseMatrix<double>(),
            const std::map<int, Eigen::MatrixXd> &A_Hrep_pos=std::map<int, Eigen::MatrixXd>(), 
            const std::map<int, Eigen::VectorXd> &b_Hrep_pos=std::map<int, Eigen::VectorXd>(),
            const Eigen::Ref<const Eigen::MatrixXd> d_mat=Eigen::MatrixXd(), const std::vector<double> &d_x02region=std::vector<double>(),
            const std::map<int, Eigen::VectorXd> &disturbance_region=std::map<int, Eigen::VectorXd>());
        void update_hybzono_Hrep_constraints(const ZonoCpp::HybZono &X_hz, 
            const std::map<int, Eigen::MatrixXd> &A_Hrep, const std::map<int, Eigen::VectorXd> &b_Hrep,
            const Eigen::Ref<const Eigen::VectorXd> region_cost_vec=Eigen::VectorXd(),
            const std::map<int, Eigen::MatrixXd> &A_Hrep_pos=std::map<int, Eigen::MatrixXd>(), 
            const std::map<int, Eigen::VectorXd> &b_Hrep_pos=std::map<int, Eigen::VectorXd>(),
            const Eigen::Ref<const Eigen::MatrixXd> d_region2region=Eigen::MatrixXd(), const std::vector<double> &d_x02region=std::vector<double>(),
            const std::map<int, Eigen::VectorXd> &disturbance_region=std::map<int, Eigen::VectorXd>());
        
        // need to clean these up
        void setup_hybzono_constraints(const ZonoCpp::HybZono &X_hz, const Eigen::SparseMatrix<double> &Pp,
            const Eigen::SparseMatrix<double> &Q_hz);
        void setup_hybzono_constraints(const ZonoCpp::HybZono &X_hz, const Eigen::SparseMatrix<double> &Pp,
            const Eigen::SparseMatrix<double> &Q_hz, double slack_max_hz, double slack_min_hz);
        void setup_hybzono_constraints(const ZonoCpp::HybZono &X_hz, const Eigen::SparseMatrix<double> &Pp);
        void setup_hybzono_constraints(const ZonoCpp::HybZono &X_hz, const Eigen::SparseMatrix<double> &Pp,
            const Eigen::SparseMatrix<double> &Q_hz,
            const Eigen::Ref<const Eigen::VectorXd> region_cost_vec);
        void setup_hybzono_constraints(const ZonoCpp::HybZono &X_hz, const Eigen::SparseMatrix<double> &Pp,
            const Eigen::SparseMatrix<double> &Q_hz, double slack_max_hz, double slack_min_hz,
            const Eigen::Ref<const Eigen::VectorXd> region_cost_vec);
        void setup_hybzono_constraints(const ZonoCpp::HybZono &X_hz, const Eigen::SparseMatrix<double> &Pp,
            const Eigen::Ref<const Eigen::VectorXd> region_cost_vec);

        void set_distance_limit(double d_lim);
        void set_dist_between_regions(const Eigen::Ref<const Eigen::MatrixXd> d_mat);
        void set_terminal_hybzono_constraints(const ZonoCpp::HybZono &X_hz_term, const Eigen::Ref<const Eigen::VectorXd> term_region_cost_vec,
            const Eigen::SparseMatrix<double> &Pp_term, const Eigen::SparseMatrix<double> &Q_term_slack_cost);

        void set_region_dependent_disturbances(const std::map<int, Eigen::VectorXd> &disturbance_region);

        // settings
        void set_MI_settings(const MI_QPIPMPC::MI_QPIPMPC_Settings &miqp_settings);

        // get methods
        int get_terminal_region_selection();
        MI_QPIPMPC::Results get_miqp_results();

        // build controller
        void build_controller();

        // control method
        std::pair<Eigen::VectorXd, bool> control(const Eigen::Ref<const Eigen::VectorXd> x, 
            const std::vector<Eigen::VectorXd> &x_ref=std::vector<Eigen::VectorXd>(), 
            const Eigen::Ref<const Eigen::VectorXd> u=Eigen::VectorXd(),
            const std::vector<Eigen::SparseMatrix<double>> &A_dyn_vec=std::vector<Eigen::SparseMatrix<double>>(), 
            const std::vector<Eigen::SparseMatrix<double>> &B_dyn_vec=std::vector<Eigen::SparseMatrix<double>>());

        // verbosity
        void set_output_stream_fcn(void (*funcptr)(const std::string &str));
        std::string get_verbose_output();

    protected:

        // flags
        bool controller_built = false;

        // problem data
        MI_QPIPMPC::QPIPMPC_Data qpipmpc_data;

        // state map matrix
        Eigen::SparseMatrix<double> Pp; // position map matrix
        Eigen::SparseMatrix<double> Pc; // constraint map matrix

        // constraint violation cost
        Eigen::SparseMatrix<double> Q_hz;
        bool softened_hybzono_constraints = false;
        double sigma_min_hz = -1e4;
        double sigma_max_hz = 1e4;

        // constraint violation costsettings_ptr
        // arry of indices for binary variables
        std::vector<int> idx_i_nom;
        std::map<int, std::vector<int>> idx_i_vec;
        std::map<int, std::vector<int>> idx_binvar;

        // original problem
        std::map<int, int> n_y_vec_orig;
        Eigen::SparseMatrix<double> P_i_nom_orig, C_i_nom_orig, D_i_nom_orig, G_i_nom_orig;
        Eigen::VectorXd crhs_i_nom_orig, w_i_nom_orig, q_i_nom_orig;
        std::map<int, Eigen::SparseMatrix<double>> P_i_vec_orig, C_i_vec_orig, D_i_vec_orig, G_i_vec_orig;
        std::map<int, Eigen::VectorXd> crhs_i_vec_orig, w_i_vec_orig, q_i_vec_orig;

        // qp solvers: projection between hybzono regions
        QP::QP_primal_dual qp_solver_pt2region;
        QP::QP_primal_dual qp_solver_region2region;

        Eigen::SparseMatrix<double> H_pt2region, A_pt2region;
        Eigen::VectorXd f_pt2region, b_pt2region;
        Eigen::SparseMatrix<double> H_region2region, A_region2region;
        Eigen::VectorXd f_region2region, b_region2region;

        // distance limit
        double d_lim = QP::inf;

        // initial region for x0
        int region_0;
        bool region_0_valid = false;
        int n_regions;

        // reachability
        std::map<int, std::vector<int>> R_x02region;
        std::map<int, std::map<int, std::vector<int>>> R_region2region;
        Eigen::MatrixXd d_mat; // better choice of data type?
        bool point2region_configured = false;
        bool region2region_configured = false;
        bool region2region_dist_computed = false;
        std::vector<int> regions_local;
        std::vector<double> d_x02region;
        bool point2region_dist_computed = false;

        // reachability solver flag
        QP::QP_settings reachability_solver_settings;

        // mixed integer solver
        MI_QPIPMPC::MI_MPC miqp_solver;
        MI_QPIPMPC::MI_QPIPMPC_Settings miqp_settings;
        MI_QPIPMPC::Results miqp_results;

        // solution info
        double qp_time, mi_time;
        int qp_iter;

        // hybzono constraints
        ZonoCpp::HybZono X_hz;
        Eigen::VectorXd region_cost_vec;
        bool region_costs = false;
        int n_pos_dim;
        bool map_updated = false;
        
        // terminal constraint
        ZonoCpp::HybZono X_hz_term;
        bool terminal_hybzono_constraint = false;
        Eigen::VectorXd term_region_cost_vec;
        Eigen::SparseMatrix<double> Q_term_slack_cost;
        Eigen::SparseMatrix<double> Pp_term;
        Eigen::SparseMatrix<double> C_term, D_term, G_term, P_term;
        Eigen::VectorXd crhs_term, w_term, q_term;
        std::vector<int> idx_term;
        std::vector<int> idx_binvar_term;
        int terminal_region_selection = 1;

        // H-rep
        bool Hrep_defined = false;
        std::map<int, Eigen::MatrixXd> A_Hrep, A_Hrep_pos;
        std::map<int, Eigen::VectorXd> b_Hrep, b_Hrep_pos;

        // region dependent disturbances
        std::map<int, Eigen::VectorXd> d_region;
        bool region_dependent_disturbances = false;

        // verbose output
        std::string verbose_output;

        // warm start
        std::vector<std::vector<int>> region_vec_warm_start;
        int term_region_warm_start = 1;

        // get problem dimensions
        void get_problem_dimensions();

        // configure solver
        void configure_solver();

        // solve optimization problem
        bool solve_optimization_problem();

        // build optimization problem
        void call_parent_setup_methods();
        void make_inequality_constraints();
        void make_equality_constraints();
        void make_terminal_hybzono_constraints();
        void make_cost_function();
        void make_terminal_region_equality_constraints();
        void make_terminal_region_cost();
        void make_terminal_region_inequality_constraints();

        // reachability calculations
        void configure_point2region();
        void configure_region2region();
        double dist_point2region(const Eigen::Ref<const Eigen::VectorXd> x, int region);
        double dist_region2region(int region1, int region2);
        void get_x0_reachability();
        void get_region2region_reachability();

        // warm start
        void generate_warm_start_nodes();
        void apply_region0_to_warm_start_nodes();
};

} // namespace HybZonoMPC

#endif