#include "MPC_MIQP.hpp"

using namespace HybZonoMPC;

// setup hybzono constraints
void MPC_MIQP::setup_hybzono_Hrep_constraints(const ZonoCpp::HybZono &X_hz, 
            const std::map<int, Eigen::MatrixXd> &A_Hrep, const std::map<int, Eigen::VectorXd> &b_Hrep,
            const Eigen::SparseMatrix<double> &Pc,
            const Eigen::SparseMatrix<double> &Q_hz, double slack_min_hz, double slack_max_hz, 
            const Eigen::Ref<const Eigen::VectorXd> region_cost_vec,
            const Eigen::SparseMatrix<double> &Pp,
            const std::map<int, Eigen::MatrixXd> &A_Hrep_pos, const std::map<int, Eigen::VectorXd> &b_Hrep_pos,
            const Eigen::Ref<const Eigen::MatrixXd> d_mat, const std::vector<double> &d_x02region,
            const std::map<int, Eigen::VectorXd> &disturbance_region)
{
    // local variables
    int region;

    // hybzono definition
    this->X_hz = X_hz;
    if (!this->X_hz.zero_one_generators)
        this->X_hz.convert_generator_range();

    // number of regions
    this->n_regions = X_hz.nGb;

    // H-rep
    this->A_Hrep = A_Hrep;
    this->b_Hrep = b_Hrep;
    this->Hrep_defined = true;

    // H-rep validity checking
    if (A_Hrep.size() != this->n_regions || b_Hrep.size() != this->n_regions)
    {
        throw std::invalid_argument("MPC_MIQP: number of regions in A_Hrep and b_Hrep do not match X_hz");
    }

    region = 1;
    for (auto const & [key, val] : A_Hrep)
    {
        if (key != region)
        {
            throw std::invalid_argument("MPC_MIQP: H-rep region indices must be contiguous starting from 1");
        }
        if (val.cols() != Pc.rows())
        {
            throw std::invalid_argument("MPC_MIQP: H-rep region state map matrix size does not match the hybzono state dimension");
        }
        if (val.rows() != b_Hrep.at(key).size())
        {
            throw std::invalid_argument("MPC_MIQP: H-rep region b vector size does not match the hybzono state dimension");
        }

        region++;
    }
    region = 1;
    for (auto const & [key, val] : b_Hrep)
    {
        if (key != region)
        {
            throw std::invalid_argument("MPC_MIQP: H-rep region indices must be contiguous starting from 1");
        }
        region++;
    }

    // state map to hybzono cons
    this->Pc = Pc;

    // constraint softening
    if (Q_hz.rows() > 0)
    {
        this->Q_hz = Q_hz;
        this->softened_hybzono_constraints = true;
        this->sigma_max_hz = slack_max_hz;
        this->sigma_min_hz = slack_min_hz;

        if (Q_hz.rows() != X_hz.n || Q_hz.cols() != X_hz.n)
        {
            throw std::invalid_argument("MPC_MIQP: Q_hz dimension must agree with hybrid zonotope dimension");
        }
    }

    // region cost
    if (region_cost_vec.size() > 0)
    {
        this->region_cost_vec = region_cost_vec;
        this->region_costs = true;
    }

    // check if position map defined, if not assume equal to hybzono cons
    if (Pp.rows() > 0)
    {
        // check H-rep validity, require this to be defined
        if (A_Hrep_pos.size() != this->n_regions || b_Hrep_pos.size() != this->n_regions)
        {
            throw std::invalid_argument("MPC_MIQP: number of regions in A_Hrep_pos and b_Hrep_pos do not match X_hz");
        }

        region = 1;
        for (auto const & [key, val] : A_Hrep_pos)
        {
            if (key != region)
            {
                throw std::invalid_argument("MPC_MIQP: position H-rep region indices must be contiguous starting from 1");
            }
            if (val.cols() != Pp.rows())
            {
                throw std::invalid_argument("MPC_MIQP: position H-rep region state map matrix size does not match the position state dimension");
            }
            if (val.rows() != b_Hrep_pos.at(key).size())
            {
                throw std::invalid_argument("MPC_MIQP: position H-rep region b vector size does not match the position state dimension");
            }

            region++;
        }
        region = 1;
        for (auto const & [key, val] : b_Hrep_pos)
        {
            if (key != region)
            {
                throw std::invalid_argument("MPC_MIQP: position H-rep region indices must be contiguous starting from 1");
            }
            region++;
        }

        this->Pp = Pp;
        this->A_Hrep_pos = A_Hrep_pos;
        this->b_Hrep_pos = b_Hrep_pos;
    }
    else
    {
        this->Pp = Pc;
        this->A_Hrep_pos = A_Hrep;
        this->b_Hrep_pos = b_Hrep;
    }

    // dimension of position space
    this->n_pos_dim = this->Pp.rows();

    // region to region distance
    this->d_mat = d_mat; // will recompute if below is false
    this->region2region_dist_computed = (d_mat.rows() == n_regions && d_mat.cols() == n_regions); 

    // point to regon reachability
    if (d_x02region.size() == n_regions)
    {
        this->d_x02region = d_x02region;
        point2region_dist_computed = true;
    }

    // region-dependent disturbances
    set_region_dependent_disturbances(disturbance_region);

    // set flags
    this->map_updated = true;
}

void MPC_MIQP::update_hybzono_Hrep_constraints(const ZonoCpp::HybZono &X_hz, 
            const std::map<int, Eigen::MatrixXd> &A_Hrep, const std::map<int, Eigen::VectorXd> &b_Hrep,
            const Eigen::Ref<const Eigen::VectorXd> region_cost_vec,
            const std::map<int, Eigen::MatrixXd> &A_Hrep_pos, const std::map<int, Eigen::VectorXd> &b_Hrep_pos,
            const Eigen::Ref<const Eigen::MatrixXd> d_region2region, const std::vector<double> &d_x02region,
            const std::map<int, Eigen::VectorXd> &disturbance_region)
{
    // call setup method
    if (A_Hrep_pos.size() == X_hz.nGb && b_Hrep_pos.size() == X_hz.nGb) // check if position H-rep is valid and update
    {
        setup_hybzono_Hrep_constraints(X_hz, A_Hrep, b_Hrep, this->Pc, this->Q_hz, this->sigma_min_hz, this->sigma_max_hz, 
            region_cost_vec, this->Pp, A_Hrep_pos, b_Hrep_pos, d_region2region, d_x02region, disturbance_region);
    }
    else
    {
        setup_hybzono_Hrep_constraints(X_hz, A_Hrep, b_Hrep, this->Pc, this->Q_hz, this->sigma_min_hz, this->sigma_max_hz, 
            region_cost_vec, this->Pp, A_Hrep, b_Hrep, d_region2region, d_x02region, disturbance_region);
    }

    // rerun the needed setup methods if the controller is already built
    if (controller_built)
    {
        // get reachability matrices
        region_0_valid = false;
        get_x0_reachability();
        get_region2region_reachability();

        // create matrices for optimization problem
        make_inequality_constraints();
        make_equality_constraints();
        make_cost_function();

        // terminal constraints
        if (terminal_hybzono_constraint)
        {
            make_terminal_region_cost();
            make_terminal_region_equality_constraints();
            make_terminal_region_inequality_constraints();
        }

        // get state dimensions at each time step
        get_problem_dimensions();

        // initial solver setup and configuration
        configure_solver();
    }

}

void MPC_MIQP::setup_hybzono_constraints(const ZonoCpp::HybZono &X_hz, const Eigen::SparseMatrix<double> &Pp,
    const Eigen::SparseMatrix<double> &Q_hz)
{
    setup_hybzono_constraints(X_hz, Pp);

    // hybzono constraint violation cost
    this->Q_hz = Q_hz;
    this->softened_hybzono_constraints = true;
}

void MPC_MIQP::setup_hybzono_constraints(const ZonoCpp::HybZono &X_hz, const Eigen::SparseMatrix<double> &Pp,
    const Eigen::SparseMatrix<double> &Q_hz, double slack_max_hz, double slack_min_hz)
{
    setup_hybzono_constraints(X_hz, Pp);

    // hybzono constraint violation cost
    this->Q_hz = Q_hz;
    this->softened_hybzono_constraints = true;
    this->sigma_max_hz = slack_max_hz;
    this->sigma_min_hz = slack_min_hz;
}

void MPC_MIQP::setup_hybzono_constraints(const ZonoCpp::HybZono &X_hz, const Eigen::SparseMatrix<double> &Pp,
    const Eigen::SparseMatrix<double> &Q_hz,
    const Eigen::Ref<const Eigen::VectorXd> region_cost_vec)
{
    setup_hybzono_constraints(X_hz, Pp);

    // hybzono constraint violation cost
    this->Q_hz = Q_hz;
    this->softened_hybzono_constraints = true;

    // region cost
    this->region_cost_vec = region_cost_vec;
    this->region_costs = true;
}

void MPC_MIQP::setup_hybzono_constraints(const ZonoCpp::HybZono &X_hz, const Eigen::SparseMatrix<double> &Pp,
    const Eigen::SparseMatrix<double> &Q_hz, double slack_max_hz, double slack_min_hz,
    const Eigen::Ref<const Eigen::VectorXd> region_cost_vec)
{
    setup_hybzono_constraints(X_hz, Pp);

    // hybzono constraint violation cost
    this->Q_hz = Q_hz;
    this->softened_hybzono_constraints = true;
    this->sigma_max_hz = slack_max_hz;
    this->sigma_min_hz = slack_min_hz;

    // region cost
    this->region_cost_vec = region_cost_vec;
    this->region_costs = true;
}

void MPC_MIQP::setup_hybzono_constraints(const ZonoCpp::HybZono &X_hz, const Eigen::SparseMatrix<double> &Pp, 
    const Eigen::Ref<const Eigen::VectorXd> region_cost_vec)
{
    setup_hybzono_constraints(X_hz, Pp);

    // region cost
    this->region_cost_vec = region_cost_vec;
    this->region_costs = true;
}

void MPC_MIQP::setup_hybzono_constraints(const ZonoCpp::HybZono &X_hz, const Eigen::SparseMatrix<double> &Pp)
{
    // hybzono definition
    this->X_hz = X_hz;
    if (!this->X_hz.zero_one_generators)
        this->X_hz.convert_generator_range();

    // number of regions
    this->n_regions = X_hz.nGb;

    // dimension of constraint space
    this->n_pos_dim = X_hz.n;

    // state map to hybzono domain
    this->Pp = Pp;

    // H-rep
    this->Hrep_defined = false;
}

// set distance limit
void MPC_MIQP::set_distance_limit(double d_lim)
{
    this->d_lim = d_lim;
}

// region dependent disturbances
void MPC_MIQP::set_region_dependent_disturbances(const std::map<int, Eigen::VectorXd> &disturbance_region)
{
    this->d_region = disturbance_region;
    region_dependent_disturbances = disturbance_region.size() > 0;
}

// terminal region constraints
void MPC_MIQP::set_terminal_hybzono_constraints(const ZonoCpp::HybZono &X_hz_term, const Eigen::Ref<const Eigen::VectorXd> term_region_cost_vec,
    const Eigen::SparseMatrix<double> &Pp_term, const Eigen::SparseMatrix<double> &Q_term_slack_cost)
{
    // hybzono definition
    this->X_hz_term = X_hz_term;
    if (!this->X_hz_term.zero_one_generators)
        this->X_hz_term.convert_generator_range();

    // state map to hybzono domain
    this->Pp_term = Pp_term;

    // terminal region cost
    this->term_region_cost_vec = term_region_cost_vec;

    // terminal hybzono constraint violation cost
    this->Q_term_slack_cost = Q_term_slack_cost;

    this->terminal_hybzono_constraint = true;

    // validity checking
    if (term_region_cost_vec.size() != X_hz_term.nGb)
    {
        throw std::invalid_argument("MPC_MIQP: terminal region cost vector size does not match the number of terminal regions");
    }
    if (Pp_term.rows() != Q_term_slack_cost.rows())
    {
        throw std::invalid_argument("MPC_MIQP: terminal region state map and constraint violation cost matrix sizes do not match");
    }
    if (Pp_term.cols() != n_state)
    {
        throw std::invalid_argument("MPC_MIQP: terminal region state map matrix size does not match the state dimension");
    }
}

// set distance between regions
void MPC_MIQP::set_dist_between_regions(const Eigen::Ref<const Eigen::MatrixXd> d_mat)
{
    this->d_mat = d_mat;
    region2region_dist_computed = true;
}

// update settings
void MPC_MIQP::set_MI_settings(const MI_QPIPMPC::MI_QPIPMPC_Settings &miqp_settings)
{
    this->miqp_settings = miqp_settings;
}

// build controller
void MPC_MIQP::build_controller()
{
    // check for constraint compatability
    if (conzono_constraints && Hrep_constraints)
    {
        throw std::invalid_argument("Mixing H-rep and conzono constraints is not supported");
        return;
    }
    if (!conzono_constraints && !Hrep_constraints)
    {
        conzono_constraints = true;
    }
    if (ref_dep_term_constraint)
    {
        if (!conzono_constraints)
        {
            throw std::invalid_argument("Reference dependent terminal constraints require conzono constraints");
            return;
        }
        ZonoCpp::ConZono * Zc_x_term_ptr = terminal_state_ineq_specified ? &Zc_x_term : &Zc_x;
        if (Zc_x_term_ptr->n != this->n_state)
        {
            throw std::invalid_argument("Terminal state inequality constraints must have same dimension as state when using reference dependent terminal constraints");
            return;
        }
    }

    // if constraints aren't specified, build them up as +/- inf
    if (Hrep_constraints)
    {
        if (!state_ineq_specified)
        {
            Ax_ineq = Eigen::SparseMatrix<double>(1, n_state);
            bx_ineq = QP::inf*Eigen::VectorXd::Ones(1);
        }
        if (!terminal_state_ineq_specified)
        {
            Ax_term_ineq = Eigen::SparseMatrix<double>(1, n_state);
            bx_term_ineq = QP::inf*Eigen::VectorXd::Ones(1);
        }
        if (!input_ineq_specified)
        {
            Au_ineq = Eigen::SparseMatrix<double>(1, n_input);
            bu_ineq = QP::inf*Eigen::VectorXd::Ones(1);
        }
    }

    // make sure that H-rep is defined
    if (!Hrep_defined)
    {
        throw std::invalid_argument("Not including H-rep representation of non-convex constraints is not currently supported");
        return;
    }

    // initialize x0 and x_ref
    if (x0.size() != n_state)
        x0 = Eigen::VectorXd::Zero(n_state);
    if (x_ref.size() != n_horizon)
        x_ref = std::vector<Eigen::VectorXd>(n_horizon, Eigen::VectorXd::Zero(n_state));
    if (u0.size() != n_input)
        u0 = Eigen::VectorXd::Zero(n_input);

    // check problem validity
    auto [valid, msg] = check_problem_validity();
    if (!valid)
    {
        throw std::invalid_argument("Problem not valid: " + msg);
        return;
    }

    // get reachability matrices
    get_x0_reachability();
    get_region2region_reachability();

    // create matrices for optimization problem
    call_parent_setup_methods();
    make_inequality_constraints();
    make_equality_constraints();
    make_cost_function();

    // terminal constraints
    if (terminal_hybzono_constraint)
    {
        make_terminal_region_cost();
        make_terminal_region_equality_constraints();
        make_terminal_region_inequality_constraints();
    }

    // get state dimensions at each time step
    get_problem_dimensions();

    // initial solver setup and configuration
    configure_solver();

    // set flag
    controller_built = true;
}

// control methods (LTI)
std::pair<Eigen::VectorXd, bool> MPC_MIQP::control(const Eigen::Ref<const Eigen::VectorXd> x, 
            const std::vector<Eigen::VectorXd> &x_ref, const Eigen::Ref<const Eigen::VectorXd> u,
            const std::vector<Eigen::SparseMatrix<double>> &A_dyn_vec, 
            const std::vector<Eigen::SparseMatrix<double>> &B_dyn_vec)
{
    // set x0
    this->x0 = x;

    // handle optional inputs
    if (x_ref.size() == 0)
        this->x_ref = std::vector<Eigen::VectorXd>(n_horizon, Eigen::VectorXd::Zero(n_state));
    else if (x_ref.size() == 1)
        this->x_ref = std::vector<Eigen::VectorXd>(n_horizon, x_ref[0]);
    else if (x_ref.size() == n_horizon)
        this->x_ref = x_ref;
    else
        throw std::invalid_argument("Invalid reference trajectory length");

    if (u.size() == n_input && u1_control)
        this->u0 = u;
    else if (u1_control)
        throw std::invalid_argument("Require valid control input be specified");
    else if (u.size() != 0)
        throw std::invalid_argument("Invalid control input size");

    if (A_dyn_vec.size() == n_horizon && B_dyn_vec.size() == n_horizon)
    {
        // set LTV flag
        this->LTV = true;

        // updated dynamics matrices
        this->A_dyn_vec = A_dyn_vec;
        this->B_dyn_vec = B_dyn_vec;
        this->A_dyn = *(A_dyn_vec.begin());
        this->B_dyn = *(B_dyn_vec.begin());   
    }

    // reachability update
    get_x0_reachability();
    region_0_valid = true;

    // constraint updates
    make_IC_constraint();
    if (ref_dep_term_constraint)
        update_ref_dep_term_constraint();

    // gradient update
    make_cost_gradient_from_reference();
    make_const_cost_from_reference();

    // update solver
    miqp_solver.update_crhs_i_vec(crhs_i_vec);
    miqp_solver.update_q_i_vec(q_i_vec);
    miqp_solver.update_b(b);
    miqp_solver.update_reachability(R_x02region, R_region2region, region_0);

    // LTV update
    if (LTV)
    {
        update_equality_constraints_LTV();
        miqp_solver.update_C_i_vec(C_i_vec);
    }

    // update terminal region constraints
    if (terminal_hybzono_constraint)
    {
        make_terminal_region_equality_constraints();
        make_terminal_region_cost();
        make_terminal_region_inequality_constraints();
        get_problem_dimensions();
        miqp_solver.update_terminal_constraint(C_term, D_term, G_term, P_term, crhs_term, w_term, q_term, idx_term, idx_binvar_term);

        // need to update problem size, this will have been done twice for local map case
        miqp_solver.update_idx_y(idx_y);
        miqp_solver.update_idx_ineq(idx_ineq);
        miqp_solver.update_idx_eq(idx_eq);
    }

    // warm start
    if (warm_start && !map_updated)
    {
        apply_region0_to_warm_start_nodes();
        miqp_solver.set_warm_start_regions(region_vec_warm_start, term_region_warm_start);
    }

    // solve
    bool feasible = solve_optimization_problem();

    // reset flags
    map_updated = false;

    // return control input for time zero (or 1 if using u1 control)
    if (u1_control)
        return std::make_pair(u_vec[1], feasible);
    else
        return std::make_pair(u_vec[0], feasible);
}

// get methods
int MPC_MIQP::get_terminal_region_selection()
{
    return terminal_region_selection;
}

MI_QPIPMPC::Results MPC_MIQP::get_miqp_results()
{
    return miqp_results;
}

void MPC_MIQP::set_output_stream_fcn(void (*funcptr)(const std::string &str))
{
    miqp_solver.set_output_stream_fcn(funcptr);
}

// configure solver
void MPC_MIQP::configure_solver()
{
    // assemble qpipmpc_data
    qpipmpc_data.P_i_nom = P_i_nom;
    qpipmpc_data.P_i_vec = P_i_vec;
    qpipmpc_data.q_i_nom = q_i_nom;
    qpipmpc_data.q_i_vec = q_i_vec;
    qpipmpc_data.b = b;
    qpipmpc_data.C_i_nom = C_i_nom;
    qpipmpc_data.C_i_vec = C_i_vec;
    qpipmpc_data.D_i_nom = D_i_nom;
    qpipmpc_data.D_i_vec = D_i_vec;
    qpipmpc_data.crhs_i_nom = crhs_i_nom;
    qpipmpc_data.crhs_i_vec = crhs_i_vec;
    qpipmpc_data.G_i_nom = G_i_nom;
    qpipmpc_data.G_i_vec = G_i_vec;
    qpipmpc_data.w_i_nom = w_i_nom;
    qpipmpc_data.w_i_vec = w_i_vec;
    qpipmpc_data.idx_i_nom = idx_i_nom;
    qpipmpc_data.idx_i_vec = idx_i_vec;
    qpipmpc_data.idx_state = idx_state;
    qpipmpc_data.idx_input = idx_input;
    qpipmpc_data.idx_binvar = idx_binvar;
    qpipmpc_data.idx_y = idx_y;
    qpipmpc_data.idx_ineq = idx_ineq;
    qpipmpc_data.idx_eq = idx_eq;
    qpipmpc_data.n_horizon = n_horizon;
    qpipmpc_data.Pc = Pc;

    // setup MIQP solver
    if (Hrep_defined)
    {
        if (terminal_hybzono_constraint)
        {
            miqp_solver.setup(qpipmpc_data, miqp_settings, qp_settings,
                R_region2region, R_x02region, region_0, A_Hrep, b_Hrep, 
                C_term, D_term, G_term, P_term, crhs_term, w_term, q_term, idx_term, idx_binvar_term);
        }
        else
        {
            miqp_solver.setup(qpipmpc_data, miqp_settings, qp_settings,
                R_region2region, R_x02region, region_0, A_Hrep, b_Hrep);
        }
    }
    else
    {
        if (terminal_hybzono_constraint)
        {
            miqp_solver.setup(qpipmpc_data, miqp_settings, qp_settings,
                R_region2region, R_x02region, region_0,
                C_term, D_term, G_term, P_term, crhs_term, w_term, q_term, idx_term, idx_binvar_term);
        }
        else
        {
            miqp_solver.setup(qpipmpc_data, miqp_settings, qp_settings,
                R_region2region, R_x02region, region_0);
        }
    }
}

// solve optimization problem
bool MPC_MIQP::solve_optimization_problem()
{
    // solve
    miqp_results = miqp_solver.solve();

    solution = miqp_results.x;
    objective = miqp_results.upper_glob;

    // check status
    bool feasible = miqp_results.status != MI_QPIPMPC::MIQP_Status::MIQP_INFEASIBLE
        && miqp_results.status != MI_QPIPMPC::MIQP_Status::MIQP_NO_SOL
        && miqp_results.status != MI_QPIPMPC::MIQP_Status::MIQP_MAX_ITER_INFEASIBLE;

    // handle infeasibility
    get_state_input_trajectories();

    // get terminal region selection
    if (terminal_hybzono_constraint && feasible)
    {
        Eigen::VectorXd term_region_binvars_relaxed = solution(idx_binvar_term);
    
        // find index of coefficient with maximum value
        int max_idx = 0;
        double max_val = term_region_binvars_relaxed.maxCoeff();
        for (int i=0; i<term_region_binvars_relaxed.size(); ++i)
        {
            if (term_region_binvars_relaxed(i) == max_val)
            {
                max_idx = i; 
                break;
            }
        }
        terminal_region_selection = max_idx+1; // 1-indexed
    }

    // warm start
    if (warm_start)
    {
        if (feasible)
        {   
            generate_warm_start_nodes();
            term_region_warm_start = terminal_region_selection;
        }
        else
        {
            region_vec_warm_start.clear();
            term_region_warm_start = 0;
        }
    }

    // verbosity logging
    if (miqp_settings.verbose)
        verbose_output = miqp_solver.get_log();

    return feasible;
}

// get problem dimensions
void MPC_MIQP::get_problem_dimensions()
{
    // reset
    n_y_vec.clear();
    n_ineq_vec.clear();
    n_eq_vec.clear();
    idx_state.clear();
    idx_input.clear();
    idx_y.clear();
    idx_ineq.clear();
    idx_eq.clear();
    idx_binvar.clear();
    idx_binvar_term.clear();

    // full problem dimensions
    n_y_tot = 0; // init
    n_ineq_tot = 0; // init
    for (int k=0; k<=n_horizon; ++k)
    {
        if (k == n_horizon && terminal_hybzono_constraint)
            n_y_vec[k] = q_term.size();
        else if (q_i_vec.count(k))
            n_y_vec[k] = q_i_vec[k].size();
        else
            n_y_vec[k] = q_i_nom.size();
        n_y_tot += n_y_vec[k];

        if (k == n_horizon && terminal_hybzono_constraint)
            n_ineq_vec[k] = w_term.size();
        else if (w_i_vec.count(k))
            n_ineq_vec[k] = w_i_vec[k].size();
        else
            n_ineq_vec[k] = w_i_nom.size();
        n_ineq_tot += n_ineq_vec[k];
    }

    n_eq_tot = 0; // init
    for (int k=0; k<n_horizon; ++k)
    {
        if (k == n_horizon-1 && terminal_hybzono_constraint)
            n_eq_vec[k] = crhs_term.size();
        else if (crhs_i_vec.count(k))
            n_eq_vec[k] = crhs_i_vec[k].size();
        else
            n_eq_vec[k] = crhs_i_nom.size();
        n_eq_tot += n_eq_vec[k];
    }

    // indexing
    int idx_offset = 0; // init
    int idx_eq_offset = 0; // init
    int idx_ineq_offset = 0; // init
    for (int k=0; k<=n_horizon; ++k)
    {
        idx_state[k] = get_indices(idx_offset, n_state);
        idx_input[k] = get_indices(idx_offset + n_state, n_input);
        idx_y[k] = get_indices(idx_offset, n_y_vec[k]);
        idx_ineq[k] = get_indices(idx_ineq_offset, n_ineq_vec[k]);
        idx_binvar[k] = get_indices(idx_offset + n_y_vec_orig[k] + X_hz.nGc, X_hz.nGb);

        if (k < n_horizon)
            idx_eq[k] = get_indices(idx_eq_offset, n_eq_vec[k]);

        if (k == n_horizon)
        {
            idx_binvar_term = idx_term; // init
            for (auto it = idx_binvar_term.begin(); it != idx_binvar_term.end(); ++it)
                *it += idx_offset;
        }

        idx_offset += n_y_vec[k];
        idx_ineq_offset += n_ineq_vec[k];
        idx_eq_offset += n_eq_vec[k];
    }
}

// call parent setup methods
void MPC_MIQP::call_parent_setup_methods()
{
    MPC_QP::make_inequality_constraints();
    MPC_QP::make_equality_constraints();
    MPC_QP::make_cost_function();

    // get original problem dimensions
    for (int k=0; k<=n_horizon; ++k)
    {
        if (q_i_vec.count(k))
            n_y_vec_orig[k] = q_i_vec[k].size();
        else
            n_y_vec_orig[k] = q_i_nom.size();
    }

    // get original problem matrices
    P_i_nom_orig = P_i_nom;
    P_i_vec_orig = P_i_vec;
    q_i_nom_orig = q_i_nom;
    q_i_vec_orig = q_i_vec;
    C_i_nom_orig = C_i_nom;
    C_i_vec_orig = C_i_vec;
    D_i_nom_orig = D_i_nom;
    D_i_vec_orig = D_i_vec;
    crhs_i_nom_orig = crhs_i_nom;
    crhs_i_vec_orig = crhs_i_vec;
    G_i_nom_orig = G_i_nom;
    G_i_vec_orig = G_i_vec;
    w_i_nom_orig = w_i_nom;
    w_i_vec_orig = w_i_vec;
}

// make inequality constraints
void MPC_MIQP::make_inequality_constraints()
{
    // reset inequality constraints to original problem
    G_i_nom = G_i_nom_orig;
    G_i_vec = G_i_vec_orig;
    w_i_nom = w_i_nom_orig;
    w_i_vec = w_i_vec_orig;

    // add hybzono box constraints to nominal constraints
    int m_G_nom = G_i_nom.rows();
    int n_G_nom = G_i_nom.cols();
    int m_sigma_xi = Q_hz.rows();

    std::vector<Eigen::Triplet<double>> G_i_nom_triplets;
    G_i_nom_triplets.reserve(G_i_nom.nonZeros() + 2*X_hz.nG + 2*m_sigma_xi);

    if (softened_hybzono_constraints)
    {
        // state vector: y = [y0, xi_c, xi_b, sigma_xi]
        get_triplets_offset(G_i_nom, G_i_nom_triplets, 0, 0);
        get_triplets_offset_identity(X_hz.nG, G_i_nom_triplets, m_G_nom, n_G_nom);
        get_triplets_offset_neg_identity(X_hz.nG, G_i_nom_triplets, m_G_nom + X_hz.nG, n_G_nom);
        get_triplets_offset_identity(m_sigma_xi, G_i_nom_triplets, m_G_nom + 2*X_hz.nG, n_G_nom + X_hz.nG);
        get_triplets_offset_neg_identity(m_sigma_xi, G_i_nom_triplets, m_G_nom + 2*X_hz.nG + m_sigma_xi, n_G_nom + X_hz.nG);

        G_i_nom.resize(m_G_nom + 2*X_hz.nG + 2*m_sigma_xi, n_G_nom + X_hz.nG + m_sigma_xi);
        G_i_nom.setFromTriplets(G_i_nom_triplets.begin(), G_i_nom_triplets.end());

        w_i_nom.conservativeResize(m_G_nom + 2*X_hz.nG + 2*m_sigma_xi);
        w_i_nom.segment(m_G_nom, X_hz.nG) = Eigen::VectorXd::Ones(X_hz.nG);
        w_i_nom.segment(m_G_nom + X_hz.nG, X_hz.nG) = Eigen::VectorXd::Zero(X_hz.nG);
        w_i_nom.segment(m_G_nom + 2*X_hz.nG, m_sigma_xi) = sigma_max_hz*Eigen::VectorXd::Ones(m_sigma_xi);
        w_i_nom.segment(m_G_nom + 2*X_hz.nG + m_sigma_xi, m_sigma_xi) = -sigma_min_hz*Eigen::VectorXd::Ones(m_sigma_xi);
    }
    else
    {
        // state vector: y = [y0, xi_c, xi_b]
        get_triplets_offset(G_i_nom, G_i_nom_triplets, 0, 0);
        get_triplets_offset_identity(X_hz.nG, G_i_nom_triplets, m_G_nom, n_G_nom);
        get_triplets_offset_neg_identity(X_hz.nG, G_i_nom_triplets, m_G_nom + X_hz.nG, n_G_nom);

        G_i_nom.resize(m_G_nom + 2*X_hz.nG, n_G_nom + X_hz.nG);
        G_i_nom.setFromTriplets(G_i_nom_triplets.begin(), G_i_nom_triplets.end());

        w_i_nom.conservativeResize(m_G_nom + 2*X_hz.nG);
        w_i_nom.segment(m_G_nom, X_hz.nG) = Eigen::VectorXd::Ones(X_hz.nG);
        w_i_nom.segment(m_G_nom + X_hz.nG, X_hz.nG) = Eigen::VectorXd::Zero(X_hz.nG);
    }

    // add hybzono box constraints to off-nominal constraints
    std::vector<Eigen::Triplet<double>> G_i_offnom_triplets;
    G_i_offnom_triplets.reserve(2*G_i_nom.nonZeros());
    int m_G_offnom, n_G_offnom;

    for (int k=0; k<=n_horizon; ++k)
    {
        if (!G_i_vec.count(k))
            continue;

        m_G_offnom = G_i_vec[k].rows();
        n_G_offnom = G_i_vec[k].cols();
        G_i_offnom_triplets.clear();

        if (softened_hybzono_constraints)
        {
            // state vector: y = [y0, xi_c, xi_b, sigma_xi]
            get_triplets_offset(G_i_vec[k], G_i_offnom_triplets, 0, 0);
            get_triplets_offset_identity(X_hz.nG, G_i_offnom_triplets, m_G_offnom, n_G_offnom);
            get_triplets_offset_neg_identity(X_hz.nG, G_i_offnom_triplets, m_G_offnom + X_hz.nG, n_G_offnom);
            get_triplets_offset_identity(m_sigma_xi, G_i_offnom_triplets, m_G_offnom + 2*X_hz.nG, n_G_offnom + X_hz.nG);
            get_triplets_offset_neg_identity(m_sigma_xi, G_i_offnom_triplets, m_G_offnom + 2*X_hz.nG + m_sigma_xi, n_G_offnom + X_hz.nG);

            G_i_vec[k].resize(m_G_offnom + 2*X_hz.nG + 2*m_sigma_xi, n_G_offnom + X_hz.nG + m_sigma_xi);
            G_i_vec[k].setFromTriplets(G_i_offnom_triplets.begin(), G_i_offnom_triplets.end());
        }
        else
        {
            // state vector: y = [y0, xi_c, xi_b]
            get_triplets_offset(G_i_vec[k], G_i_offnom_triplets, 0, 0);
            get_triplets_offset_identity(X_hz.nG, G_i_offnom_triplets, m_G_offnom, n_G_offnom);
            get_triplets_offset_neg_identity(X_hz.nG, G_i_offnom_triplets, m_G_offnom + X_hz.nG, n_G_offnom);

            G_i_vec[k].resize(m_G_offnom + 2*X_hz.nG, n_G_offnom + X_hz.nG);
            G_i_vec[k].setFromTriplets(G_i_offnom_triplets.begin(), G_i_offnom_triplets.end());
        }
    }

    int w_size_offnom;
    for (int k=0; k<=n_horizon; ++k)
    {
        if (!w_i_vec.count(k))
            continue;

        w_size_offnom = w_i_vec[k].size();

        if (softened_hybzono_constraints)
        {
            w_i_vec[k].conservativeResize(w_i_vec[k].size() + 2*X_hz.nG + 2*m_sigma_xi);
            w_i_vec[k].segment(w_size_offnom, X_hz.nG) = Eigen::VectorXd::Ones(X_hz.nG);
            w_i_vec[k].segment(w_size_offnom + X_hz.nG, X_hz.nG) = Eigen::VectorXd::Zero(X_hz.nG);
            w_i_vec[k].segment(w_size_offnom + 2*X_hz.nG, m_sigma_xi) = sigma_max_hz*Eigen::VectorXd::Ones(m_sigma_xi);
            w_i_vec[k].segment(w_size_offnom + 2*X_hz.nG + m_sigma_xi, m_sigma_xi) = -sigma_min_hz*Eigen::VectorXd::Ones(m_sigma_xi);
        }
        else
        {
            w_i_vec[k].conservativeResize(w_i_vec[k].size() + 2*X_hz.nG);
            w_i_vec[k].segment(w_size_offnom, X_hz.nG) = Eigen::VectorXd::Ones(X_hz.nG);
            w_i_vec[k].segment(w_size_offnom + X_hz.nG, X_hz.nG) = Eigen::VectorXd::Zero(X_hz.nG);
        }
    }
}

// make equality constraints
void MPC_MIQP::make_equality_constraints()
{
    // reset equality constraints to original problem
    C_i_nom = C_i_nom_orig;
    C_i_vec = C_i_vec_orig;
    D_i_nom = D_i_nom_orig;
    D_i_vec = D_i_vec_orig;
    crhs_i_nom = crhs_i_nom_orig;
    crhs_i_vec = crhs_i_vec_orig;
    idx_i_nom.clear();
    idx_i_vec.clear();

    // LTV update if applicable
    if (LTV)
        update_equality_constraints_LTV();

    // add hybzono constraints to nominal constraints
    int m_C_nom = C_i_nom.rows();
    int n_C_nom = C_i_nom.cols();
    int m_D_nom = D_i_nom.rows();
    int m_sigma_xi = Q_hz.rows();

    // indices for binary vars
    idx_i_nom = get_indices(n_C_nom + X_hz.nGc, X_hz.nGb);

    // declare
    std::vector<Eigen::Triplet<double>> C_i_nom_triplets;
    C_i_nom_triplets.reserve(C_i_nom.nonZeros() + Pc.nonZeros() + X_hz.G.nonZeros() + m_sigma_xi + X_hz.A.nonZeros());

    // common additions to C_i_nom
    get_triplets_offset(C_i_nom, C_i_nom_triplets, 0, 0);
    get_triplets_offset(Pc, C_i_nom_triplets, m_C_nom, 0);
    get_triplets_offset(-1*X_hz.G, C_i_nom_triplets, m_C_nom, n_C_nom);
    get_triplets_offset(X_hz.A, C_i_nom_triplets, m_C_nom + X_hz.n, n_C_nom);

    // additions that depend on whether disturbances are included
    int i_r, j_r; // declare
    if (region_dependent_disturbances)
    {
        for (auto it = d_region.begin(); it != d_region.end(); ++it)
        {
            // TO DO: get indexing figured out
            i_r = idx_dyn_nom[0]; // just need first element
            j_r = idx_i_nom[it->first - 1]; // 0-indexed
            get_triplets_offset(it->second, C_i_nom_triplets, i_r, j_r);
        }
    }

    // additions that depend on softening
    if (softened_hybzono_constraints)
    {
        get_triplets_offset_identity(m_sigma_xi, C_i_nom_triplets, m_C_nom, n_C_nom + X_hz.nG);

        C_i_nom.resize(m_C_nom + X_hz.n + X_hz.nC, n_C_nom + X_hz.nG + X_hz.n);
        D_i_nom.conservativeResize(m_D_nom + X_hz.n + X_hz.nC, n_C_nom + X_hz.nG + X_hz.n);
    }
    else
    {
        C_i_nom.resize(m_C_nom + X_hz.n + X_hz.nC, n_C_nom + X_hz.nG);
        D_i_nom.conservativeResize(m_D_nom + X_hz.n + X_hz.nC, n_C_nom + X_hz.nG);
    }
    C_i_nom.setFromTriplets(C_i_nom_triplets.begin(), C_i_nom_triplets.end());

    crhs_i_nom.conservativeResize(m_C_nom + X_hz.n + X_hz.nC);
    crhs_i_nom.segment(m_C_nom, X_hz.n) = -X_hz.c;
    crhs_i_nom.segment(m_C_nom + X_hz.n, X_hz.nC) = -X_hz.b;

    // add hybzono inequality constraints to off-nominal constraints
    std::vector<Eigen::Triplet<double>> C_i_offnom_triplets;
    C_i_offnom_triplets.reserve(2*C_i_nom.nonZeros());

    for (int k=0; k<n_horizon; ++k)
    {
        if (!C_i_vec.count(k))
            continue;

        int m_C_offnom = C_i_vec[k].rows();
        int n_C_offnom = C_i_vec[k].cols();
        C_i_offnom_triplets.clear();

        // indexing
        idx_i_vec[k] = get_indices(n_C_offnom + X_hz.nGc, X_hz.nGb);

        // common additions
        get_triplets_offset(C_i_vec[k], C_i_offnom_triplets, 0, 0);
        get_triplets_offset(Pc, C_i_offnom_triplets, m_C_offnom, 0);
        get_triplets_offset(-1*X_hz.G, C_i_offnom_triplets, m_C_offnom, n_C_offnom);
        get_triplets_offset(X_hz.A, C_i_offnom_triplets, m_C_offnom + X_hz.n, n_C_offnom);

        // additions that depend on whether disturbances are included
        int i_r_k, j_r_k; // declare
        if (region_dependent_disturbances)
        {
            for (auto it = d_region.begin(); it != d_region.end(); ++it)
            {
                if (idx_dyn_vec.count(k))
                    i_r_k = idx_dyn_vec[k][0]; // just need first element
                else
                    i_r_k = idx_dyn_nom[0]; // just need first element
                if (idx_i_vec.count(k))
                    j_r_k = idx_i_vec[k][it->first - 1]; // 0-indexed
                else
                    j_r_k = idx_i_nom[it->first - 1]; // 0-indexed
                get_triplets_offset(it->second, C_i_offnom_triplets, i_r, j_r);
            }
        }

        // additions that depend on softening
        if (softened_hybzono_constraints)
        {
            get_triplets_offset_identity(m_sigma_xi, C_i_offnom_triplets, m_C_offnom, n_C_offnom + X_hz.nG);
            C_i_vec[k].resize(m_C_offnom + X_hz.n + X_hz.nC, n_C_offnom + X_hz.nG + X_hz.n);
        }
        else
        {
            C_i_vec[k].resize(m_C_offnom + X_hz.n + X_hz.nC, n_C_offnom + X_hz.nG);
        }
        C_i_vec[k].setFromTriplets(C_i_offnom_triplets.begin(), C_i_offnom_triplets.end());
    }
    for (int k=1; k<=n_horizon; ++k)
    {
        if (!D_i_vec.count(k))
            continue;

        int m_D_offnom = D_i_vec[k].rows();
        int n_D_offnom = D_i_vec[k].cols();

        if (softened_hybzono_constraints)
            D_i_vec[k].conservativeResize(m_D_offnom + X_hz.n + X_hz.nC, n_D_offnom + X_hz.nG + X_hz.n);
        else
            D_i_vec[k].conservativeResize(m_D_offnom + X_hz.n + X_hz.nC, n_D_offnom + X_hz.nG);

        idx_i_vec[k] = get_indices(n_D_offnom + X_hz.nGc, X_hz.nGb);
    }
    for (int k=0; k<n_horizon; ++k)
    {
        if (!crhs_i_vec.count(k))
            continue;

        int m_crhs_offnom = crhs_i_vec[k].size();
        crhs_i_vec[k].conservativeResize(m_crhs_offnom + X_hz.n + X_hz.nC);
        crhs_i_vec[k].segment(m_crhs_offnom, X_hz.n) = -X_hz.c;
        crhs_i_vec[k].segment(m_crhs_offnom + X_hz.n, X_hz.nC) = -X_hz.b;
    }

    // make terminal hybzono constraints
    make_terminal_hybzono_constraints();
}

// make terminal hybzono constraints
void MPC_MIQP::make_terminal_hybzono_constraints()
{
    // get D matrix at n_horizon
    if (!D_i_vec.count(n_horizon))
        D_i_vec[n_horizon] = D_i_nom;

    // matrix dimensions
    int m_D_N = D_i_vec[n_horizon].rows();
    int n_D_N = D_i_vec[n_horizon].cols();

    std::vector<Eigen::Triplet<double>> D_N_triplets;
    D_N_triplets.reserve(D_i_vec[n_horizon].nonZeros() + Pc.nonZeros() + X_hz.G.nonZeros() + X_hz.n + X_hz.A.nonZeros());

    if (softened_hybzono_constraints)
    {
        get_triplets_offset(D_i_vec[n_horizon], D_N_triplets, 0, 0);
        get_triplets_offset(Pc, D_N_triplets, m_D_N, 0);
        get_triplets_offset(-1*X_hz.G, D_N_triplets, m_D_N, n_y_vec_orig[n_horizon]);
        get_triplets_offset_identity(X_hz.n, D_N_triplets, m_D_N, n_y_vec_orig[n_horizon] + X_hz.nG);
        get_triplets_offset(X_hz.A, D_N_triplets, m_D_N + X_hz.n, n_y_vec_orig[n_horizon]);

        D_i_vec[n_horizon].resize(m_D_N + X_hz.n + X_hz.nC, n_D_N);
        D_i_vec[n_horizon].setFromTriplets(D_N_triplets.begin(), D_N_triplets.end());
    }
    else
    {
        get_triplets_offset(D_i_vec[n_horizon], D_N_triplets, 0, 0);
        get_triplets_offset(Pc, D_N_triplets, m_D_N, 0);
        get_triplets_offset(-1*X_hz.G, D_N_triplets, m_D_N, n_y_vec_orig[n_horizon]);
        get_triplets_offset(X_hz.A, D_N_triplets, m_D_N + X_hz.n, n_y_vec_orig[n_horizon]);

        D_i_vec[n_horizon].resize(m_D_N + X_hz.n + X_hz.nC, n_D_N);
        D_i_vec[n_horizon].setFromTriplets(D_N_triplets.begin(), D_N_triplets.end());
    }

    // get crhs vector at n_horizon-1
    if (!crhs_i_vec.count(n_horizon-1))
        crhs_i_vec[n_horizon-1] = crhs_i_nom;
    
    int m_crhs_Nm1 = crhs_i_vec[n_horizon-1].size();
    crhs_i_vec[n_horizon-1].conservativeResize(m_crhs_Nm1 + X_hz.n + X_hz.nC);
    crhs_i_vec[n_horizon-1].segment(m_crhs_Nm1, X_hz.n) = -X_hz.c;
    crhs_i_vec[n_horizon-1].segment(m_crhs_Nm1 + X_hz.n, X_hz.nC) = -X_hz.b;

    // get C at n_horizon-1
    if (!C_i_vec.count(n_horizon-1))
        C_i_vec[n_horizon-1] = C_i_nom;

    int m_C_Nm1 = C_i_vec[n_horizon-1].rows();
    int n_C_Nm1 = C_i_vec[n_horizon-1].cols();
    C_i_vec[n_horizon-1].conservativeResize(m_C_Nm1 + X_hz.n + X_hz.nC, n_C_Nm1);
}

// make cost function matrices
void MPC_MIQP::make_cost_function()
{
    // reset cost function matrices to original problem
    P_i_nom = P_i_nom_orig;
    P_i_vec = P_i_vec_orig;
    q_i_nom = q_i_nom_orig;
    q_i_vec = q_i_vec_orig;

    // add hybzono cost to nominal cost
    int n_P_nom = P_i_nom.cols();
    int n_q_nom = q_i_nom.size();
    int n_sigma_xi = Q_hz.rows();

    std::vector<Eigen::Triplet<double>> P_i_nom_triplets;
    P_i_nom_triplets.reserve(P_i_nom.nonZeros() + Q_hz.nonZeros());

    if (softened_hybzono_constraints)
    {
        get_triplets_offset(P_i_nom, P_i_nom_triplets, 0, 0);
        get_triplets_offset(Q_hz, P_i_nom_triplets, n_P_nom + X_hz.nG, n_P_nom + X_hz.nG);
        P_i_nom.resize(n_P_nom + X_hz.nG + n_sigma_xi, n_P_nom + X_hz.nG + n_sigma_xi);
        P_i_nom.setFromTriplets(P_i_nom_triplets.begin(), P_i_nom_triplets.end());

        q_i_nom.conservativeResize(n_q_nom + X_hz.nG + n_sigma_xi);
        q_i_nom.segment(n_q_nom, X_hz.nG + n_sigma_xi) = Eigen::VectorXd::Zero(X_hz.nG + n_sigma_xi);
    }
    else
    {
        get_triplets_offset(P_i_nom, P_i_nom_triplets, 0, 0);
        P_i_nom.resize(n_P_nom + X_hz.nG, n_P_nom + X_hz.nG);
        P_i_nom.setFromTriplets(P_i_nom_triplets.begin(), P_i_nom_triplets.end());

        q_i_nom.conservativeResize(n_q_nom + X_hz.nG);
        q_i_nom.segment(n_q_nom, X_hz.nG) = Eigen::VectorXd::Zero(X_hz.nG);
    }

    // region costs
    if (region_costs)
    {
        q_i_nom.segment(n_q_nom + X_hz.nGc, X_hz.nGb) = region_cost_vec;
    }

    // add hybzono cost to off-nominal cost
    std::vector<Eigen::Triplet<double>> P_i_offnom_triplets;
    P_i_offnom_triplets.reserve(2*P_i_nom.nonZeros());
    for (int k=0; k<=n_horizon; ++k)
    {
        if (!P_i_vec.count(k))
            continue;

        int n_P_offnom = P_i_vec[k].cols();

        P_i_offnom_triplets.clear();

        if (softened_hybzono_constraints)
        {
            get_triplets_offset(P_i_vec[k], P_i_offnom_triplets, 0, 0);
            get_triplets_offset(Q_hz, P_i_offnom_triplets, n_P_offnom + X_hz.nG, n_P_offnom + X_hz.nG);
            P_i_vec[k].resize(n_P_offnom + X_hz.nG + n_sigma_xi, n_P_offnom + X_hz.nG + n_sigma_xi);
            P_i_vec[k].setFromTriplets(P_i_offnom_triplets.begin(), P_i_offnom_triplets.end());
        }
        else
        {
            get_triplets_offset(P_i_vec[k], P_i_offnom_triplets, 0, 0);
            P_i_vec[k].resize(n_P_offnom + X_hz.nG, n_P_offnom + X_hz.nG);
            P_i_vec[k].setFromTriplets(P_i_offnom_triplets.begin(), P_i_offnom_triplets.end());
        }
    }

    for (int k=0; k<=n_horizon; ++k)
    {
        if (!q_i_vec.count(k))
            continue;

        int n_q_offnom = q_i_vec[k].size();

        if (softened_hybzono_constraints)
        {
            q_i_vec[k].conservativeResize(n_q_offnom + X_hz.nG + n_sigma_xi);
            q_i_vec[k].segment(n_q_offnom, X_hz.nG + n_sigma_xi) = Eigen::VectorXd::Zero(X_hz.nG + n_sigma_xi);
        }
        else
        {
            q_i_vec[k].conservativeResize(n_q_offnom + X_hz.nG);
            q_i_vec[k].segment(n_q_offnom, X_hz.nG) = Eigen::VectorXd::Zero(X_hz.nG);
        }

        if (region_costs)
        {
            q_i_vec[k].segment(n_q_offnom + X_hz.nGc, X_hz.nGb) = region_cost_vec;
        }
    }
}

void MPC_MIQP::make_terminal_region_equality_constraints()
{
    // check if terminal region is set
    if (!terminal_hybzono_constraint) return;

    // get original equality constraints at last time step
    if (C_i_vec.count(n_horizon-1))
        C_term = C_i_vec[n_horizon-1];
    else
        C_term = C_i_nom;
    if (D_i_vec.count(n_horizon))
        D_term = D_i_vec[n_horizon];
    else
        D_term = D_i_nom;
    if (crhs_i_vec.count(n_horizon-1))
        crhs_term = crhs_i_vec[n_horizon-1];
    else
        crhs_term = crhs_i_nom;

    // get dimensions
    int m_C = C_term.rows();
    int n_C = C_term.cols();
    int m_D = D_term.rows();
    int n_D = D_term.cols();

    // add terminal region constraints

    // state vector: y = [y0, xi_c, xi_b, sigma_xi]

    // equality constraints
    std::vector<Eigen::Triplet<double>> D_term_triplets;
    D_term_triplets.reserve(D_term.nonZeros() + Pp_term.nonZeros() + X_hz_term.G.nonZeros() + X_hz_term.n + X_hz_term.A.nonZeros());
    get_triplets_offset(D_term, D_term_triplets, 0, 0);
    get_triplets_offset(Pp_term, D_term_triplets, m_D, 0);
    get_triplets_offset(-1*X_hz_term.G, D_term_triplets, m_D, n_D);
    get_triplets_offset_identity(X_hz_term.n, D_term_triplets, m_D, n_D + X_hz_term.nG);
    get_triplets_offset(X_hz_term.A, D_term_triplets, m_D + X_hz_term.n, n_D);
    D_term.resize(m_D + X_hz_term.n + X_hz_term.nC, n_D + X_hz_term.nG + X_hz_term.n);
    D_term.setFromTriplets(D_term_triplets.begin(), D_term_triplets.end());

    C_term.conservativeResize(m_C + X_hz_term.n + X_hz_term.nC, n_C);

    crhs_term.conservativeResize(m_C + X_hz_term.n + X_hz_term.nC);
    crhs_term.segment(m_C, X_hz_term.n) = -X_hz_term.c;
    crhs_term.segment(m_C + X_hz_term.n, X_hz_term.nC) = -X_hz_term.b;

    // indexing
    idx_term = get_indices(n_D + X_hz_term.nGc, X_hz_term.nGb);
}

// terminal region cost
void MPC_MIQP::make_terminal_region_cost()
{
    // check if terminal region is set
    if (!terminal_hybzono_constraint) return;

    // get original cost matrices at last time step
    if (P_i_vec.count(n_horizon))
        P_term = P_i_vec[n_horizon];
    else    
        P_term = P_i_nom;
    if (q_i_vec.count(n_horizon))
        q_term = q_i_vec[n_horizon];
    else
        q_term = q_i_nom;

    // get dimensions
    int n_P = P_term.rows();
    int n_q = q_term.size();

    // add terminal region cost

    // state vector: y = [y0, xi_c, xi_b, sigma_xi]

    // cost function
    std::vector<Eigen::Triplet<double>> P_term_triplets;
    P_term_triplets.reserve(P_term.nonZeros() + Q_term_slack_cost.nonZeros());
    get_triplets_offset(P_term, P_term_triplets, 0, 0);
    get_triplets_offset(Q_term_slack_cost, P_term_triplets, n_P + X_hz_term.nG, n_P + X_hz_term.nG);
    P_term.resize(n_P + X_hz_term.nG + X_hz_term.n, n_P + X_hz_term.nG + X_hz_term.n);
    P_term.setFromTriplets(P_term_triplets.begin(), P_term_triplets.end());

    q_term.conservativeResize(n_q + X_hz_term.nG + X_hz_term.n);
    q_term.segment(n_q, X_hz_term.nGc) = Eigen::VectorXd::Zero(X_hz_term.nGc);
    q_term.segment(n_q + X_hz_term.nGc, X_hz_term.nGb) = term_region_cost_vec;
    q_term.segment(n_q + X_hz_term.nG, X_hz_term.n) = Eigen::VectorXd::Zero(X_hz_term.n);
}

// make terminal inequality constraints
void MPC_MIQP::make_terminal_region_inequality_constraints()
{
    // check if terminal region is set
    if (!terminal_hybzono_constraint) return;

    // get original inequality constraints at last time step
    if (G_i_vec.count(n_horizon))
        G_term = G_i_vec[n_horizon];
    else
        G_term = G_i_nom;
    if (w_i_vec.count(n_horizon))
        w_term = w_i_vec[n_horizon];
    else
        w_term = w_i_nom;

    // get dimensions
    int m_G = G_term.rows();
    int n_G = G_term.cols();
    int m_w = w_term.size();

    // hybzono inequality constraints
    std::vector<Eigen::Triplet<double>> G_term_triplets;
    G_term_triplets.reserve(G_term.nonZeros() + 2*X_hz_term.nG + 2*X_hz_term.n);
    get_triplets_offset(G_term, G_term_triplets, 0, 0);
    get_triplets_offset_identity(X_hz_term.nG, G_term_triplets, m_G, n_G);
    get_triplets_offset_neg_identity(X_hz_term.nG, G_term_triplets, m_G + X_hz_term.nG, n_G);
    get_triplets_offset_identity(X_hz_term.n, G_term_triplets, m_G + 2*X_hz_term.nG, n_G + X_hz_term.nG);
    get_triplets_offset_neg_identity(X_hz_term.n, G_term_triplets, m_G + 2*X_hz_term.nG + X_hz_term.n, n_G + X_hz_term.nG);
    G_term.resize(m_G + 2*X_hz_term.nG + 2*X_hz_term.n, n_G + X_hz_term.nG + X_hz_term.n);
    G_term.setFromTriplets(G_term_triplets.begin(), G_term_triplets.end());

    // right hand side
    w_term.conservativeResize(m_w + 2*X_hz_term.nG + 2*X_hz_term.n);
    w_term.segment(m_w, X_hz_term.nGc) = Eigen::VectorXd::Ones(X_hz_term.nGc);
    w_term.segment(m_w + X_hz_term.nGc, X_hz_term.nGb) = Eigen::VectorXd::Ones(X_hz_term.nGb);
    w_term.segment(m_w + X_hz_term.nG, X_hz_term.nGc) = Eigen::VectorXd::Ones(X_hz_term.nGc);
    w_term.segment(m_w + X_hz_term.nG + X_hz_term.nGc, X_hz_term.nGb) = Eigen::VectorXd::Zero(X_hz_term.nGb);
    w_term.segment(m_w + 2*X_hz_term.nG, X_hz_term.n) = sigma_max_hz*Eigen::VectorXd::Ones(X_hz_term.n);
    w_term.segment(m_w + 2*X_hz_term.nG + X_hz_term.n, X_hz_term.n) = -sigma_min_hz*Eigen::VectorXd::Ones(X_hz_term.n);
}

// reachability calculations

// comfigure point 2 region
void MPC_MIQP::configure_point2region()
{
    // region
    int region = 1; // init

    // declare
    Eigen::SparseMatrix<double> P, A, G;
    Eigen::VectorXd x, q, b, w, c_i, xi_b_i;
    std::vector<Eigen::Triplet<double>> triplets;

    // problem definition
    if (Hrep_defined)
    {
        P.resize(n_pos_dim, n_pos_dim);
        P.setIdentity();
        x.resize(n_state);
        x.setZero();
        q = -1*Pp*x;

        A.resize(0, 0);
        b.resize(0);

        G = A_Hrep_pos[region].sparseView();
        w = b_Hrep_pos[region];
    }
    else
    {
        P.resize(X_hz.nGc, X_hz.nGc);
        P = X_hz.Gc.transpose()*X_hz.Gc;
        x.resize(n_state);
        x.setZero();
        
        xi_b_i.resize(X_hz.nGb);
        xi_b_i.setZero();
        xi_b_i(0) = 1;
        c_i = X_hz.Gb*xi_b_i + X_hz.c;

        q = ((c_i - Pp*x).transpose())*X_hz.Gc;

        A = X_hz.Ac;
        b = X_hz.b - X_hz.Ab*xi_b_i;

        G.resize(2*X_hz.nGc, X_hz.nGc);
        triplets.reserve(2*X_hz.nGc);
        for (int i=0; i<X_hz.nGc; ++i)
        {
            triplets.push_back(Eigen::Triplet<double>(i, i, 1));
            triplets.push_back(Eigen::Triplet<double>(i + X_hz.nGc, i, -1));
        }
        G.setFromTriplets(triplets.begin(), triplets.end());
        w = Eigen::VectorXd::Ones(2*X_hz.nGc);
    }

    // setup solver
    qp_solver_pt2region.setup(P, q, A, b, G, w);

    // set flag
    point2region_configured = true;
}

// get distance of point x to specified region
double MPC_MIQP::dist_point2region(const Eigen::Ref<const Eigen::VectorXd> x, int region)
{
    // configure solver if necessary
    if (!point2region_configured)
        configure_point2region();

    // declare
    Eigen::VectorXd q, w, x_region, xi_region, xi_b_i, c_i, b;
    Eigen::SparseMatrix<double> G;
    QP::QP_results results;

    if (Hrep_defined)
    {
        // update gradient
        q = -1*Pp*x;
        qp_solver_pt2region.set_q(q);

        // update constraints
        G = A_Hrep_pos[region].sparseView();
        w = b_Hrep_pos[region];
        qp_solver_pt2region.set_G(G);
        qp_solver_pt2region.set_w(w);

        // solve
        results = qp_solver_pt2region.solve();
        x_region = results.x;
    }
    else
    {
        // update gradient
        xi_b_i.resize(X_hz.nGb);
        xi_b_i.setZero();
        xi_b_i(region-1) = 1; // 1-based to 0-based indexing
        c_i = X_hz.Gb*xi_b_i + X_hz.c;

        q = ((c_i - Pp*x).transpose())*X_hz.Gc;
        qp_solver_pt2region.set_q(q);

        // update constraints
        b = X_hz.b - X_hz.Ab*xi_b_i;
        qp_solver_pt2region.set_b(b);

        // solve
        results = qp_solver_pt2region.solve();
        xi_region = results.x;

        // get x_region
        x_region = X_hz.Gc*xi_region + X_hz.Gb*xi_b_i + X_hz.c;   
    }

    // get distance from x to region
    return (Pp*x - x_region).norm();
}

// configure region to region calculations
void MPC_MIQP::configure_region2region()
{
    // region
    int region1 = 1; // init
    int region2 = 1; // init

    // declare
    Eigen::SparseMatrix<double> P, A, G, ImI;
    Eigen::VectorXd x, q, b, w, c_i, xi_b_i_1, xi_b_i_2;
    std::vector<Eigen::Triplet<double>> triplets;
    Eigen::MatrixXd P_full, A_full, Gc_full, GcTGc;

    if (Hrep_defined)
    {
        // problem definition
        P.resize(2*n_pos_dim, 2*n_pos_dim);
        triplets.clear();
        get_triplets_offset_identity(n_pos_dim, triplets, 0, 0);
        get_triplets_offset_neg_identity(n_pos_dim, triplets, 0, n_pos_dim);
        get_triplets_offset_neg_identity(n_pos_dim, triplets, n_pos_dim, 0);
        get_triplets_offset_identity(n_pos_dim, triplets, n_pos_dim, n_pos_dim);
        P.setFromTriplets(triplets.begin(), triplets.end());

        q = Eigen::VectorXd::Zero(2*n_pos_dim);

        A.resize(0, 0);
        b.resize(0);

        G.resize(A_Hrep_pos[region1].rows() + A_Hrep_pos[region2].rows(), 2*n_pos_dim);
        triplets.clear();
        get_triplets_offset(A_Hrep_pos[region1].sparseView(), triplets, 0, 0);
        get_triplets_offset(A_Hrep_pos[region2].sparseView(), triplets, A_Hrep_pos[region1].rows(), n_pos_dim);
        G.setFromTriplets(triplets.begin(), triplets.end());

        w.resize(b_Hrep_pos[region1].rows() + b_Hrep_pos[region2].rows());
        w.segment(0, b_Hrep_pos[region1].rows()) = b_Hrep_pos[region1];
        w.segment(b_Hrep_pos[region1].rows(), b_Hrep_pos[region2].rows()) = b_Hrep_pos[region2];
    }
    else
    {
        // cost
        P_full.resize(2*X_hz.nGc, 2*X_hz.nGc);
        Gc_full = Eigen::MatrixXd(X_hz.Gc);
        GcTGc = Gc_full.transpose()*Gc_full;
        P_full.block(0, 0, X_hz.nGc, X_hz.nGc) = GcTGc;
        P_full.block(0, X_hz.nGc, X_hz.nGc, X_hz.nGc) = -1*GcTGc;
        P_full.block(X_hz.nGc, 0, X_hz.nGc, X_hz.nGc) = -1*GcTGc;
        P_full.block(X_hz.nGc, X_hz.nGc, X_hz.nGc, X_hz.nGc) = GcTGc;
        P = P_full.sparseView();

        xi_b_i_1.resize(X_hz.nGb);
        xi_b_i_1.setZero();
        xi_b_i_1(0) = 1;
        xi_b_i_2 = xi_b_i_1;

        triplets.reserve(2*X_hz.nGc);
        triplets.clear();
        ImI.resize(X_hz.nGc, 2*X_hz.nGc);
        get_triplets_offset_identity(X_hz.nGc, triplets, 0, 0);
        get_triplets_offset_neg_identity(X_hz.nGc, triplets, 0, X_hz.nGc);
        ImI.setFromTriplets(triplets.begin(), triplets.end());

        q = (xi_b_i_1 - xi_b_i_2).transpose()*X_hz.Gb.transpose()*X_hz.Gc*ImI;

        // constraints
        A_full.resize(2*X_hz.nC, 2*X_hz.nGc);
        A_full.setZero();
        A_full.block(0, 0, X_hz.nC, X_hz.nGc) = Eigen::MatrixXd(X_hz.Ac);
        A_full.block(X_hz.nC, X_hz.nGc, X_hz.nC, X_hz.nGc) = Eigen::MatrixXd(X_hz.Ac);
        A = A_full.sparseView();

        b.resize(2*X_hz.nC);
        b.segment(0, X_hz.nC) = X_hz.b - X_hz.Ab*xi_b_i_1;
        b.segment(X_hz.nC, X_hz.nC) = X_hz.b - X_hz.Ab*xi_b_i_2;

        G.resize(4*X_hz.nGc, 2*X_hz.nGc);
        triplets.clear();
        for (int i=0; i<2*X_hz.nGc; ++i)
        {
            triplets.push_back(Eigen::Triplet<double>(i, i, 1));
            triplets.push_back(Eigen::Triplet<double>(i + 2*X_hz.nGc, i, -1));
        }
        G.setFromTriplets(triplets.begin(), triplets.end());
        w = Eigen::VectorXd::Ones(4*X_hz.nGc);
    }

    // setup solver
    qp_solver_region2region.setup(P, q, A, b, G, w);

    // set flag
    region2region_configured = true;
}

// get distance between regions
double MPC_MIQP::dist_region2region(int region1, int region2)
{
    // configure solver if necessary
    if (!region2region_configured)
        configure_region2region();

    // declare
    Eigen::SparseMatrix<double> G, ImI;
    Eigen::VectorXd q, w, b, x_region1, x_region2, xi_region1, xi_region2, xi_b_i_1, xi_b_i_2;
    std::vector<Eigen::Triplet<double>> triplets;
    QP::QP_results results;

    if (Hrep_defined)
    {
        // update constraints
        G.resize(A_Hrep_pos[region1].rows() + A_Hrep_pos[region2].rows(), 2*n_pos_dim);
        triplets.reserve(G.rows()*G.cols());
        triplets.clear();
        get_triplets_offset(A_Hrep_pos[region1].sparseView(), triplets, 0, 0);
        get_triplets_offset(A_Hrep_pos[region2].sparseView(), triplets, A_Hrep_pos[region1].rows(), n_pos_dim);
        G.setFromTriplets(triplets.begin(), triplets.end());

        w.resize(b_Hrep_pos[region1].rows() + b_Hrep_pos[region2].rows());
        w.segment(0, b_Hrep_pos[region1].rows()) = b_Hrep_pos[region1];
        w.segment(b_Hrep_pos[region1].rows(), b_Hrep_pos[region2].rows()) = b_Hrep_pos[region2];
        qp_solver_region2region.set_G(G);
        qp_solver_region2region.set_w(w);

        // solve
        results = qp_solver_region2region.solve();
        x_region1 = results.x.segment(0, n_pos_dim);
        x_region2 = results.x.segment(n_pos_dim, n_pos_dim);        
    }
    else
    {
        // update gradient
        xi_b_i_1.resize(X_hz.nGb);
        xi_b_i_1.setZero();
        xi_b_i_2 = xi_b_i_1;
        xi_b_i_1(region1-1) = 1; // one-based indexing to zero-based
        xi_b_i_2(region2-1) = 1; // one-based indexing to zero-based

        triplets.reserve(2*X_hz.nGc);
        triplets.clear();
        ImI.resize(X_hz.nGc, 2*X_hz.nGc);
        get_triplets_offset_identity(X_hz.nGc, triplets, 0, 0);
        get_triplets_offset_neg_identity(X_hz.nGc, triplets, 0, X_hz.nGc);
        ImI.setFromTriplets(triplets.begin(), triplets.end());

        q = (xi_b_i_1 - xi_b_i_2).transpose()*X_hz.Gb.transpose()*X_hz.Gc*ImI;
        qp_solver_region2region.set_q(q);

        // update constraints
        b.resize(2*X_hz.nC);
        b.segment(0, X_hz.nC) = X_hz.b - X_hz.Ab*xi_b_i_1;
        b.segment(X_hz.nC, X_hz.nC) = X_hz.b - X_hz.Ab*xi_b_i_2;

        qp_solver_region2region.set_b(b);

        // solve
        results = qp_solver_region2region.solve();
        xi_region1 = results.x.segment(0, X_hz.nGc);
        xi_region2 = results.x.segment(X_hz.nGc, X_hz.nGc);

        // get x_region
        x_region1 = X_hz.Gc*xi_region1 + X_hz.Gb*xi_b_i_1 + X_hz.c;   
        x_region2 = X_hz.Gc*xi_region2 + X_hz.Gb*xi_b_i_2 + X_hz.c;   
    }

    // get distance from x to region
    return (x_region1 - x_region2).norm();
}

// get x0 reachability
void MPC_MIQP::get_x0_reachability()
{
    // declare
    std::vector<int> possible_regions;

    // check if distances specified
    if (!point2region_dist_computed)
    {
        // init
        d_x02region.clear();
        for (int i=0; i<n_regions; i++)
            d_x02region.push_back(QP::inf);

        // if already initialized, only explore regions that are known to be reachable from possible current regions
        if (region_0_valid)
        {
            // possible regions at current time step
            std::vector<int> possible_start_regions = R_x02region[1];

            // get all potentially reachable regions given possible current regions
            for (int region : possible_start_regions)
            {
                for (auto it = R_region2region[n_horizon][region].begin(); it != R_region2region[n_horizon][region].end(); ++it)
                {
                    if (std::find(possible_regions.begin(), possible_regions.end(), *it) == possible_regions.end())
                        possible_regions.push_back(*it);
                }
            }

            // get distance to reachable regions only
            for (auto it = possible_regions.begin(); it != possible_regions.end(); ++it)
                d_x02region[*it-1] = dist_point2region(x0, *it); // indexing correction
        }
        else
        {
            // get distance to all regions
            for (int region=1; region<=n_regions; ++region)
                d_x02region[region-1] = dist_point2region(x0, region); // indexing correction
        }
    }

    // check for system currently violating constraints
    // find min of d_x02region
    double d_min = *std::min_element(d_x02region.begin(), d_x02region.end());
    double d_max;
    if (d_min > d_lim)
        d_max = 1.1*d_min;
    else
        d_max = d_lim;

    // build reachability matrix
    for (int k=1; k<=n_horizon; ++k)
    {
        // reset
        R_x02region[k].clear();

        // get elements of d_x02region that are less than d_max
        for (int region=1; region<=n_regions; ++region)
        {
            if (d_x02region[region-1] <= k*d_max) // indexing correction
                R_x02region[k].push_back(region);
        }
    }

    // get initial region (region with minimum distance)
    region_0 = std::find(d_x02region.begin(), d_x02region.end(), d_min) - d_x02region.begin() + 1; // one-based indexing of regions

    // reset flag
    point2region_dist_computed = false;
}

// get reachability between regions as a function of the number of time steps
void MPC_MIQP::get_region2region_reachability()
{
    if (!region2region_dist_computed)
    {
        // first get distance between each pair of regions
        d_mat = Eigen::MatrixXd::Zero(n_regions, n_regions); // init
        for (int region1=1; region1<=n_regions; ++region1)
        {
            for (int region2=region1+1; region2<=n_regions; ++region2)
            {
                double dist = dist_region2region(region1, region2);
                d_mat(region1-1, region2-1) = dist;
                d_mat(region2-1, region1-1) = dist;
            }
        }
        region2region_dist_computed = true;
    }

    // build reachability matrix
    R_region2region.clear(); // reset
    for (int k=1; k<=n_horizon; ++k)
    {
        // get elements of d_mat_global that are less than d_lim
        for (int region1=1; region1<=n_regions; ++region1)
        {
            for (int region2=1; region2<=n_regions; ++region2)
            {
                if (d_mat(region1-1, region2-1) <= k*d_lim)
                {
                    R_region2region[k][region1].push_back(region2);
                }
            }
        }
    }
}

// verbose output
std::string MPC_MIQP::get_verbose_output()
{
    return verbose_output;
}

// generate warm start nodes
void MPC_MIQP::generate_warm_start_nodes()
{
    // get all regions that are one-step reachable from terminal region in region_vec_global
    std::vector<int> reachable_regions;
    for (auto it = R_region2region[1][miqp_results.region_vec.back()].begin(); it != R_region2region[1][miqp_results.region_vec.back()].end(); ++it)
        reachable_regions.push_back(*it);

    // construct warm start nodes
    region_vec_warm_start.clear();
    region_vec_warm_start.reserve(reachable_regions.size());

    std::vector<int> common_regions;
    auto it = miqp_results.region_vec.begin();
    it++; // skip first element
    while (it != miqp_results.region_vec.end())
    {
        common_regions.push_back(*it);
        ++it;
    }

    for (auto it = reachable_regions.begin(); it != reachable_regions.end(); ++it)
    {
        region_vec_warm_start.push_back(common_regions);
        region_vec_warm_start.back().push_back(*it);
    }
}

void MPC_MIQP::apply_region0_to_warm_start_nodes()
{
    for (auto it = region_vec_warm_start.begin(); it != region_vec_warm_start.end(); ++it)
        *(it->begin()) = region_0;
}