#include "MPC_MIQP_Hrep.hpp"

using namespace HybZonoMPC;

// constructor
MPC_MIQP_Hrep::MPC_MIQP_Hrep()
{
    Hrep_defined = true;
}

// setup H-rep constraints
void MPC_MIQP_Hrep::setup_Hrep_constraints(const std::map<int, Eigen::MatrixXd> &A_Hrep, const std::map<int, Eigen::VectorXd> &b_Hrep,
    const Eigen::SparseMatrix<double> &Pc, double big_M, const Eigen::SparseMatrix<double> &Q_hr)
{
    // constraint softening
    this->Q_hr = Q_hr;
    this->softened_hrep_constraints = true;

    // call other setup method
    setup_Hrep_constraints(A_Hrep, b_Hrep, Pc, big_M);
}

void MPC_MIQP_Hrep::setup_Hrep_constraints(const std::map<int, Eigen::MatrixXd> &A_Hrep, const std::map<int, Eigen::VectorXd> &b_Hrep,
    const Eigen::SparseMatrix<double> &Pc, double big_M, const Eigen::SparseMatrix<double> &Q_hr, const Eigen::Ref<const Eigen::VectorXd> region_cost_vec)
{
    // constraint softening
    this->Q_hr = Q_hr;
    this->softened_hrep_constraints = true;

    // region costs
    this->region_cost_vec = region_cost_vec;
    this->region_costs = true;

    // call other setup method
    setup_Hrep_constraints(A_Hrep, b_Hrep, Pc, big_M);
}

void MPC_MIQP_Hrep::setup_Hrep_constraints(const std::map<int, Eigen::MatrixXd> &A_Hrep, const std::map<int, Eigen::VectorXd> &b_Hrep,
        const Eigen::SparseMatrix<double> &Pp, double big_M,
        const Eigen::Ref<const Eigen::VectorXd> region_cost_vec)
{
    // region costs
    this->region_cost_vec = region_cost_vec;
    this->region_costs = true;

    // call other setup method
    setup_Hrep_constraints(A_Hrep, b_Hrep, Pc, big_M);
}

void MPC_MIQP_Hrep::setup_Hrep_constraints(const std::map<int, Eigen::MatrixXd> &A_Hrep, const std::map<int, Eigen::VectorXd> &b_Hrep,
    const Eigen::SparseMatrix<double> &Pc, double big_M)
{
    // set up H-rep constraints
    this->A_Hrep = A_Hrep;
    this->b_Hrep = b_Hrep;

    // number of regions
    this->n_regions = A_Hrep.size();

    // number of H-rep constraints
    n_cons = 0;
    for (auto const& [key, val] : A_Hrep)
    {
        n_cons += val.rows();
    }

    // state map to constraint domain
    this->Pc = Pc;

    // big M
    this->M = big_M;

    // position info
    this->Pp = Pc;
    this->A_Hrep_pos = A_Hrep;
    this->b_Hrep_pos = b_Hrep;

    // dimension of constraint space
    this->n_pos_dim = this->Pp.rows();

    // set flag
    Hrep_setup = true;
}

// build controller
void MPC_MIQP_Hrep::build_controller()
{
    // check for invalid problem setup for H-rep case
    if (!Hrep_setup)
        throw std::invalid_argument("H-rep constraints not set up");

    if (terminal_hybzono_constraint)
        throw std::invalid_argument("terminal hybzono constraints not supported for H-rep");

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

    // initialize x0 and x_ref
    x0 = Eigen::VectorXd::Zero(n_state);
    x_ref = std::vector<Eigen::VectorXd>(n_horizon, Eigen::VectorXd::Zero(n_state));

    // check problem validity
    auto [valid, msg] = check_problem_validity();
    if (!valid)
    {
        throw std::invalid_argument("Problem not valid: " + msg);
        return;
    }

    // create dummy empty hybzono constraint
    X_hz.set(Eigen::SparseMatrix<double>(0, 0), Eigen::SparseMatrix<double>(0, 0), Eigen::VectorXd(0), 
        Eigen::SparseMatrix<double>(0, 0), Eigen::SparseMatrix<double>(0, 0), Eigen::VectorXd(0));

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
}

// get problem dimensions
void MPC_MIQP_Hrep::get_problem_dimensions()
{
    // call base class method
    MPC_MIQP::get_problem_dimensions();

    // binvars are in wrong spot so need to redefine
    int idx_offset = 0;
    for (int k=0; k<=n_horizon; k++)
    {
        idx_binvar[k] = get_indices(idx_offset + n_y_vec_orig[k], n_regions);
        idx_offset += n_y_vec[k];
    }
}

// make inequality constraints
void MPC_MIQP_Hrep::make_inequality_constraints()
{
    // declare
    std::vector<Eigen::Triplet<double>> tripvec;
    int i_offset, j_offset;

    // reset inequality constraints to original problem
    G_i_nom = G_i_nom_orig;
    G_i_vec = G_i_vec_orig;
    w_i_nom = w_i_nom_orig;
    w_i_vec = w_i_vec_orig;

    // concatenated constraint matrices
    tripvec.clear();
    i_offset = 0;
    for (auto const& [key, val] : A_Hrep)
    {
        Eigen::SparseMatrix<double> A_sp = val.sparseView();
        get_triplets_offset(A_sp*Pp, tripvec, i_offset, 0);
        i_offset += val.rows();
    }
    Eigen::SparseMatrix<double> A_cat (n_cons, n_state);
    A_cat.setFromTriplets(tripvec.begin(), tripvec.end());

    Eigen::VectorXd b_cat (n_cons);
    i_offset = 0;
    for (auto const& [key, val] : b_Hrep)
    {
        b_cat.segment(i_offset, val.size()) = val + M*Eigen::VectorXd::Ones(val.size());
        i_offset += val.size();
    }

    // selector matrices with big M
    tripvec.clear();
    i_offset = 0;
    for (auto const& [key, val] : A_Hrep)
    {
        for (int i=0; i<val.rows(); i++)
        {
            tripvec.push_back(Eigen::Triplet<double>(i_offset+i, key-1, M)); // one-based to zero-based indexing
        }
        i_offset += val.rows();
    }
    Eigen::SparseMatrix<double> M_s (n_cons, n_regions);
    M_s.setFromTriplets(tripvec.begin(), tripvec.end());

    // lambda functions to update G_i, w_i
    auto update_G_w_soft = [&](Eigen::SparseMatrix<double> &G_i, Eigen::VectorXd &w_i) 
    {
        tripvec.clear();
        i_offset = 0;

        // update G_i
        get_triplets_offset(G_i, tripvec, 0, 0);
        get_triplets_offset(A_cat, tripvec, G_i.rows(), 0);
        get_triplets_offset(M_s, tripvec, G_i.rows(), G_i.cols());
        get_triplets_offset_neg_identity(n_cons, tripvec, G_i.rows(), G_i.cols() + n_regions);
        get_triplets_offset_neg_identity(n_regions, tripvec, G_i.rows() + n_cons, G_i.cols());
        get_triplets_offset_identity(n_regions, tripvec, G_i.rows() + n_cons + n_regions, G_i.cols());
        get_triplets_offset_neg_identity(n_cons, tripvec, G_i.rows()+n_cons+2*n_regions, G_i.cols() + n_regions);
        get_triplets_offset_identity(n_cons, tripvec, G_i.rows()+n_cons+2*n_regions+n_cons, G_i.cols() + n_regions);

        G_i.resize(G_i.rows() + n_cons + 2*n_regions + 2*n_cons, G_i.cols() + n_regions + n_cons);
        G_i.setFromTriplets(tripvec.begin(), tripvec.end());

        // update w_i
        int w_i_size_orig = w_i.size();
        w_i.conservativeResize(G_i.rows());
        w_i.segment(w_i_size_orig, n_cons) = b_cat;
        w_i.segment(w_i_size_orig + n_cons, n_regions) = Eigen::VectorXd::Zero(n_regions);
        w_i.segment(w_i_size_orig + n_cons + n_regions, n_regions) = Eigen::VectorXd::Ones(n_regions);
        w_i.segment(w_i_size_orig + n_cons + 2*n_regions, n_cons) = Eigen::VectorXd::Zero(n_cons);
        w_i.segment(w_i_size_orig + n_cons + 2*n_regions + n_cons, n_cons) = sigma_max_hr*Eigen::VectorXd::Ones(n_cons);
    };

    auto update_G_w_hard = [&](Eigen::SparseMatrix<double> &G_i, Eigen::VectorXd &w_i) 
    {
        tripvec.clear();
        i_offset = 0;

        // update G_i
        get_triplets_offset(G_i, tripvec, 0, 0);
        get_triplets_offset(A_cat, tripvec, G_i.rows(), 0);
        get_triplets_offset(M_s, tripvec, G_i.rows(), G_i.cols());
        get_triplets_offset_neg_identity(n_regions, tripvec, G_i.rows() + n_cons, G_i.cols());
        get_triplets_offset_identity(n_regions, tripvec, G_i.rows() + n_cons + n_regions, G_i.cols());

        G_i.resize(G_i.rows() + n_cons + 2*n_regions, G_i.cols() + n_regions);
        G_i.setFromTriplets(tripvec.begin(), tripvec.end());

        // update w_i
        int w_i_size_orig = w_i.size();
        w_i.conservativeResize(G_i.rows());
        w_i.segment(w_i_size_orig, n_cons) = b_cat;
        w_i.segment(w_i_size_orig + n_cons, n_regions) = Eigen::VectorXd::Zero(n_regions);
        w_i.segment(w_i_size_orig + n_cons + n_regions, n_regions) = Eigen::VectorXd::Ones(n_regions);
    };

    // nominal inequality constraints
    if (softened_hrep_constraints)
        update_G_w_soft(G_i_nom, w_i_nom);
    else
        update_G_w_hard(G_i_nom, w_i_nom);

    // off-nominal constraints
    for (int k=0; k<=n_horizon; k++)
    {
        int G_cnt = G_i_vec.count(k);
        int w_cnt = w_i_vec.count(k);

        if (G_cnt || w_cnt)
        {
            if (!G_cnt)
                G_i_vec[k] = G_i_nom_orig;;
            if (!w_cnt)
                w_i_vec[k] = w_i_nom_orig;

            if (softened_hrep_constraints)
                update_G_w_soft(G_i_vec[k], w_i_vec[k]);
            else
                update_G_w_hard(G_i_vec[k], w_i_vec[k]);
        }
    }
}

// make equality constraints
void MPC_MIQP_Hrep::make_equality_constraints()
{
    // declare
    std::vector<Eigen::Triplet<double>> tripvec;
    int i_offset, j_offset;

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

    // lambda functions to update C_i, D_i, crhs_i
    auto update_C = [&](Eigen::SparseMatrix<double> &C_i, bool soft) 
    {
        tripvec.clear();

        // update C_i
        get_triplets_offset(C_i, tripvec, 0, 0);
        for (int i=0; i<n_regions; i++)
        {
            tripvec.push_back(Eigen::Triplet<double>(C_i.rows(), C_i.cols() + i, 1));
        }

        if (soft)
            C_i.resize(C_i.rows()+1, C_i.cols()+n_regions+n_cons);
        else
            C_i.resize(C_i.rows()+1, C_i.cols()+n_regions);

        C_i.setFromTriplets(tripvec.begin(), tripvec.end());
    };
    
    auto update_D = [&](Eigen::SparseMatrix<double> &D_i, bool soft)
    {
        // update D_i
        if (soft)
            D_i.conservativeResize(D_i.rows()+1, D_i.cols()+n_regions+n_cons);
        else
            D_i.conservativeResize(D_i.rows()+1, D_i.cols()+n_regions);
    };

    auto update_crhs = [&](Eigen::VectorXd &crhs_i)
    {
        // update crhs_i
        crhs_i.conservativeResize(crhs_i.size()+1);
        crhs_i(crhs_i.size()-1) = -1;
    };

    // nominal equality constraints
    update_C(C_i_nom, softened_hrep_constraints);
    update_D(D_i_nom, softened_hrep_constraints);
    update_crhs(crhs_i_nom);
    idx_i_nom = get_indices(C_i_nom_orig.cols(), n_regions);

    // off-nominal constraints
    for (int k=0; k<=n_horizon; k++)
    {
        if (C_i_vec.count(k))
        {
            int n_C_orig = C_i_vec[k].cols();
            update_C(C_i_vec[k], softened_hrep_constraints);
            idx_i_vec[k] = get_indices(n_C_orig, n_regions);
        }
        if (D_i_vec.count(k))
        {
            int n_D_orig = D_i_vec[k].cols();
            update_D(D_i_vec[k], softened_hrep_constraints);
            idx_i_vec[k] = get_indices(n_D_orig, n_regions);
        }
        if (crhs_i_vec.count(k))
        {
            update_crhs(crhs_i_vec[k]);
        }
    }

    // terminal H-rep constraints
    if (!D_i_vec.count(n_horizon))
        D_i_vec[n_horizon] = D_i_nom;

    // D_N
    tripvec.clear();
    get_triplets_offset(D_i_vec[n_horizon], tripvec, 0, 0);
    for (int i=0; i<n_regions; i++)
    {
        tripvec.push_back(Eigen::Triplet<double>(D_i_vec[n_horizon].rows(), D_i_vec_orig[n_horizon].cols() + i, 1));
    }
    D_i_vec[n_horizon].resize(D_i_vec[n_horizon].rows()+1, D_i_vec[n_horizon].cols());
    D_i_vec[n_horizon].setFromTriplets(tripvec.begin(), tripvec.end());

    // C_Nm1
    if (!C_i_vec.count(n_horizon-1))
        C_i_vec[n_horizon-1] = C_i_nom;

    C_i_vec[n_horizon-1].conservativeResize(C_i_vec[n_horizon-1].rows()+1, C_i_vec[n_horizon-1].cols());

    // crhs_Nm1
    if (!crhs_i_vec.count(n_horizon-1))
        crhs_i_vec[n_horizon-1] = crhs_i_nom;
    
    crhs_i_vec[n_horizon-1].conservativeResize(crhs_i_vec[n_horizon-1].size()+1);
    crhs_i_vec[n_horizon-1](crhs_i_vec[n_horizon-1].size()-1) = -1;
}

// make cost function
void MPC_MIQP_Hrep::make_cost_function()
{
    // declare
    std::vector<Eigen::Triplet<double>> tripvec;
    int ij_offset;

    // reset cost function matrices to original problem
    P_i_nom = P_i_nom_orig;
    P_i_vec = P_i_vec_orig;
    q_i_nom = q_i_nom_orig;
    q_i_vec = q_i_vec_orig;

    // lambdas to update P, q
    auto update_P = [&] (Eigen::SparseMatrix<double> &P_i, bool soft)
    {
        tripvec.clear();

        if (soft)
        {
            get_triplets_offset(P_i, tripvec, 0, 0);
            ij_offset = P_i.rows()+n_regions;
            get_triplets_offset(Q_hr, tripvec, ij_offset, ij_offset);
            ij_offset += n_cons;
        }
        else
        {
            get_triplets_offset(P_i, tripvec, 0, 0);
            ij_offset = P_i.rows()+n_regions;
        }
        
        P_i.resize(ij_offset, ij_offset);
        P_i.setFromTriplets(tripvec.begin(), tripvec.end());
    };

    auto update_q = [&] (Eigen::VectorXd &q_i, bool soft, bool region_costs)
    {
        if (soft)
        {
            ij_offset = q_i.size()+n_regions+n_cons;    
        }
        else
        {
            ij_offset = q_i.size()+n_regions;   
        }

        int q_i_size_init = q_i.size();
        q_i.conservativeResize(ij_offset);
        q_i.segment(q_i_size_init, ij_offset-q_i_size_init) = Eigen::VectorXd::Zero(ij_offset-q_i_size_init);
        
        if (region_costs)
        {
            q_i.segment(q_i_size_init, n_regions) = region_cost_vec;
        }
    };

    // nominal cost
    update_P(P_i_nom, softened_hrep_constraints);
    update_q(q_i_nom, softened_hrep_constraints, region_costs);

    // off-nominal cost
    for (int k=0; k<=n_horizon; k++)
    {
        if (P_i_vec.count(k))
        {
            update_P(P_i_vec[k], softened_hrep_constraints);
        }
        if (q_i_vec.count(k))
        {
            update_q(q_i_vec[k], softened_hrep_constraints, region_costs);
        }
    }
}

// unsupported inherited methods
void MPC_MIQP_Hrep::set_terminal_hybzono_constraints(const ZonoCpp::HybZono &X_hz_term, const Eigen::Ref<const Eigen::VectorXd> term_region_cost_vec,
            const Eigen::SparseMatrix<double> &Pp_term, const Eigen::SparseMatrix<double> &Q_term_slack_cost)
{
    throw std::runtime_error("Terminal hybzono constraints not supported for H-rep constraints");
}

int MPC_MIQP_Hrep::get_terminal_region_selection()
{
    throw std::runtime_error("Terminal hybzono constraints not supported for H-rep constraints");
}