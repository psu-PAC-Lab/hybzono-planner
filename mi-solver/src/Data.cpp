#include "Data.hpp"

using namespace MI_QPIPMPC;

// constructors
Data::Data(const QPIPMPC_Data &qpipmpc_data, 
            const std::map<int, std::vector<int>> &R_x02region, 
            const std::map<int, std::map<int, std::vector<int>>> &R_region2region,
            int region_0, const MI_QPIPMPC_Settings &settings, 
            const QP_IP_MPC::QP_settings &qp_settings,
            const std::map<int, Eigen::MatrixXd> &A_Hrep,
            const std::map<int, Eigen::VectorXd> &b_Hrep): 
            qpipmpc_data(qpipmpc_data), R_x02region(R_x02region),
            R_region2region(R_region2region), region_0(region_0), 
            settings(settings), qp_settings(qp_settings), A_Hrep(A_Hrep), b_Hrep(b_Hrep)
{
    n_regions = this->R_region2region[1].size();
    n_horizon = this->qpipmpc_data.n_horizon;
    Hrep_defined = true;
    get_problem_dimension();
    get_pos_vars();

    solver.setup(this->qpipmpc_data.P_i_nom, this->qpipmpc_data.q_i_nom,
        this->qpipmpc_data.C_i_nom, this->qpipmpc_data.D_i_nom,
        this->qpipmpc_data.crhs_i_nom, this->qpipmpc_data.G_i_nom,
        this->qpipmpc_data.w_i_nom, this->n_horizon);
    solver.set_P_vec(this->qpipmpc_data.P_i_vec);
    solver.set_q_vec(this->qpipmpc_data.q_i_vec);
    solver.set_C_vec(this->qpipmpc_data.C_i_vec);
    solver.set_D_vec(this->qpipmpc_data.D_i_vec);
    solver.set_crhs_vec(this->qpipmpc_data.crhs_i_vec);
    solver.set_G_vec(this->qpipmpc_data.G_i_vec);
    solver.set_w_vec(this->qpipmpc_data.w_i_vec);
    solver.set_settings(this->qp_settings);
}

Data::Data(const QPIPMPC_Data &qpipmpc_data, 
            const std::map<int, std::vector<int>> &R_x02region, 
            const std::map<int, std::map<int, std::vector<int>>> &R_region2region,
            int region_0, const MI_QPIPMPC_Settings &settings, 
            const QP_IP_MPC::QP_settings &qp_settings): 
            qpipmpc_data(qpipmpc_data), R_x02region(R_x02region),
            R_region2region(R_region2region), region_0(region_0), 
            settings(settings), qp_settings(qp_settings)
{
    n_regions = this->R_region2region[1].size();
    n_horizon = this->qpipmpc_data.n_horizon;
    Hrep_defined = false;
    get_problem_dimension();
    get_pos_vars();

    solver.setup(this->qpipmpc_data.P_i_nom, this->qpipmpc_data.q_i_nom,
        this->qpipmpc_data.C_i_nom, this->qpipmpc_data.D_i_nom,
        this->qpipmpc_data.crhs_i_nom, this->qpipmpc_data.G_i_nom,
        this->qpipmpc_data.w_i_nom, this->n_horizon);
    solver.set_P_vec(this->qpipmpc_data.P_i_vec);
    solver.set_q_vec(this->qpipmpc_data.q_i_vec);
    solver.set_C_vec(this->qpipmpc_data.C_i_vec);
    solver.set_D_vec(this->qpipmpc_data.D_i_vec);
    solver.set_crhs_vec(this->qpipmpc_data.crhs_i_vec);
    solver.set_G_vec(this->qpipmpc_data.G_i_vec);
    solver.set_w_vec(this->qpipmpc_data.w_i_vec);
    solver.set_settings(this->qp_settings);
}

Data::Data()
{
    n_regions = 0;
    region_0 = 0;
    n_horizon = 0;
    n_term_regions = 0;
}

Data::Data(Data &&other)
{
    qpipmpc_data = std::move(other.qpipmpc_data);
    qp_settings = std::move(other.qp_settings);
    R_x02region = std::move(other.R_x02region);
    R_region2region = std::move(other.R_region2region);
    region_0 = other.region_0;
    n_regions = other.n_regions;
    n = other.n;
    Hrep_defined = other.Hrep_defined;
    A_Hrep = std::move(other.A_Hrep);
    b_Hrep = std::move(other.b_Hrep);
    settings = std::move(other.settings);
    n_horizon = other.n_horizon;
    solver = std::move(other.solver);
    terminal_constraint = other.terminal_constraint;
    C_term = std::move(other.C_term);
    D_term = std::move(other.D_term);
    G_term = std::move(other.G_term);
    P_term = std::move(other.P_term);
    crhs_term = std::move(other.crhs_term);
    w_term = std::move(other.w_term);
    q_term = std::move(other.q_term);
    idx_term = std::move(other.idx_term);
    idx_binvar_term = std::move(other.idx_binvar_term);
    n_term_regions = other.n_term_regions;
    warm_start = other.warm_start;
    region_vec_warm_start = std::move(other.region_vec_warm_start);
    idx_pos_vars = std::move(other.idx_pos_vars);
}

Data::Data(const Data &other)
{
    qpipmpc_data = other.qpipmpc_data;
    qp_settings = other.qp_settings;
    R_x02region = other.R_x02region;
    R_region2region = other.R_region2region;
    region_0 = other.region_0;
    n_regions = other.n_regions;
    n = other.n;
    Hrep_defined = other.Hrep_defined;
    A_Hrep = other.A_Hrep;
    b_Hrep = other.b_Hrep;
    settings = other.settings;
    n_horizon = other.n_horizon;
    solver = other.solver;
    terminal_constraint = other.terminal_constraint;
    C_term = other.C_term;
    D_term = other.D_term;
    G_term = other.G_term;
    P_term = other.P_term;
    crhs_term = other.crhs_term;
    w_term = other.w_term;
    q_term = other.q_term;
    idx_term = other.idx_term;
    idx_binvar_term = other.idx_binvar_term;
    n_term_regions = other.n_term_regions;
    warm_start = other.warm_start;
    region_vec_warm_start = other.region_vec_warm_start;
    idx_pos_vars = other.idx_pos_vars;
}

// copy assignment
Data& Data::operator=(const Data &other)
{
    if (this != &other)
    {
        qpipmpc_data = other.qpipmpc_data;
        qp_settings = other.qp_settings;
        R_x02region = other.R_x02region;
        R_region2region = other.R_region2region;
        region_0 = other.region_0;
        n_regions = other.n_regions;
        n = other.n;
        Hrep_defined = other.Hrep_defined;
        A_Hrep = other.A_Hrep;
        b_Hrep = other.b_Hrep;
        settings = other.settings;
        n_horizon = other.n_horizon;    
        solver = other.solver;
        terminal_constraint = other.terminal_constraint;
        C_term = other.C_term;
        D_term = other.D_term;
        G_term = other.G_term;
        P_term = other.P_term;
        crhs_term = other.crhs_term;
        w_term = other.w_term;
        q_term = other.q_term;
        idx_term = other.idx_term;
        idx_binvar_term = other.idx_binvar_term;
        n_term_regions = other.n_term_regions;
        warm_start = other.warm_start;
        region_vec_warm_start = other.region_vec_warm_start;
        idx_pos_vars = other.idx_pos_vars;
    }
    return *this;
}

// get problem dimension
void Data::get_problem_dimension()
{
    n = 0; // init
    for (auto it = qpipmpc_data.idx_y.begin(); it != qpipmpc_data.idx_y.end(); ++it)
    {
        n += it->second.size();
    }
}

// update/set methods
void Data::set_qp_settings(QP_IP_MPC::QP_settings &qp_settings_in)
{
    qp_settings = qp_settings;
}

void Data::set_Hrep(const std::map<int, Eigen::MatrixXd> &A_Hrep,
            const std::map<int, Eigen::VectorXd> &b_Hrep)
{
    this->A_Hrep = A_Hrep;
    this->b_Hrep = b_Hrep;
    this->n_regions = A_Hrep.size();
}

void Data::set_terminal_constraint(const Eigen::SparseMatrix<double> &C_term,
            const Eigen::SparseMatrix<double> &D_term,
            const Eigen::SparseMatrix<double> &G_term,
            const Eigen::SparseMatrix<double> &P_term,
            const Eigen::Ref<const Eigen::VectorXd> crhs_term,
            const Eigen::Ref<const Eigen::VectorXd> w_term,
            const Eigen::Ref<const Eigen::VectorXd> q_term,
            const std::vector<int> &idx_term,
            const std::vector<int> &idx_binvar_term)
{
    this->terminal_constraint = true;
    this->C_term = C_term;
    this->D_term = D_term;
    this->G_term = G_term;
    this->P_term = P_term;
    this->crhs_term = crhs_term;
    this->w_term = w_term;
    this->q_term = q_term;
    this->idx_term = idx_term;
    this->idx_binvar_term = idx_binvar_term;
    this->n_term_regions = idx_term.size();
    this->apply_terminal_constraints_to_problem_data();
}

void Data::update_reachability(const std::map<int, std::vector<int>> &R_x02region,
            const std::map<int, std::map<int, std::vector<int>>> &R_region2region,
            int region_0)
{
    this->R_x02region = R_x02region;
    this->R_region2region = R_region2region;
    this->region_0 = region_0;
    this->n_horizon = R_x02region.size();
}

void Data::set_P_i_nom(const Eigen::SparseMatrix<double> &P_i_nom)
{
    qpipmpc_data.P_i_nom = P_i_nom;
}

void Data::set_q_i_nom(const Eigen::Ref<const Eigen::VectorXd> q_i_nom)
{
    qpipmpc_data.q_i_nom = q_i_nom;
}

void Data::set_b(double b)
{
    qpipmpc_data.b = b;
}   

void Data::set_P_i_vec(const std::map<int, Eigen::SparseMatrix<double>> &P_i_vec)
{
    qpipmpc_data.P_i_vec = P_i_vec;
    if (terminal_constraint)
        qpipmpc_data.P_i_vec[n_horizon] = P_term;
}

void Data::set_q_i_vec(const std::map<int, Eigen::VectorXd> &q_i_vec)
{
    qpipmpc_data.q_i_vec = q_i_vec;
    if (terminal_constraint)
        qpipmpc_data.q_i_vec[n_horizon] = q_term;
}

void Data::set_C_i_nom(const Eigen::SparseMatrix<double> &C_i_nom)
{
    qpipmpc_data.C_i_nom = C_i_nom;
}        

void Data::set_D_i_nom(const Eigen::SparseMatrix<double> &D_i_nom)
{
    qpipmpc_data.D_i_nom = D_i_nom;
}

void Data::set_crhs_i_nom(const Eigen::Ref<const Eigen::VectorXd> crhs_i_nom)
{
    qpipmpc_data.crhs_i_nom = crhs_i_nom;
}

void Data::set_C_i_vec(const std::map<int, Eigen::SparseMatrix<double>> &C_i_vec)
{
    qpipmpc_data.C_i_vec = C_i_vec;
    if (terminal_constraint)
        qpipmpc_data.C_i_vec[n_horizon-1] = C_term;
}

void Data::set_D_i_vec(const std::map<int, Eigen::SparseMatrix<double>> &D_i_vec)
{
    qpipmpc_data.D_i_vec = D_i_vec;
    if (terminal_constraint)
        qpipmpc_data.D_i_vec[n_horizon] = D_term;
}

void Data::set_crhs_i_vec(const std::map<int, Eigen::VectorXd> crhs_i_vec)
{
    qpipmpc_data.crhs_i_vec = crhs_i_vec;
    if (terminal_constraint)
        qpipmpc_data.crhs_i_vec[n_horizon-1] = crhs_term;
}

void Data::set_G_i_nom(const Eigen::SparseMatrix<double> &G_i_nom)
{
    qpipmpc_data.G_i_nom = G_i_nom; 
    throw std::runtime_error("TO DO: Need to implement a problem update method and call get_pos_vars() only once");
}

void Data::set_w_i_nom(const Eigen::Ref<const Eigen::VectorXd> w_i_nom)
{
    qpipmpc_data.w_i_nom = w_i_nom;
    throw std::runtime_error("TO DO: Need to implement a problem update method and call get_pos_vars() only once");
}

void Data::set_G_i_vec(const std::map<int, Eigen::SparseMatrix<double>> &G_i_vec)
{
    qpipmpc_data.G_i_vec = G_i_vec;
    if (terminal_constraint)
        qpipmpc_data.G_i_vec[n_horizon] = G_term;
    throw std::runtime_error("TO DO: Need to implement a problem update method and call get_pos_vars() only once");
}

void Data::set_w_i_vec(const std::map<int, Eigen::VectorXd> &w_i_vec)
{
    qpipmpc_data.w_i_vec = w_i_vec;
    if (terminal_constraint)
        qpipmpc_data.w_i_vec[n_horizon] = w_term;
    throw std::runtime_error("TO DO: Need to implement a problem update method and call get_pos_vars() only once");
}

void Data::set_idx_i_nom(const std::vector<int> &idx_i_nom)
{
    qpipmpc_data.idx_i_nom = idx_i_nom;
}

void Data::set_idx_i_vec(const std::map<int, std::vector<int>> &idx_i_vec)
{
    qpipmpc_data.idx_i_vec = idx_i_vec;
}

void Data::set_idx_state(const std::map<int, std::vector<int>> &idx_state)
{
    qpipmpc_data.idx_state = idx_state;
}

void Data::set_idx_input(const std::map<int, std::vector<int>> &idx_input)
{
    qpipmpc_data.idx_input = idx_input;
}

void Data::set_idx_binvar(const std::map<int, std::vector<int>> &idx_binvar)
{
    qpipmpc_data.idx_binvar = idx_binvar;
}

void Data::set_idx_y(const std::map<int, std::vector<int>> &idx_y)
{
    qpipmpc_data.idx_y = idx_y;
    get_problem_dimension();
}

void Data::set_idx_eq(const std::map<int, std::vector<int>> &idx_eq)
{
    qpipmpc_data.idx_eq = idx_eq;
}

void Data::set_idx_ineq(const std::map<int, std::vector<int>> &idx_ineq)
{
    qpipmpc_data.idx_ineq = idx_ineq;
}

// terminal constraints
void Data::apply_terminal_constraints_to_problem_data()
{
    if (terminal_constraint)
    {
        qpipmpc_data.P_i_vec[n_horizon] = P_term;
        qpipmpc_data.q_i_vec[n_horizon] = q_term;
        qpipmpc_data.C_i_vec[n_horizon-1] = C_term;
        qpipmpc_data.D_i_vec[n_horizon] = D_term;
        qpipmpc_data.crhs_i_vec[n_horizon-1] = crhs_term;
        qpipmpc_data.G_i_vec[n_horizon] = G_term;
        qpipmpc_data.w_i_vec[n_horizon] = w_term;
    }
}

// warm start regions
void Data::set_warm_start_regions(const std::vector<std::vector<int>> &region_vec_warm_start, int term_region_warm_start)
{
    this->region_vec_warm_start = region_vec_warm_start;
    this->warm_start = true;
    this->term_region_warm_start = term_region_warm_start;
}

// get indices of variables that are guaranteed positive
void Data::get_pos_vars()
{
    // initializations / declarations
    Eigen::SparseMatrix<double> * G_ptr = nullptr; // init
    Eigen::VectorXd * w_ptr = nullptr; // init
    std::vector<int> elem_in_row, zero_var_k; // declare

    // loop through by time step
    for (int k=0; k<=n_horizon; k++)
    {
        // get inequality constraints
        G_ptr = qpipmpc_data.G_i_vec.count(k) ? &qpipmpc_data.G_i_vec[k] : &qpipmpc_data.G_i_nom;
        w_ptr = qpipmpc_data.w_i_vec.count(k) ? &qpipmpc_data.w_i_vec[k] : &qpipmpc_data.w_i_nom;

        // reset
        elem_in_row.clear();
        zero_var_k.clear();
        for (int i=0; i<w_ptr->size(); i++)
        {
            elem_in_row.push_back(0);
            zero_var_k.push_back(-1);
        }

        // find positive variables
        for (int i=0; i<G_ptr->outerSize(); i++)
        {
            for (Eigen::SparseMatrix<double>::InnerIterator it(*G_ptr, i); it; ++it)
            {
                elem_in_row[it.row()]++;
                if ((*w_ptr)[it.row()] == 0)
                {
                    if (it.value() < 0)
                        zero_var_k[it.row()] = it.col();
                    else
                        zero_var_k[it.row()] = -1; // invalid
                }
            }
        }

        // update idx_pos_vars
        for (int i=0; i<w_ptr->size(); i++)
        {
            if (elem_in_row[i] == 1 && zero_var_k[i] != -1 && (*w_ptr)[i] == 0)
                idx_pos_vars[k].push_back(zero_var_k[i]);
        }
    }
}