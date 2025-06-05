#include "QP_primal_dual_mpc.hpp"

using namespace QP_IP_MPC;

// constructor
QP_primal_dual_mpc::QP_primal_dual_mpc() {}

// destructor
QP_primal_dual_mpc::~QP_primal_dual_mpc() = default;

// set settings
void QP_primal_dual_mpc::set_settings(const QP_settings &settings_in)
{
    settings = settings_in;
}

// setup
void QP_primal_dual_mpc::setup(const Eigen::SparseMatrix<double> &P_i_nom_in, const Eigen::VectorXd &q_i_nom_in, 
            const Eigen::SparseMatrix<double> &C_i_nom_in, const Eigen::SparseMatrix<double> &D_i_nom_in,
            const Eigen::VectorXd &crhs_i_nom_in, const Eigen::SparseMatrix<double> &G_i_nom_in,
            const Eigen::VectorXd &w_i_nom_in, int n_horizon_in, double b)
{
    // copy in matrices / vectors
    P_i_nom = P_i_nom_in;
    q_i_nom = q_i_nom_in;
    C_i_nom = C_i_nom_in;
    D_i_nom = D_i_nom_in;
    crhs_i_nom = crhs_i_nom_in;
    G_i_nom = G_i_nom_in;
    w_i_nom = w_i_nom_in;
    n_horizon = n_horizon_in;

    // optional arguments
    this->b = b;

    // check if inequality constraints are box constraints
    isBox_i_nom = checkIfBox(&G_i_nom);

    // check if quadratic cost is diagonal
    isDiag_i_nom = checkIfDiag(&P_i_nom);

    // update flags
    P_updated = true;
    q_updated = true; 
    C_updated = true; 
    crhs_updated = true;
    G_updated = true;
    w_updated = true;
    b_updated = true;
}

void QP_primal_dual_mpc::set_P_vec(const std::map<int, Eigen::SparseMatrix<double>> &P_vec_in)
{
    P_i_vec = P_vec_in;
    for (auto it = P_i_vec.begin(); it != P_i_vec.end(); it++)
        isDiag_i_vec[it->first] = checkIfDiag(&(it->second));
    P_updated = true;
}

void QP_primal_dual_mpc::set_q_vec(const std::map<int, Eigen::VectorXd> &q_vec_in)
{
    q_i_vec = q_vec_in;
    q_updated = true;
}

void QP_primal_dual_mpc::set_b(double b)
{
    this->b = b;
    b_updated = true;
}

void QP_primal_dual_mpc::set_C_vec(const std::map<int, Eigen::SparseMatrix<double>> &C_vec_in)
{
    C_i_vec = C_vec_in;
    C_updated = true;
}

void QP_primal_dual_mpc::set_D_vec(const std::map<int, Eigen::SparseMatrix<double>> &D_vec_in)
{
    D_i_vec = D_vec_in;
    C_updated = true;
}

void QP_primal_dual_mpc::set_crhs_vec(const std::map<int, Eigen::VectorXd> &crhs_vec_in)
{
    crhs_i_vec = crhs_vec_in;
    crhs_updated = true;
}

void QP_primal_dual_mpc::set_G_vec(const std::map<int, Eigen::SparseMatrix<double>> &G_vec_in)
{
    G_i_vec = G_vec_in;
    for (auto it = G_i_vec.begin(); it != G_i_vec.end(); it++)
        isBox_i_vec[it->first] = checkIfBox(&(it->second));
    G_updated = true;
}

void QP_primal_dual_mpc::set_w_vec(const std::map<int, Eigen::VectorXd> &w_vec_in)
{
    w_i_vec = w_vec_in;
    w_updated = true;
}

// solve method
QP_results QP_primal_dual_mpc::solve()
{
    // declare
    double mu, h, mu_pred, sigma;
    Eigen::VectorXd s_pred, lambda_pred;
    bool numerical_limit, feas;
    std::pair<double, bool> line_search_out;

    // start timer
    auto timer_init = std::chrono::high_resolution_clock::now();

    // running timer
    bool enforce_max_time = (settings.T_max > 0.0);
    double running_timer = 0; // init
    double T_max; // declare
    if (enforce_max_time)
        T_max = settings.T_max;
    else        
        T_max = inf;

    // compute problem dimensions
    computeProblemDimensions();

    // construct full system matrices
    if (P_updated) build_P();
    if (q_updated) build_q();
    if (C_updated) build_C();
    if (crhs_updated) build_crhs();
    if (G_updated) build_G();
    if (w_updated) build_w();

    // initialize working matrices
    initializeWorkingMatrices();

    // initialize start var if necessary
    y0 = Eigen::VectorXd::Zero(n_y_tot);
    double zeta = sqrt(settings.mu_init);
    nu0 = Eigen::VectorXd::Zero(n_nu_tot);
    lambda0 = zeta*Eigen::VectorXd::Ones(n_lambda_s_tot);
    s0 = zeta*Eigen::VectorXd::Ones(n_lambda_s_tot);

    // initialize primal and dual vars
    y = y0;
    nu = nu0;
    lambda = lambda0;
    s = s0;

    // make s_inv_lambda
    numerical_limit = make_S_inv_Lambda();

    // outer loop init
    int k = 0;

    // compute duality measure
    mu = computeMu(s, lambda);

    // loop
    while ((mu > settings.mu_term) && (k < settings.iter_max) && !numerical_limit && (running_timer < T_max) && (mu <= settings.mu_max))
    {       
        // generate system matrix and decompose
        factorizeSystemMatrix();

        // predictor step
        generateRHS();
        solveYdelnu_eq_beta();
        computeFullStepDirectionFrom_del_nu();

        // line search for step size
        line_search_out = lineSearch(del_s, del_lambda);
        h = line_search_out.first;
        if (line_search_out.second)
        {
            numerical_limit = true;
            break;
        }

        // predicted duality measure
        s_pred = s + h*del_s;
        lambda_pred = lambda + h*del_lambda;
        mu_pred = computeMu(s_pred, lambda_pred);

        // centering parameter
        sigma = pow(mu_pred/mu, 3);

        // calculate corrected nu and recompute search direction
        Del_s.diagonal() = del_s;
        nu_centering = sigma*mu*Eigen::VectorXd::Ones(n_lambda_s_tot) - Del_s*del_lambda;
        updateRHS();

        solveYdelnu_eq_beta();
        computeFullStepDirectionFrom_del_nu();

        // line search for step size
        line_search_out = lineSearch(del_s, del_lambda);
        h = line_search_out.first;
        if (line_search_out.second)
        {
            numerical_limit = true;
            break;
        }

        // update
        y += settings.gamma*h*del_y;
        nu += settings.gamma*h*del_nu;
        lambda += settings.gamma*h*del_lambda;
        s += settings.gamma*h*del_s;

        // compute duality measure
        mu = computeMu(s, lambda);

        // compute s_inv_lambda
        numerical_limit = make_S_inv_Lambda();

        // update running timer
        if (enforce_max_time)
        {
            auto timer_running = std::chrono::high_resolution_clock::now();
            auto duration_running = std::chrono::duration_cast<std::chrono::microseconds>(timer_running - timer_init);
            running_timer = 1e-6 * ((double) duration_running.count());
        }

        // iterate
        k++;
    }

    // timing
    auto timer_final = std::chrono::high_resolution_clock::now();
    auto duration_final = std::chrono::duration_cast<std::chrono::microseconds>(timer_final - timer_init);
    double time = 1e-6 * ((double) duration_final.count());

    // assemble results
    QP_results results;
    results.y = y;
    results.nu = nu;
    results.lambda = lambda;
    results.s = s;
    results.objective = objective(y);
    results.feas = mu <= settings.mu_feas; // divergence check
    results.converged = !numerical_limit && (k < settings.iter_max) && results.feas;
    results.num_iter = k;
    results.sol_time = time;

    return results;
}

// factorization
void QP_primal_dual_mpc::factorizeSystemMatrix()
{
    // loop
    for (int i=0; i<=n_horizon; i++)
    {
        // get matrices needed to construct phi_i
        Eigen::SparseMatrix<double> * P_i_ptr = nullptr; // declare
        Eigen::SparseMatrix<double> * G_i_ptr = nullptr; // declare
        Eigen::SparseMatrix<double> * G_i_ptr_T = nullptr; // declare
        bool P_diag, G_box;
        if (P_i_vec.count(i)) // step-specific value specified
        {
            P_i_ptr = &(P_i_vec[i]);
            P_diag = isDiag_i_vec[i];   
        }
        else
        {
            P_i_ptr = &(P_i_nom);
            P_diag = isDiag_i_nom;
        }
        if (G_i_vec.count(i)) // step-specific value specified
        {
            G_i_ptr = &(G_i_vec[i]);
            G_i_ptr_T = &(G_i_vec_T[i]);
            G_box = isBox_i_vec[i];
        }
        else
        {
            G_i_ptr = &(G_i_nom);
            G_i_ptr_T = &(G_i_nom_T);
            G_box = isBox_i_nom;
        }

        // construct phi_i
        Eigen::MatrixXd phi_i_dense; // declare
        Eigen::SparseMatrix<double> phi_i = *P_i_ptr + (*G_i_ptr_T)*S_inv_Lambda_vec[i]*(*G_i_ptr);

        // declare C_i and D_i pointers
        Eigen::SparseMatrix<double> * C_i_ptr = nullptr;
        Eigen::SparseMatrix<double> * D_i_ptr = nullptr;
        Eigen::SparseMatrix<double> * C_i_ptr_T = nullptr;
        Eigen::SparseMatrix<double> * D_i_ptr_T = nullptr;
        Eigen::MatrixXd * C_i_ptr_T_dense = nullptr;
        Eigen::MatrixXd * D_i_ptr_T_dense = nullptr;

        // skip cholesky decomposition if phi_i is diagonal
        isDiag_phi_i_vec[i] = P_diag && G_box;
        if (isDiag_phi_i_vec[i])
        {
            // get inverse of phi_i
            phi_i_inv_vec[i].diagonal() = phi_i.diagonal();
            phi_i_inv_vec[i] = phi_i_inv_vec[i].inverse();

            // get CPC
            if (i != n_horizon)
            {
                // get C_i
                if (C_i_vec.count(i)) // step-specific value specified
                {
                    C_i_ptr = &(C_i_vec[i]);
                    C_i_ptr_T = &(C_i_vec_T[i]);
                }
                else
                {
                    C_i_ptr = &(C_i_nom);
                    C_i_ptr_T = &(C_i_nom_T);
                }
                CPC[i] = Eigen::MatrixXd((*C_i_ptr)*phi_i_inv_vec[i]*(*C_i_ptr_T));
            }

            // get DPD
            if (i != 0)
            {
                // get D_i
                if (D_i_vec.count(i)) // step-specific value specified
                {
                    D_i_ptr = &(D_i_vec[i]);
                    D_i_ptr_T = &(D_i_vec_T[i]);
                }
                else
                {
                    D_i_ptr = &(D_i_nom);
                    D_i_ptr_T = &(D_i_nom_T);
                }
                DPD[i] = Eigen::MatrixXd((*D_i_ptr)*phi_i_inv_vec[i]*(*D_i_ptr_T));
            }

            // get DPC
            if ((i != n_horizon) && (i != 0))
                DPC[i] = Eigen::MatrixXd((*D_i_ptr)*phi_i_inv_vec[i]*(*C_i_ptr_T));
        }
        else // need to take Cholesky decomposition
        {
            // take cholesky decomposition
            phi_i_dense = Eigen::MatrixXd(phi_i);
            LLT_solver.compute(phi_i_dense);
            L_phi_i_vec[i] = LLT_solver.matrixL();

            // get C*Phi_inv*C', etc. via intermediate V, W matrices
            Eigen::MatrixXd V_i, W_i; // declare
            if (i != n_horizon) // V not defined at end of horizon
            {
                // get C_i
                if (C_i_vec.count(i)) // step-specific value specified
                    C_i_ptr_T_dense = &(C_i_vec_T_dense[i]);
                else
                    C_i_ptr_T_dense = &(C_i_nom_T_dense);

                // compute CPC by way of V_i
                V_i = (L_phi_i_vec[i].triangularView<Eigen::Lower>().solve(*C_i_ptr_T_dense)).transpose();
                CPC[i] = V_i*V_i.transpose();
            }

            if (i != 0)
            {
                // get D_i
                if (D_i_vec.count(i)) // step-specific value specified
                    D_i_ptr_T_dense = &(D_i_vec_T_dense[i]);
                else
                    D_i_ptr_T_dense = &(D_i_nom_T_dense);

                // compute DPD by way of W_i
                W_i = (L_phi_i_vec[i].triangularView<Eigen::Lower>().solve(*D_i_ptr_T_dense)).transpose();
                DPD[i] = W_i*W_i.transpose();
            }

            // compute DPC by way of V_i and W_i
            if ((i != n_horizon) && (i != 0))
                DPC[i] = W_i*V_i.transpose();
        
        } // end else: need to take Cholesky decomposition

        // got everything we need for step 0, continue
        if (i == 0)
            continue;
        
        // Cholesky factorization
        if (i == 1)
            LLT_solver.compute(CPC[i-1] + DPD[i]);
        else
        {
            Eigen::MatrixXd Y_mat = CPC[i-1] + DPD[i] - L_ip1_i_vec[i-1]*L_ip1_i_vec[i-1].transpose();
            LLT_solver.compute(Y_mat);
        }
        L_i_i_vec[i] = LLT_solver.matrixL();

        // forward substitution
        if (i < n_horizon)
            L_ip1_i_vec[i] = (L_i_i_vec[i].triangularView<Eigen::Lower>().solve(DPC[i])).transpose();

    } // end for i=0:n_horizon

}

// generate RHS of linear system
void QP_primal_dual_mpc::generateRHS()
{
    // compute r_C term
    r_C = P*y + q + C_T*nu + G_T*lambda;

    // compute r_E term
    r_E = C*y + crhs; // equality constraint: C*y + crhs = 0

    // compute r_I term
    r_I = (G*y - w) + s;

    // compute r_S term
    r_S = S*lambda; // no centering term

    // compute r_d term
    r_d = r_C + G_T*(S_inv_Lambda)*r_I - G_T*S_inv*r_S;

    // compute beta (RHS)
    solve_PHIx_eq_b(r_d, phixb);
    beta = r_E - C*phixb;
}

// update RHS
void QP_primal_dual_mpc::updateRHS()
{
    // update r_s term
    r_S = r_S - nu_centering;

    // update r_d
    r_d = r_C + G_T*(S_inv_Lambda)*r_I - G_T*S_inv*r_S;

    // upcate beta
    solve_PHIx_eq_b(r_d, phixb);
    beta = r_E - C*phixb;
}

// compute full step direction
void QP_primal_dual_mpc::computeFullStepDirectionFrom_del_nu()
{
    solve_PHIx_eq_b(-r_d - C_T*del_nu, del_y);
    del_s = -r_I - G*del_y;
    del_lambda = S_inv*(-r_S - Lambda*del_s);
}

// solve linear system PHI*x = b
void QP_primal_dual_mpc::solve_PHIx_eq_b(const Eigen::Ref<const Eigen::VectorXd> b, Eigen::Ref<Eigen::VectorXd> x)
{
    // declare
    int i_offset = 0;

    // solve by step
    for (int i=0; i <= n_horizon; i++)
    {

        if (isDiag_phi_i_vec[i]) 
        {
            // solve directly via inverse
            x.segment(i_offset, n_y_vec[i]) = phi_i_inv_vec[i]*(b.segment(i_offset, n_y_vec[i]));
        }
        else // solve via forward and backward substitution
        {
            // PHI_i = L_i*L_i'
            // solve PHI_i*x_i = b_i
            solveCholeskyLinSys(L_phi_i_vec[i], b.segment(i_offset, n_y_vec[i]), x.segment(i_offset, n_y_vec[i]));
        }

        // update offset
        i_offset += n_y_vec[i];
    }
}

// solve linear system Y*del_nu = beta
void QP_primal_dual_mpc::solveYdelnu_eq_beta()
{
    // solve by time step
    int i_st = 0;
    int im1_st = 0;
    Eigen::VectorXd del = Eigen::VectorXd::Zero(n_nu_tot); // intermediate solution  

    // forward substitution
    for (int k=1; k<=n_horizon; k++)
    {
        if (k == 1)
        {
            del.segment(i_st, n_nu_vec[k-1]) = L_i_i_vec[k].triangularView<Eigen::Lower>().solve(
                beta.segment(i_st, n_nu_vec[k-1]));
        }
        else
        {
            del.segment(i_st, n_nu_vec[k-1]) = L_i_i_vec[k].triangularView<Eigen::Lower>().solve(
                beta.segment(i_st, n_nu_vec[k-1])-L_ip1_i_vec[k-1]*del.segment(im1_st, n_nu_vec[k-2]));
            im1_st += n_nu_vec[k-2];
        }
        i_st += n_nu_vec[k-1];
    }

    // backward substitution
    i_st = n_nu_tot - n_nu_vec[n_horizon-1];
    int ip1_st = i_st;
    for (int k=n_horizon; k>=1; k--)
    {
        if (k == n_horizon)
        {
            del_nu.segment(i_st, n_nu_vec[k-1]) = L_i_i_vec[k].transpose().triangularView<Eigen::Upper>().solve(
                del.segment(i_st, n_nu_vec[k-1]));
        }
        else
        {
            del_nu.segment(i_st, n_nu_vec[k-1]) = L_i_i_vec[k].transpose().triangularView<Eigen::Upper>().solve(
                del.segment(i_st, n_nu_vec[k-1])-L_ip1_i_vec[k].transpose()*del_nu.segment(ip1_st, n_nu_vec[k]));
            ip1_st -= n_nu_vec[k-1];
        }
        i_st -= n_nu_vec[k-2];
    }
}

// solve linear system
void QP_primal_dual_mpc::solveCholeskyLinSys(const Eigen::SparseMatrix<double> * L_lowerTriangular_ptr, const Eigen::Ref<const Eigen::VectorXd> b, Eigen::Ref<Eigen::VectorXd> x)
{
    // solves (L*L')*x = b

    // forward and backward substitution
    x = L_lowerTriangular_ptr->transpose().triangularView<Eigen::Upper>().solve((L_lowerTriangular_ptr->triangularView<Eigen::Lower>().solve(b)));
}

void QP_primal_dual_mpc::solveCholeskyLinSys(const Eigen::Ref<const Eigen::MatrixXd> L_lowerTriangular, const Eigen::Ref<const Eigen::VectorXd> b, Eigen::Ref<Eigen::VectorXd> x)
{
    // solves (L*L')*x = b

    // forward and backward substitution
    x = L_lowerTriangular.transpose().triangularView<Eigen::Upper>().solve((L_lowerTriangular.triangularView<Eigen::Lower>().solve(b)));
}


// other methods
std::pair<double, bool> QP_primal_dual_mpc::lineSearch(const Eigen::Ref<const Eigen::VectorXd> del_s_in, const Eigen::Ref<const Eigen::VectorXd> del_lambda_in)
{
    // declare
    Eigen::VectorXd s_ds, lamda_dlambda;

    // init
    double h = 1;
    bool valid = false;
    int cnt_max = 10000;
    int cnt = 0;

    // loop
    while (!valid && (cnt < cnt_max) && (h > eps))
    {
        // updated s, u
        s_ds = s + h*del_s_in;
        lamda_dlambda = lambda + h*del_lambda_in;

        // check for validity of step
        if ((s_ds.minCoeff() > 0) && (lamda_dlambda.minCoeff() > 0))
            valid = true;
        else
            h = settings.t_ls*h;

        // increment
        cnt++;
    }

    // check for validity
    bool numerical_issue = ((cnt == cnt_max) || (h <= eps));

    // output
    return std::pair<double, bool> {h, numerical_issue};
}

// compute problem dimensions
void QP_primal_dual_mpc::computeProblemDimensions()
{
    // get problem dimensions at each time step
    for (int i=0; i <= n_horizon; i++)
    {
        if (P_i_vec.count(i))
            n_y_vec[i] = P_i_vec[i].rows();
        else
            n_y_vec[i] = P_i_nom.rows();

        if (G_i_vec.count(i))
            n_lambda_s_vec[i] = G_i_vec[i].rows();
        else
            n_lambda_s_vec[i] = G_i_nom.rows();

        if (i != n_horizon)
        {
            if (C_i_vec.count(i))
                n_nu_vec[i] = C_i_vec[i].rows();
            else
                n_nu_vec[i] = C_i_nom.rows();
        }
    }

    // add up dimension to compute total
    n_y_tot = 0; // init
    for (auto it = n_y_vec.begin(); it != n_y_vec.end(); it++)
        n_y_tot += it->second;

    n_lambda_s_tot = 0; // init
    for (auto it = n_lambda_s_vec.begin(); it != n_lambda_s_vec.end(); it++)
        n_lambda_s_tot += it->second;
    
    n_nu_tot = 0; // init
    for (auto it = n_nu_vec.begin(); it != n_nu_vec.end(); it++)
        n_nu_tot += it->second;
}

// initialize working matrices
void QP_primal_dual_mpc::initializeWorkingMatrices()
{   
    // initialize constant size matrices and triplet vectors
    del_nu = Eigen::VectorXd::Zero(n_nu_tot);
    del_y = Eigen::VectorXd::Zero(n_y_tot);
    phixb = Eigen::VectorXd::Zero(n_y_tot);
}

// build matrices
void QP_primal_dual_mpc::build_P()
{
    // make P
    std::vector<Eigen::Triplet<double>> tripvec_P;
    int nnz_P = 0;
    for (int i=0; i <= n_horizon; i++)
    {
        if (P_i_vec.count(i))
            nnz_P += P_i_vec[i].nonZeros();
        else
            nnz_P += P_i_nom.nonZeros();
    }
    tripvec_P.reserve(nnz_P);

    Eigen::SparseMatrix<double> * P_i_ptr=nullptr;
    int i_offset_P = 0;
    for (int i=0; i <= n_horizon; i++)
    {
        if (P_i_vec.count(i))
            P_i_ptr = &(P_i_vec[i]);
        else
            P_i_ptr = &(P_i_nom);
        getTripletsForMatrix(P_i_ptr, tripvec_P, i_offset_P, i_offset_P);
        i_offset_P += P_i_ptr->rows();
    }
    P.resize(i_offset_P, i_offset_P);
    P.setFromTriplets(tripvec_P.begin(), tripvec_P.end());

    // update flag
    P_updated = false;
}

void QP_primal_dual_mpc::build_q()
{
    // make q
    q.resize(n_y_tot);
    Eigen::VectorXd * q_i_ptr=nullptr;
    int i_offset_q = 0;
    for (int i=0; i <= n_horizon; i++)
    {
        if (q_i_vec.count(i))
            q_i_ptr = &(q_i_vec[i]);
        else
            q_i_ptr = &(q_i_nom);
        q.segment(i_offset_q, q_i_ptr->rows()) = *q_i_ptr;
        i_offset_q += q_i_ptr->rows();
    }

    // update flag
    q_updated = false;
}

void QP_primal_dual_mpc::build_C()
{
    // make C
    std::vector<Eigen::Triplet<double>> tripvec_C;
    int nnz_C = 0;
    for (int i=0; i < n_horizon; i++) // not inclusive of n_horizon
    {
        if (C_i_vec.count(i))
            nnz_C += C_i_vec[i].nonZeros();
        else
            nnz_C += C_i_nom.nonZeros();
        
        if (D_i_vec.count(i+1))
            nnz_C += D_i_vec[i+1].nonZeros();
        else
            nnz_C += D_i_nom.nonZeros();
    }
    tripvec_C.reserve(nnz_C);

    Eigen::SparseMatrix<double> * C_i_ptr = nullptr;
    Eigen::SparseMatrix<double> * D_i_ptr = nullptr;
    int i_offset_C = 0;
    int j_offset_C = 0;
    int j_offset_D = 0;
    for (int i=0; i < n_horizon; i++) // not inclusive of n_horizon
    {
        if (C_i_vec.count(i))
            C_i_ptr = &(C_i_vec[i]);
        else
            C_i_ptr = &(C_i_nom);

        if (D_i_vec.count(i+1))
            D_i_ptr = &(D_i_vec[i+1]);
        else
            D_i_ptr = &(D_i_nom);

        getTripletsForMatrix(C_i_ptr, tripvec_C, i_offset_C, j_offset_C);

        j_offset_D += C_i_ptr->cols();
        getTripletsForMatrix(D_i_ptr, tripvec_C, i_offset_C, j_offset_D);

        i_offset_C += C_i_ptr->rows();
        j_offset_C = j_offset_D;
    }
    C.resize(i_offset_C, j_offset_C + D_i_ptr->cols());

    C.setFromTriplets(tripvec_C.begin(), tripvec_C.end());

    // precompute transposes
    C_i_nom_T = C_i_nom.transpose();
    D_i_nom_T = D_i_nom.transpose();
    for (auto it=C_i_vec.begin(); it!=C_i_vec.end(); it++)
        C_i_vec_T[it->first] = it->second.transpose();
    for (auto it=D_i_vec.begin(); it!=D_i_vec.end(); it++)
        D_i_vec_T[it->first] = it->second.transpose();
    C_T = C.transpose();

    // precompute sparse->dense matrix conversions
    C_i_nom_T_dense = Eigen::MatrixXd(C_i_nom_T);
    D_i_nom_T_dense = Eigen::MatrixXd(D_i_nom_T);
    for (auto it=C_i_vec_T.begin(); it!=C_i_vec_T.end(); it++)
        C_i_vec_T_dense[it->first] = Eigen::MatrixXd(it->second);
    for (auto it=D_i_vec_T.begin(); it!=D_i_vec_T.end(); it++)
        D_i_vec_T_dense[it->first] = Eigen::MatrixXd(it->second);

    // update flag
    C_updated = false;
}

void QP_primal_dual_mpc::build_crhs()
{
    // make crhs
    crhs.resize(n_nu_tot);
    Eigen::VectorXd * crhs_i_ptr=nullptr;
    int i_offset_crhs = 0;
    for (int i=0; i < n_horizon; i++) // not inclusive of n_horizon
    {
        if (crhs_i_vec.count(i))
            crhs_i_ptr = &(crhs_i_vec[i]);
        else
            crhs_i_ptr = &(crhs_i_nom);
        crhs.segment(i_offset_crhs, crhs_i_ptr->rows()) = *crhs_i_ptr;
        i_offset_crhs += crhs_i_ptr->rows();
    }

    // update flag
    crhs_updated = false;
}

void QP_primal_dual_mpc::build_G()
{
    // make G
    std::vector<Eigen::Triplet<double>> tripvec_G;
    int nnz_G = 0;
    for (int i=0; i <= n_horizon; i++)
    {
        if (G_i_vec.count(i))
            nnz_G += G_i_vec[i].nonZeros();
        else
            nnz_G += G_i_nom.nonZeros();
    }
    tripvec_G.reserve(nnz_G);

    Eigen::SparseMatrix<double> * G_i_ptr=nullptr;
    int i_offset_G = 0;
    int j_offset_G = 0;
    for (int i=0; i <= n_horizon; i++)
    {
        if (G_i_vec.count(i))
            G_i_ptr = &(G_i_vec[i]);
        else
            G_i_ptr = &(G_i_nom);
        getTripletsForMatrix(G_i_ptr, tripvec_G, i_offset_G, j_offset_G);
        i_offset_G += G_i_ptr->rows();
        j_offset_G += G_i_ptr->cols();
    }
    G.resize(i_offset_G, j_offset_G);
    G.setFromTriplets(tripvec_G.begin(), tripvec_G.end());

    // precompute transposes
    G_i_nom_T = G_i_nom.transpose();
    for (auto it=G_i_vec.begin(); it!=G_i_vec.end(); it++)
        G_i_vec_T[it->first] = it->second.transpose();
    G_T = G.transpose();

    // update flag
    G_updated = false;
}

void QP_primal_dual_mpc::build_w()
{
    // make w
    w.resize(n_lambda_s_tot);
    Eigen::VectorXd * w_i_ptr;
    int i_offset_w = 0;
    for (int i=0; i <= n_horizon; i++)
    {
        if (w_i_vec.count(i))
            w_i_ptr = &(w_i_vec[i]);
        else
            w_i_ptr = &(w_i_nom);
        w.segment(i_offset_w, w_i_ptr->rows()) = *w_i_ptr;
        i_offset_w += w_i_ptr->rows();
    }

    // update flag
    w_updated = false;
}

// duality measure 
double QP_primal_dual_mpc::computeMu(const Eigen::Ref<const Eigen::VectorXd> s_in, const Eigen::Ref<const Eigen::VectorXd> lambda_in)
{
    return (s_in.dot(lambda_in))/n_lambda_s_tot;
}

// objective function
double QP_primal_dual_mpc::objective(const Eigen::Ref<const Eigen::VectorXd> x_in)
{
    double J = 0.5*x_in.dot(P*x_in) + q.dot(x_in) + b;
    return J;
}

// make S_inv_lambda
bool QP_primal_dual_mpc::make_S_inv_Lambda()
{
    // compute quantities of interest as arrays
    Eigen::ArrayXd s_arr(s);
    Eigen::ArrayXd lambda_arr(lambda);
    Eigen::ArrayXd s_inv_arr = 1.0/s_arr;
    Eigen::ArrayXd s_inv_lambda_arr = s_inv_arr*lambda_arr;

    // convert to diagonal matrices
    S_inv_Lambda.diagonal() = s_inv_lambda_arr;
    S_inv.diagonal() = s_inv_arr;
    Lambda.diagonal() = lambda_arr;
    S.diagonal() = s_arr;

    // get s_inv_lambda by time step
    int i_st = 0;
    for (int k=0; k<=n_horizon; k++)
    {
        S_inv_Lambda_vec[k].diagonal() = s_inv_lambda_arr.segment(i_st, n_lambda_s_vec[k]);
        i_st += n_lambda_s_vec[k];
    }

    // check for numerical limit
    bool numerical_limit = (s_inv_lambda_arr.minCoeff() <= eps);
    return numerical_limit;
}


// utilities
void inline QP_primal_dual_mpc::getTripletsForMatrix(const Eigen::Ref<const Eigen::MatrixXd> mat, std::vector<Eigen::Triplet<double>> &tripvec,
      int rowOffset, int colOffset)
{
    int m = mat.rows();
    int n = mat.cols();
    for (int i=0; i<m; i++)
    {
        for (int j=0; j<n; j++)
        {
            if (mat(i,j) != 0)
                tripvec.push_back(Eigen::Triplet<double>(i+rowOffset, j+colOffset, mat(i,j)));
        }
    }
}

void inline QP_primal_dual_mpc::getTripletsForMatrix(const Eigen::SparseMatrix<double> * mat_ptr, std::vector<Eigen::Triplet<double>> &tripvec,
      int rowOffset, int colOffset)
{
    int m = mat_ptr->rows();
    int n = mat_ptr->cols();
    
    for (int i=0; i<mat_ptr->outerSize(); i++)
    {
        for (typename Eigen::SparseMatrix<double>::InnerIterator it(*mat_ptr,i); it; ++it)
            tripvec.push_back(Eigen::Triplet<double>(it.row()+rowOffset, it.col()+colOffset, it.value()));
    }
}

bool QP_primal_dual_mpc::checkIfBox(const Eigen::SparseMatrix<double> * M_ptr)
{
    // optional bypass for sim study
    #if NODIAG
    return false;
    #endif

    bool isBox = true; // init
    bool * boxVec = new bool [M_ptr->rows()] {}; // initialize array of zeros

    for (int k=0; k<M_ptr->outerSize(); ++k)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator it(*M_ptr,k); it; ++it)
        {
            if (!boxVec[it.row()])
                boxVec[it.row()] = true;
            else
                isBox = false;
        }
    }

    delete [] boxVec;

    return isBox;
}

bool QP_primal_dual_mpc::checkIfDiag(const Eigen::SparseMatrix<double> * M_ptr)
{
    // optional bypass for sim study
    #if NODIAG
    return false;
    #endif

    for (int k=0; k<M_ptr->outerSize(); ++k)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator it(*M_ptr,k); it; ++it)
        {
            if (it.row() != it.col())
                return false;
        }
    }
    return true;
}