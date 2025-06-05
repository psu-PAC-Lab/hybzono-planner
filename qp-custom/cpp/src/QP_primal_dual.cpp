#include "QP_primal_dual.hpp"

using namespace QP;

// constructor
QP_primal_dual::QP_primal_dual() {}

QP_primal_dual::QP_primal_dual(const Eigen::SparseMatrix<double> &P_in,
            const Eigen::VectorXd &q_in,
            const Eigen::SparseMatrix<double> &A_in,
            const Eigen::VectorXd &b_in,
            const Eigen::SparseMatrix<double> &G_in,
            const Eigen::VectorXd &w_in)
{
    // setup
    setup(P_in, q_in, A_in, b_in, G_in, w_in);
}

// destructor
QP_primal_dual::~QP_primal_dual() = default;

// setup
void QP_primal_dual::setup(const Eigen::SparseMatrix<double> &P_in,
            const Eigen::VectorXd &q_in,
            const Eigen::SparseMatrix<double> &A_in,
            const Eigen::VectorXd &b_in,
            const Eigen::SparseMatrix<double> &G_in,
            const Eigen::VectorXd &w_in)
{
    // copy in vars
    P = P_in;
    q = q_in;
    A = A_in;
    b = b_in;
    G = G_in;
    w = w_in;

    // preprocessing (make sure A is full rank)
    if (settings.preprocessing_enable)
    {
        getValidEqualityConstraints();
        makeValid_A();
        makeValid_b();

        A_updated = false;
        b_updated = false;
    }

    // compute problem dimensions
    computeProblemDimensions();
}

// set settings
void QP_primal_dual::set_settings(const QP_settings *settings_in)
{
    // copy over settings
    settings = *settings_in;
}

// set methods
void QP_primal_dual::set_P(const Eigen::SparseMatrix<double> &P_in)
{
    P = P_in;
    if (P.rows() != n) // change in problem dimension
    {
        x0 = Eigen::VectorXd::Zero(P.rows()); // reset IC
        n = P.rows();
    }
}

void QP_primal_dual::set_q(const Eigen::VectorXd &q_in)
{
    q = q_in;
}

void QP_primal_dual::set_A(const Eigen::SparseMatrix<double> &A_in)
{    
    A = A_in;
    A_updated = true;
}

void QP_primal_dual::set_b(const Eigen::VectorXd &b_in)
{
    b = b_in;
    b_updated = true;
}

void QP_primal_dual::set_G(const Eigen::SparseMatrix<double> &G_in)
{
    G = G_in;
}

void QP_primal_dual::set_w(const Eigen::VectorXd &w_in)
{
    w = w_in;
}

// solve
QP_results QP_primal_dual::solve()
{
    // declare
    double mu, h, mu_pred, sigma;
    Eigen::VectorXd del, del_x, del_v, del_u, del_s;
    Eigen::VectorXd s_pred, u_pred;
    bool numerical_issue;
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

    // pre-processing if matrices were updated
    if (settings.preprocessing_enable)
    {
        if (A_updated)
        {
            getValidEqualityConstraints();
            makeValid_A();
            makeValid_b();
            computeProblemDimensions();
        }
        else if (b_updated)
        {
            makeValid_b();
            computeProblemDimensions();
        }

        A_updated = false; // reset flag
        b_updated = false; // reset flag
    }
        
    // make sure problem dimensions are correct in case matrices changed
    computeProblemDimensions();

    // initialize start vars
    x0 = Eigen::VectorXd::Zero(n);
    double zeta = sqrt(settings.mu_init);
    v0 = Eigen::VectorXd::Zero(m_eq);
    u0 = zeta*Eigen::VectorXd::Ones(m_ineq);
    s0 = zeta*Eigen::VectorXd::Ones(m_ineq);

    // initialize primal and dual vars
    x = x0;
    v = v0;
    u = u0;
    s = s0;

    // initialize working matrices
    initializeWorkingMatrices();

    // outer loop init
    int k = 0;
    numerical_issue = false;

    // compute duality measure
    mu = computeMu(s, u);

    // loop
    while ((mu > settings.mu_term) && (k < settings.iter_max) && !numerical_issue && (running_timer < T_max) && (mu <= settings.mu_max))
    {
        // generate system matrix and decompose
        generateSystemMatrix();

        // check for numerical issues and terminate if necessary
        if (LU_status != Eigen::ComputationInfo::Success)
        {
            numerical_issue = true;
            break;
        } 

        // predictor step
        generateRHS();
        del = LU_solver.solve(bm);
        del_u = del.segment(n+m_eq, m_ineq);
        del_s = del.segment(n+m_eq+m_ineq, m_ineq);

        // line search for step size
        line_search_out = lineSearch(del_s, del_u);
        h = line_search_out.first;
        if (line_search_out.second)
        {
            numerical_issue = true;
            break;
        }

        // predicted duality measure
        s_pred = s + h*del_s;
        u_pred = u + h*del_u;
        mu_pred = computeMu(s_pred, u_pred);

        // centering parameter
        sigma = pow(mu_pred/mu, 3);

        // calculate corrected nu and recompute search direction
        Del_S.diagonal() = del_s;
        nu = (sigma*mu)*Eigen::VectorXd::Ones(m_ineq) - Del_S*del_u;
        updateRHS();

        del = LU_solver.solve(bm);
        del_x = del.segment(0, n);
        del_v = del.segment(n, m_eq);
        del_u = del.segment(n+m_eq, m_ineq);
        del_s = del.segment(n+m_eq+m_ineq, m_ineq);

        // line search for step size
        line_search_out = lineSearch(del_s, del_u);
        h = line_search_out.first;
        if (line_search_out.second)
        {
            numerical_issue = true;
            break;
        }

        // update
        x = x + settings.gamma*h*del_x;
        v = v + settings.gamma*h*del_v;
        u = u + settings.gamma*h*del_u;
        s = s + settings.gamma*h*del_s;

        // update running timer
        if (enforce_max_time)
        {
            auto timer_running = std::chrono::high_resolution_clock::now();
            auto duration_running = std::chrono::duration_cast<std::chrono::microseconds>(timer_running - timer_init);
            running_timer = 1e-6 * ((double) duration_running.count());
        }

        // compute duality measure
        mu = computeMu(s, u);

        // iterate 
        k++;
    }

    // timing
    auto timer_final = std::chrono::high_resolution_clock::now();
    auto duration_final = std::chrono::duration_cast<std::chrono::microseconds>(timer_final - timer_init);
    double time = 1e-6 * ((double) duration_final.count());

    // assemble results
    QP_results results;
    results.x = x;
    results.v = v;
    results.u = u;
    results.s = s;
    results.objective = objective(x);
    results.feas = mu <= settings.mu_feas; // divergence check
    results.converged = !numerical_issue && (k < settings.iter_max) && results.feas;
    results.num_iter = k;
    results.sol_time = time;

    // return
    return results;
}

/* helper methods */

// objective function
double QP_primal_dual::objective(const Eigen::Ref<const Eigen::VectorXd> x_in)
{
    double J;
    J = 0.5*x_in.dot(P*x_in) + q.dot(x_in);
    return J;
}

// initialize working matrices
void QP_primal_dual::initializeWorkingMatrices()
{
    // precompute transposes
    A_T = A.transpose();
    G_T = G.transpose();

    // initialize working matrices
    bm.resize(n+m_eq+m_ineq+m_ineq);

    // precompute constant part of M matrix
    // M = [P, A', G', 0;
    //      A, 0, 0, 0;
    //      G, 0, 0, I;
    //      0, 0, S, Z]
    std::vector<Eigen::Triplet<double>> tripvec_M0;
    int nnz_M = P.nonZeros() + 2*A.nonZeros() + 2*G.nonZeros() + m_ineq;
    tripvec_M0.reserve(nnz_M);
    M0.resize(n+m_eq+2*m_ineq, n+m_eq+2*m_ineq);

    getTripletsForMatrix(&P, tripvec_M0, 0, 0);
    if (equalityConstrained)
        getTripletsForMatrix(&A_T, tripvec_M0, 0, n);
    getTripletsForMatrix(&G_T, tripvec_M0, 0, n+m_eq);

    if (equalityConstrained)
        getTripletsForMatrix(&A, tripvec_M0, n, 0);

    getTripletsForMatrix(&G, tripvec_M0, n+m_eq, 0);
    getTripletsForMatrixDiagonal(Eigen::VectorXd::Ones(m_ineq), tripvec_M0, n+m_eq, n+m_eq+m_ineq);

    M0.setFromTriplets(tripvec_M0.begin(), tripvec_M0.end());

    // pre-allocate and initialize dM
    tripvec_dM.clear();
    tripvec_dM.reserve(2*m_ineq);
    dM.resize(n+m_eq+2*m_ineq, n+m_eq+2*m_ineq);
    getTripletsForMatrixDiagonal(s, tripvec_dM, n+m_eq+m_ineq, n+m_eq);
    getTripletsForMatrixDiagonal(u, tripvec_dM, n+m_eq+m_ineq, n+m_eq+m_ineq);
    dM.setFromTriplets(tripvec_dM.begin(), tripvec_dM.end());

    // analyze pattern
    M = M0 + dM;
    LU_solver.analyzePattern(M);
}

// generate system matrix
void QP_primal_dual::generateSystemMatrix()
{
    // update changing part of system matrix
    // update dM based on known sparsity pattern
    int i_s = 0;
    int i_u = 0;
    for (int k=0; k<dM.outerSize(); k++)
    {
        for (typename Eigen::SparseMatrix<double>::InnerIterator it(dM,k); it; ++it)
        {
            if (i_s < s.rows())
            {
                it.valueRef() = s(i_s);
                i_s++;
            }
            else if (i_u < u.rows())
            {
                it.valueRef() = u(i_u);
                i_u++;
            }
        }
    }
    
    // update M
    M = M0 + dM;

    // LU decomposition
    LU_solver.factorize(M);

    // get status
    LU_status = LU_solver.info();

}

// generate right hand side of linear system
void QP_primal_dual::generateRHS()
{
    // compute r_E term
    if (equalityConstrained)
    {
        // compute r_C term
        r_C = P*x + q + A_T*v + G_T*u;

        // compute r_E term
        r_E = A*x - b;
    }
    else
        r_C = P*x + q + G_T*u;

    // compute r_I term
    r_I = (G*x - w) + s;

    // compute r_S term
    S.diagonal() = s;
    r_S = S*u; // no centering term

    // RHS
    bm.segment(0, n) = -r_C;
    if (equalityConstrained)
        bm.segment(n, m_eq) = -r_E;
    bm.segment(n+m_eq, m_ineq) = -r_I;
    bm.segment(n+m_eq+m_ineq, m_ineq) = -r_S;
}

// update RHS
void QP_primal_dual::updateRHS()
{
    // update r_S term
    r_S = r_S - nu;
    
    // update RHS
    bm.segment(n+m_eq+m_ineq, m_ineq) = -r_S;
}

// line search
std::pair<double, bool> QP_primal_dual::lineSearch(const Eigen::Ref<const Eigen::VectorXd> del_s, const Eigen::Ref<const Eigen::VectorXd> del_u)
{
    // declare
    Eigen::VectorXd s_ds, u_du;

    // init
    double h = 1;
    bool valid = false;
    int cnt_max = 10000;
    int cnt = 0;

    // loop
    while (!valid && (cnt < cnt_max) && (h > eps))
    {
        // updated s, u
        s_ds = s + h*del_s;
        u_du = u + h*del_u;

        // check for validity of step
        if ((s_ds.minCoeff() > 0) && (u_du.minCoeff() > 0))
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

// duality measure 
double QP_primal_dual::computeMu(const Eigen::Ref<const Eigen::VectorXd> s_in, const Eigen::Ref<const Eigen::VectorXd> u_in)
{
    return (s_in.dot(u_in))/m_ineq;
}

// pre-processing
void QP_primal_dual::getLinDepPermuteAndChopMatrix(const Eigen::SparseMatrix<double> * mat_in, Eigen::SparseMatrix<double> * mat_out)
{
    // check for empty input matrix
    if (mat_in->rows() == 0 || mat_in->cols() == 0)
    {
        mat_out->resize(0, 0);
        return;
    }

    // matrix transpose
    Eigen::SparseMatrix<double> mat_T = mat_in->transpose();

    // compute QR decomposition
    QR_solver.analyzePattern(mat_T);
    QR_solver.factorize(mat_T);

    // get permutation matrix and its indices
    Eigen::PermutationMatrix<-1, -1> P_full = QR_solver.colsPermutation();
    Eigen::VectorXi ind_full = P_full.indices();

    // construct permute and chop matrix directly
    std::vector<Eigen::Triplet<double>> tripvec;
    tripvec.reserve(ind_full.size());
    for (int i=0; i<QR_solver.rank(); i++) // QR solver automatically puts linearly dependent rows at end
        tripvec.push_back(Eigen::Triplet<double>(ind_full[i], i, 1));

    mat_out->resize(mat_T.cols(), QR_solver.rank());
    mat_out->setFromTriplets(tripvec.begin(), tripvec.end());
}

void QP_primal_dual::computeProblemDimensions()
{
    n = P.rows();
    m_eq = A.rows();
    m_ineq = G.rows();
    equalityConstrained = (m_eq == 0) ? false : true;
}

void QP_primal_dual::getValidEqualityConstraints()
{
    // remove any invalid constraints (linearly dependent)
    getLinDepPermuteAndChopMatrix(&A, &P_eq);
}

void QP_primal_dual::makeValid_A()
{
    A = (A.transpose()*P_eq).transpose();
}

void QP_primal_dual::makeValid_b()
{
    b = (b.transpose()*P_eq).transpose();
}

// utilities
void inline QP_primal_dual::getTripletsForMatrix(const Eigen::Ref<const Eigen::MatrixXd> mat, std::vector<Eigen::Triplet<double>> &tripvec,
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

void inline QP_primal_dual::getTripletsForMatrix(const Eigen::SparseMatrix<double> * mat_ptr, std::vector<Eigen::Triplet<double>> &tripvec,
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

void inline QP_primal_dual::getTripletsForMatrixDiagonal(const Eigen::Ref<const Eigen::VectorXd> d, std::vector<Eigen::Triplet<double>> &tripvec,
            int rowOffset, int colOffset)
{
    int m = d.rows();
    for (int i=0; i<m; i++)
    {
        tripvec.push_back(Eigen::Triplet<double>(i+rowOffset, i+colOffset, d(i)));
    }
}