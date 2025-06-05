#include "mex_set_methods.hpp"

// SET METHODS
void MEXSET::getFullMatrixInput(const mxArray* M_ptr, Eigen::MatrixXd& M_out)
{
    // get size of matrix
    const mwSize * M_dim_ptr = mxGetDimensions(M_ptr);
    int m_M = (int) M_dim_ptr[0];
    int n_M = (int) M_dim_ptr[1];

    // check for empty G matrix
    if (!mxIsEmpty(M_ptr))
    {
        // get matrix data
        double * M_dbl = (double *) mxCalloc(m_M*n_M, sizeof(double));
        copyToDoubleVector(mxGetPr(M_ptr), m_M*n_M, M_dbl);
        if (M_out.rows() != m_M || M_out.cols() != n_M)
            M_out.resize(m_M, n_M);
        buildEigenMatrixFromDblPtr(M_dbl, m_M, n_M, M_out);
        mxFree(M_dbl);
    }
    else
    {
        M_out.resize(m_M, n_M); // one of these should be 0
    }
}

void MEXSET::getFullMatrixMap(const mxArray* M_i_cell_ptr, const mxArray* ind_M_i_ptr, std::map<int, Eigen::MatrixXd> &M_i_vec)
{
    // get number of entries
    const mwSize * ind_M_dim_ptr = mxGetDimensions(ind_M_i_ptr);
    int m_M_vec = (int) ind_M_dim_ptr[0];
    int n_M_vec = (int) ind_M_dim_ptr[1];
    if (m_M_vec > n_M_vec)
        n_M_vec = m_M_vec; // in case vector comes in as a col vector
    
    // copy indices to array
    double * ind_M_i_array = (double *) mxCalloc(n_M_vec, sizeof(double));
    copyToDoubleVector(mxGetPr(ind_M_i_ptr), n_M_vec, ind_M_i_array);

    // loop through and add elements
    int ind_M;
    mxArray * M_i_cell;
    double * M_i_dbl_cell;
    const mwSize * M_cell_dim_ptr;
    int m_M_cell, n_M_cell;
    Eigen::MatrixXd M_cell;

    for (int i=0; i<n_M_vec; i++)
    {
        // get MPC step index for off-nominal matrix
        ind_M = (int) ind_M_i_array[i];

        // get pointer to cell array data
        M_i_cell = mxGetCell(M_i_cell_ptr, i);

        // get dimensions of matrix
        M_cell_dim_ptr = mxGetDimensions(M_i_cell);
        m_M_cell = (int) M_cell_dim_ptr[0];
        n_M_cell = (int) M_cell_dim_ptr[1];

        // get data
        M_i_dbl_cell = (double *) mxCalloc(m_M_cell*n_M_cell, sizeof(double)); // init
        copyToDoubleVector(mxGetPr(M_i_cell), m_M_cell*n_M_cell, M_i_dbl_cell);

        // copy to matrix
        if (M_cell.rows() != m_M_cell || M_cell.cols() != n_M_cell)
            M_cell.resize(m_M_cell, n_M_cell);
        buildEigenMatrixFromDblPtr(M_i_dbl_cell, m_M_cell, n_M_cell, M_cell);
            
        // free memory
        mxFree(M_i_dbl_cell);

        // add to map
        M_i_vec[ind_M] = M_cell;        
    }

    // free memory
    mxFree(ind_M_i_array);
}



void MEXSET::getSparseMatrixMap(const mxArray* M_i_cell_ptr, const mxArray* ind_M_i_ptr, std::map<int, Eigen::SparseMatrix<double>> &M_i_vec)
{
    // get number of entries
    const mwSize * ind_M_dim_ptr = mxGetDimensions(ind_M_i_ptr);
    int m_M_vec = (int) ind_M_dim_ptr[0];
    int n_M_vec = (int) ind_M_dim_ptr[1];
    if (m_M_vec > n_M_vec)
        n_M_vec = m_M_vec; // in case vector comes in as a col vector

    // copy indices to array
    double * ind_M_i_array = (double *) mxCalloc(n_M_vec, sizeof(double));
    copyToDoubleVector(mxGetPr(ind_M_i_ptr), n_M_vec, ind_M_i_array);

    // loop through and add elements
    int ind_M, nnz_M_cell;
    mxArray * M_i_cell;
    int * Mp_cell;
    int * Mi_cell;
    double * Mx_cell;
    const mwSize * M_cell_dim_ptr;
    int m_M_cell, n_M_cell;
    Eigen::SparseMatrix<double> M_cell;
    std::vector<Eigen::Triplet<double>> tripvec_M_cell;

    for (int i=0; i<n_M_vec; i++)
    {
        // get MPC step index for off-nominal matrix
        ind_M = (int) ind_M_i_array[i];

        // get pointer to cell array data
        M_i_cell = mxGetCell(M_i_cell_ptr, i);

        // get dimensions of matrix
        M_cell_dim_ptr = mxGetDimensions(M_i_cell);
        m_M_cell = (int) M_cell_dim_ptr[0];
        n_M_cell = (int) M_cell_dim_ptr[1];

        // get matrix data
        Mp_cell = (int *) mxCalloc(n_M_cell+1, sizeof(int));
        copyToIntVector(mxGetJc(M_i_cell), n_M_cell+1, Mp_cell);
        nnz_M_cell = Mp_cell[n_M_cell]; // number of nonzeros

        Mi_cell = (int *) mxCalloc(nnz_M_cell, sizeof(int));
        copyToIntVector(mxGetIr(M_i_cell), nnz_M_cell, Mi_cell);

        Mx_cell = (double *) mxCalloc(nnz_M_cell, sizeof(double));
        copyToDoubleVector(mxGetPr(M_i_cell), nnz_M_cell, Mx_cell);

        // build P matrix
        tripvec_M_cell.clear();
        tripvec_M_cell.reserve(n_M_cell);
        buildSparseTripletsFromJcIrPr(Mp_cell, Mi_cell, Mx_cell, n_M_cell+1, tripvec_M_cell);
        M_cell.resize(m_M_cell, n_M_cell);
        M_cell.setFromTriplets(tripvec_M_cell.begin(), tripvec_M_cell.end());

        // free memory
        mxFree(Mp_cell);
        mxFree(Mi_cell);
        mxFree(Mx_cell);

        // add to map
        M_i_vec[ind_M] = M_cell;
            
    } // end for

    // free memory
    mxFree(ind_M_i_array);
}

void MEXSET::getVectorMap(const mxArray* x_i_cell_ptr, const mxArray* ind_x_i_ptr, std::map<int, Eigen::VectorXd> &x_i_vec)
{
    // get number of entries
    const mwSize * ind_x_dim_ptr = mxGetDimensions(ind_x_i_ptr);
    int m_x_vec = (int) ind_x_dim_ptr[0];
    int n_x_vec = (int) ind_x_dim_ptr[1];
    if (m_x_vec > n_x_vec)
        n_x_vec = m_x_vec; // in case vector comes in as a row vector

    // copy indices to array
    double * ind_x_i_array = (double *) mxCalloc(n_x_vec, sizeof(double));
    copyToDoubleVector(mxGetPr(ind_x_i_ptr), n_x_vec, ind_x_i_array);

    // loop through and add elements
    int ind_x;
    mxArray * x_i_cell;
    double * x_i_dbl_cell;
    const mwSize * x_cell_dim_ptr;
    int m_x_cell, n_x_cell;
    Eigen::VectorXd x_cell;

    for (int i=0; i<n_x_vec; i++)
    {
        // get MPC step index for off-nominal
        ind_x = (int) ind_x_i_array[i];

        // get pointer to cell array data
        x_i_cell = mxGetCell(x_i_cell_ptr, i);

        // get dimensions of matrix
        x_cell_dim_ptr = mxGetDimensions(x_i_cell);
        m_x_cell = (int) x_cell_dim_ptr[0];
        n_x_cell = (int) x_cell_dim_ptr[1];
        if (n_x_cell > m_x_cell)
            m_x_cell = n_x_cell; // protection against row vectors

        // get data
        x_i_dbl_cell = (double *) mxCalloc(m_x_cell, sizeof(double)); // init
        copyToDoubleVector(mxGetPr(x_i_cell), m_x_cell, x_i_dbl_cell);

        // copy to vector
        if (x_cell.rows() != m_x_cell)
            x_cell.resize(m_x_cell);
        buildEigenVectorFromDblPtr(x_i_dbl_cell, m_x_cell, x_cell);
            
        // free memory
        mxFree(x_i_dbl_cell);

        // add to map
        x_i_vec[ind_x] = x_cell;        
    }

    // free memory
    mxFree(ind_x_i_array);
}

void MEXSET::getVectorInput(const mxArray* x_ptr, Eigen::VectorXd &x_out)
{
    // get size of array
    const mwSize * x_dim_ptr = mxGetDimensions(x_ptr);
    int m_x = (int) x_dim_ptr[0];
    int n_x = (int) x_dim_ptr[1];
    if (m_x > n_x)
        n_x = m_x; // in case input comes in as a col vector

    // make sure array not empty
    if (!mxIsEmpty(x_ptr))
    {
        // make x vector
        double * x_ptr_dbl = (double *) mxCalloc(n_x, sizeof(double));
        copyToDoubleVector(mxGetPr(x_ptr), n_x, x_ptr_dbl);
        if (x_out.rows() != n_x)
            x_out.resize(n_x);
        buildEigenVectorFromDblPtr(x_ptr_dbl, n_x, x_out);
        mxFree(x_ptr_dbl);
    }
    else
    {
        x_out.resize(0);
    }
}

void MEXSET::getSparseMatrixInput(const mxArray* M_ptr, std::vector<Eigen::Triplet<double>> &tripvec, Eigen::SparseMatrix<double> * M_out)
{
    // clear triplet vector
    tripvec.clear();

    // get size of matrix
    const mwSize * M_dim_ptr = mxGetDimensions(M_ptr);
    int m_M = (int) M_dim_ptr[0];
    int n_M = (int) M_dim_ptr[1];

    // check for empty G matrix
    if (!mxIsEmpty(M_ptr))
    {
        // get matrix data
        int * Mp = (int *) mxCalloc(n_M+1, sizeof(int));
        copyToIntVector(mxGetJc(M_ptr), n_M+1, Mp);
        int nnz_M = Mp[n_M]; // number of nonzeros

        int * Mi = (int *) mxCalloc(nnz_M, sizeof(int));
        copyToIntVector(mxGetIr(M_ptr), nnz_M, Mi);

        double * Mx = (double *) mxCalloc(nnz_M, sizeof(double));
        copyToDoubleVector(mxGetPr(M_ptr), nnz_M, Mx);

        // get matrix triplets
        buildSparseTripletsFromJcIrPr(Mp, Mi, Mx, n_M+1, tripvec);

        // free memory
        mxFree(Mp);
        mxFree(Mi);
        mxFree(Mx);
    }
    else
    {
        m_M = 0;
        n_M = 0;
    }

    // resize if necessary
    if (M_out->rows() != m_M || M_out->cols() != n_M)
        M_out->resize(m_M, n_M);
    
    // fill matrix
    if (m_M != 0 && n_M != 0)
        M_out->setFromTriplets(tripvec.begin(), tripvec.end());
}

void MEXSET::getIntVectorInput(const mxArray* x_ptr, std::vector<int> &x_out)
{
    // get size of array
    const mwSize * x_dim_ptr = mxGetDimensions(x_ptr);
    int m_x = (int) x_dim_ptr[0];
    int n_x = (int) x_dim_ptr[1];
    if (m_x > n_x)
        n_x = m_x; // in case input comes in as a col vector

    // make sure array not empty
    if (!mxIsEmpty(x_ptr))
    {
        // make x vector
        double * x_ptr_dbl = (double *) mxCalloc(n_x, sizeof(double));
        copyToDoubleVector(mxGetPr(x_ptr), n_x, x_ptr_dbl);
        for (int i=0; i<n_x; i++)
            x_out.push_back((int) x_ptr_dbl[i]);
        mxFree(x_ptr_dbl);
    }
    else
    {
        x_out.clear();
    }
}

void MEXSET::getIntVectorInput_0_1_corrected(const mxArray* x_ptr, std::vector<int> &x_out)
{
    // get size of array
    const mwSize * x_dim_ptr = mxGetDimensions(x_ptr);
    int m_x = (int) x_dim_ptr[0];
    int n_x = (int) x_dim_ptr[1];
    if (m_x > n_x)
        n_x = m_x; // in case input comes in as a col vector

    // make sure array not empty
    if (!mxIsEmpty(x_ptr))
    {
        // make x vector
        double * x_ptr_dbl = (double *) mxCalloc(n_x, sizeof(double));
        copyToDoubleVector(mxGetPr(x_ptr), n_x, x_ptr_dbl);
        for (int i=0; i<n_x; i++)
            x_out.push_back((int) x_ptr_dbl[i]-1); // correct for 1-based indexing
        mxFree(x_ptr_dbl);
    }
    else
    {
        x_out.clear();
    }
}

void MEXSET::getIntVectorMap(const mxArray* x_i_cell_ptr, const mxArray* ind_x_i_ptr, std::map<int, std::vector<int>> &x_i_vec)
{
    // get number of entries
    const mwSize * ind_x_dim_ptr = mxGetDimensions(ind_x_i_ptr);
    int m_x_vec = (int) ind_x_dim_ptr[0];
    int n_x_vec = (int) ind_x_dim_ptr[1];
    if (m_x_vec > n_x_vec)
        n_x_vec = m_x_vec; // in case vector comes in as a row vector

    // copy indices to array
    double * ind_x_i_array = (double *) mxCalloc(n_x_vec, sizeof(double));
    copyToDoubleVector(mxGetPr(ind_x_i_ptr), n_x_vec, ind_x_i_array);

    // loop through and add elements
    int ind_x;
    mxArray * x_i_cell;
    double * x_i_dbl_cell;
    const mwSize * x_cell_dim_ptr;
    int m_x_cell, n_x_cell;
    std::vector<int> x_cell;

    for (int i=0; i<n_x_vec; i++)
    {
        // get MPC step index for off-nominal
        ind_x = (int) ind_x_i_array[i];

        // get pointer to cell array data
        x_i_cell = mxGetCell(x_i_cell_ptr, i);

        // get dimensions of matrix
        x_cell_dim_ptr = mxGetDimensions(x_i_cell);
        m_x_cell = (int) x_cell_dim_ptr[0];
        n_x_cell = (int) x_cell_dim_ptr[1];
        if (n_x_cell > m_x_cell)
            m_x_cell = n_x_cell; // protection against row vectors

        // protect against empty cell
        if (!mxIsEmpty(x_i_cell))
        {
            // get data
            x_i_dbl_cell = (double *) mxCalloc(m_x_cell, sizeof(double)); // init
            copyToDoubleVector(mxGetPr(x_i_cell), m_x_cell, x_i_dbl_cell);

            // copy to vector
            x_cell.clear();
            for (int j=0; j<m_x_cell; j++)
                x_cell.push_back((int) x_i_dbl_cell[j]);
            
            // free memory
            mxFree(x_i_dbl_cell);

            // add to map
            x_i_vec[ind_x] = x_cell;       

        } // end if not empty 
    }

    // free memory
    mxFree(ind_x_i_array);   
}

void MEXSET::getIntVectorMap_0_1_corrected(const mxArray* x_i_cell_ptr, const mxArray* ind_x_i_ptr, std::map<int, std::vector<int>> &x_i_vec)
{
    // get number of entries
    const mwSize * ind_x_dim_ptr = mxGetDimensions(ind_x_i_ptr);
    int m_x_vec = (int) ind_x_dim_ptr[0];
    int n_x_vec = (int) ind_x_dim_ptr[1];
    if (m_x_vec > n_x_vec)
        n_x_vec = m_x_vec; // in case vector comes in as a row vector

    // copy indices to array
    double * ind_x_i_array = (double *) mxCalloc(n_x_vec, sizeof(double));
    copyToDoubleVector(mxGetPr(ind_x_i_ptr), n_x_vec, ind_x_i_array);

    // loop through and add elements
    int ind_x;
    mxArray * x_i_cell;
    double * x_i_dbl_cell;
    const mwSize * x_cell_dim_ptr;
    int m_x_cell, n_x_cell;
    std::vector<int> x_cell;

    for (int i=0; i<n_x_vec; i++)
    {
        // get MPC step index for off-nominal
        ind_x = (int) ind_x_i_array[i];

        // get pointer to cell array data
        x_i_cell = mxGetCell(x_i_cell_ptr, i);

        // get dimensions of matrix
        x_cell_dim_ptr = mxGetDimensions(x_i_cell);
        m_x_cell = (int) x_cell_dim_ptr[0];
        n_x_cell = (int) x_cell_dim_ptr[1];
        if (n_x_cell > m_x_cell)
            m_x_cell = n_x_cell; // protection against row vectors

        // protect against empty cell
        if (!mxIsEmpty(x_i_cell))
        {
            // get data
            x_i_dbl_cell = (double *) mxCalloc(m_x_cell, sizeof(double)); // init
            copyToDoubleVector(mxGetPr(x_i_cell), m_x_cell, x_i_dbl_cell);

            // copy to vector
            x_cell.clear();
            for (int j=0; j<m_x_cell; j++)
                x_cell.push_back((int) x_i_dbl_cell[j]-1); // correct for 1-based indexing
            
            // free memory
            mxFree(x_i_dbl_cell);

            // add to map
            x_i_vec[ind_x] = x_cell;       

        } // end if not empty 
    }

    // free memory
    mxFree(ind_x_i_array);   
}

void MEXSET::R_x02region_processing(const mxArray* R_mat, std::map<int, std::vector<int>>& R_x02region)
{
    // get dimensions
    const mwSize* R_dim_ptr = mxGetDimensions(R_mat);
    int m_R = (int)R_dim_ptr[0];
    int n_R = (int)R_dim_ptr[1];

    // check for empty matrix
    if (!mxIsEmpty(R_mat))
    {
        // get matrix data
        double* R_dbl = (double*)mxCalloc(m_R * n_R, sizeof(double));
        copyToDoubleVector(mxGetPr(R_mat), m_R * n_R, R_dbl);

        // loop through and add R_dbl data to R_x02region
        for (int i = 0; i < m_R; i++)
        {
            for (int j = 0; j < n_R; j++)
            {
                if (R_dbl[i + j * m_R] != 0)
                {
                    // time step = i+1, region = j+1;
                    R_x02region[i + 1].push_back(j + 1);
                }
            }
        }
        mxFree(R_dbl);
    }
    else
    {
        R_x02region.clear();
    }
}

void MEXSET::R_region2region_processing(const mxArray* R_mat, std::map<int, std::map<int, std::vector<int>>>& R_region2region)
{
    // get dimensions
    const mwSize* R_dim_ptr = mxGetDimensions(R_mat);
    int m_R = (int)R_dim_ptr[0];
    int n_R = (int)R_dim_ptr[1];
    int k_R = (int)R_dim_ptr[2];

    // check for empty matrix
    if (!mxIsEmpty(R_mat))
    {
        // get matrix data
        double* R_dbl = (double*)mxCalloc(m_R * n_R * k_R, sizeof(double));
        copyToDoubleVector(mxGetPr(R_mat), m_R * n_R * k_R, R_dbl);

        // loop through and add R_dbl data to R_region2region
        for (int i = 0; i < m_R; i++)
        {
            for (int j = 0; j < n_R; j++)
            {
                for (int k = 0; k < k_R; k++)
                {
                    // add to R_region2region
                    if (R_dbl[i + j * m_R + k * m_R * n_R] != 0)
                    {
                        // time step = k+1, region = i+1, next region = j+1;
                        R_region2region[k + 1][i + 1].push_back(j + 1);
                    }
                }
            }
        }
        mxFree(R_dbl);
    }
    else
    {
        R_region2region.clear();
    }
}




// UTILITIES

// build sparse Eigen matrix triplets from jc, ir, pr
// ref https://www.mathworks.com/help/matlab/apiref/mxsetjc.html
// ref https://www.mathworks.com/help/matlab/apiref/mxgetir.html explore.c example

int MEXSET::getNNZ(const mxArray* M_ptr)
{
    // get size of matrix
    const mwSize * M_dim_ptr = mxGetDimensions(M_ptr);
    int m_M = (int) M_dim_ptr[0];
    int n_M = (int) M_dim_ptr[1];

    // get matrix data
    int * Mp = (int *) mxCalloc(n_M+1, sizeof(int));
    copyToIntVector(mxGetJc(M_ptr), n_M+1, Mp);
    int nnz_M = Mp[n_M]; // number of nonzeros
    
    // free memory
    mxFree(Mp);

    // return number of nonzeros
    return nnz_M;
}

int MEXSET::getVectorLength(const mxArray* x_ptr)
{
    // get size of array
    const mwSize * x_dim_ptr = mxGetDimensions(x_ptr);
    int m_x = (int) x_dim_ptr[0];
    int n_x = (int) x_dim_ptr[1];
    if (m_x > n_x)
        n_x = m_x; // in case input comes in as a col vector
    
    return n_x;
}

void MEXSET::buildSparseTripletsFromJcIrPr(int * jc, int * ir, double * pr, int n_jc, 
    std::vector<Eigen::Triplet<double>> &tripvec)
{
    // declare
    int row;
    double val;

    // loop
    for (int col=0; col<n_jc-1; col++)
    {
       for (int k=jc[col]; k<jc[col+1]; k++)
       {
        row = ir[k];
        val = pr[k];
        tripvec.push_back(Eigen::Triplet<double>(row, col, val));
       }
    }
}

void MEXSET::buildEigenMatrixFromDblPtr(double * M, int m, int n, Eigen::MatrixXd & M_out)
{
    for (int i=0; i<m; i++)
    {
        for (int j=0; j<n; j++)
        {
            M_out(i, j) = M[i + j*m];
        }
    }
}

void MEXSET::buildEigenVectorFromDblPtr(double * x, int numel, Eigen::VectorXd & x_out)
{
    for (int i=0; i<numel; i++)
        x_out[i] = x[i];
}

void MEXSET::copyEigenVector2DblPtr(const Eigen::VectorXd& vec, double * arr)
{
    for (int i=0; i<vec.rows(); i++)
        arr[i] = (double) vec[i];
}

void MEXSET::copyEigenMatrix2DblPtr(const Eigen::MatrixXd& mat, double * arr)
{
    for (int i=0; i<mat.rows(); i++)
    {
        for (int j=0; j<mat.cols(); j++)
        {
            arr[i + j*mat.rows()] = (double) mat(i, j);
        }
    }
}

// adapting from osqp-matlab on github
void MEXSET::copyToDoubleVector(double* vecData, int numel, double* out)
{
    //copy data
    for(int i=0; i < numel; i++)
        out[i] = (double)vecData[i];
}

void MEXSET::copyToIntVector(mwIndex* vecData, int numel, int* out)
{
    //copy data
    for(int i=0; i < numel; i++)
        out[i] = (int)vecData[i];
}
