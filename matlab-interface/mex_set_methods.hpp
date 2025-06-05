#ifndef _MEX_SET_METHODS_HPP_
#define _MEX_SET_METHODS_HPP_

#include "Eigen/Sparse"
#include "Eigen/Dense"
#include "mex.h"
#include <map>
#include <vector>

namespace MEXSET
{

// set methods
void getFullMatrixInput(const mxArray* M_ptr, Eigen::MatrixXd & M_out);
void getFullMatrixMap(const mxArray* M_i_cell_ptr, const mxArray* ind_M_i_ptr, std::map<int, Eigen::MatrixXd> &M_i_vec);
void getSparseMatrixMap(const mxArray* M_i_cell_ptr, const mxArray* ind_M_i_ptr, std::map<int, Eigen::SparseMatrix<double>> &M_i_vec);
void getVectorMap(const mxArray* x_i_cell_ptr, const mxArray* ind_x_i_ptr, std::map<int, Eigen::VectorXd> &x_i_vec);
void getVectorInput(const mxArray* x_ptr, Eigen::VectorXd &x_out);
void getSparseMatrixInput(const mxArray* M_ptr, std::vector<Eigen::Triplet<double>> &tripvec, Eigen::SparseMatrix<double> * M_out);
void getIntVectorInput(const mxArray* x_ptr, std::vector<int> &x_out);
void getIntVectorMap(const mxArray* x_i_cell_ptr, const mxArray* ind_x_i_ptr, std::map<int, std::vector<int>> &x_i_vec);
void getIntVectorInput_0_1_corrected(const mxArray* x_ptr, std::vector<int> &x_out);
void getIntVectorMap_0_1_corrected(const mxArray* x_i_cell_ptr, const mxArray* ind_x_i_ptr, std::map<int, std::vector<int>> &x_i_vec);
void R_x02region_processing(const mxArray* R_mat, std::map<int, std::vector<int>>& R_x02region);
void R_region2region_processing(const mxArray* R_mat, std::map<int, std::map<int, std::vector<int>>>& R_region2region);

// utilities
int getNNZ(const mxArray* M_ptr);
int getVectorLength(const mxArray* x_ptr);
void buildSparseTripletsFromJcIrPr(int * jc, int * ir, double * pr, int n_jc, 
    std::vector<Eigen::Triplet<double>> &tripvec);
void buildEigenMatrixFromDblPtr(double * M, int m, int n, Eigen::MatrixXd & M_out);
void buildEigenVectorFromDblPtr(double * x, int numel, Eigen::VectorXd & x_out);
void copyEigenVector2DblPtr(const Eigen::VectorXd & vec, double * arr);
void copyEigenMatrix2DblPtr(const Eigen::MatrixXd & mat, double * arr);
void copyToDoubleVector(double* vecData, int numel, double* out);
void copyToIntVector(mwIndex* vecData, int numel, int* out);

template <typename T>
void copyStdVector2DblPtr(const std::vector<T> & vec, double * arr)
{
    for (int i = 0; i < vec.size(); i++)
    {
        arr[i] = (double) vec[i];
    }
}

} // end namespace MEXSET

#endif