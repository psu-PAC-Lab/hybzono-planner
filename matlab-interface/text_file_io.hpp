#ifndef _TEXT_FILE_IO_HPP_
#define _TEXT_FILE_IO_HPP_

#include <iostream>
#include <fstream>
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include <map>
#include <vector>
#include "MI_MPC.hpp"
#include <string>

namespace TextIO
{

// primitive functions
void int_2_txt(int val, std::ofstream& stream);
void txt_2_int(std::ifstream& stream, int& val);
void eigen_mat_2_txt(const Eigen::MatrixXd& mat, std::ofstream& stream);
void txt_2_eigen_mat(std::ifstream& stream, Eigen::MatrixXd& mat);
void eigen_sparse_2_txt(const Eigen::SparseMatrix<double>& mat, std::ofstream& stream);
void txt_2_eigen_sparse(std::ifstream& stream, Eigen::SparseMatrix<double>& mat);
void eigen_vec_2_txt(const Eigen::VectorXd& vec, std::ofstream& stream);
void txt_2_eigen_vec(std::ifstream& stream, Eigen::VectorXd& vec);
void sparse_map_2_txt(const std::map<int, Eigen::SparseMatrix<double>>& sparse_map, std::ofstream& stream);
void txt_2_sparse_map(std::ifstream& stream, std::map<int, Eigen::SparseMatrix<double>>& sparse_map);
void vec_map_2_txt(const std::map<int, Eigen::VectorXd>& vec_map, std::ofstream& stream);
void txt_2_vec_map(std::ifstream& stream, std::map<int, Eigen::VectorXd>& vec_map);
void int_vec_2_txt(const std::vector<int>& vec, std::ofstream& stream);
void txt_2_int_vec(std::ifstream& stream, std::vector<int>& vec);
void int_map_2_txt(const std::map<int, std::vector<int>>& int_map, std::ofstream& stream);
void txt_2_int_map(std::ifstream& stream, std::map<int, std::vector<int>>& int_map);

// MI_MPC functions
void qpipmpc_2_file(const MI_QPIPMPC::QPIPMPC_Data& qpipmpc_data, const std::string& filename);
void file_2_qpipmpc(const std::string& filename, MI_QPIPMPC::QPIPMPC_Data& qpipmpc_data);
void mi_settings_2_file(const MI_QPIPMPC::MI_QPIPMPC_Settings& mi_settings, const std::string& filename);
void file_2_mi_settings(const std::string& filename, MI_QPIPMPC::MI_QPIPMPC_Settings& mi_settings);
void qp_settings_2_file(const QP_IP_MPC::QP_settings& qp_settings, const std::string& filename);
void file_2_qp_settings(const std::string& filename, QP_IP_MPC::QP_settings& qp_settings);
void R_region2region_2_file(const std::map<int, std::map<int, std::vector<int>>>& R_region2region, const std::string& filename);
void file_2_R_region2region(const std::string& filename, std::map<int, std::map<int, std::vector<int>>>& R_region2region);
void R_x02region_2_file(const std::map<int, std::vector<int>>& R_x02region, const std::string& filename);
void file_2_R_x02region(const std::string& filename, std::map<int, std::vector<int>>& R_x02region);
void int_2_file(int val, const std::string& filename);
void file_2_int(const std::string& filename, int& val);
void mat_map_2_file(const std::map<int, Eigen::MatrixXd>& mat_map, const std::string& filename);
void file_2_mat_map(const std::string& filename, std::map<int, Eigen::MatrixXd>& mat_map);
void vec_map_2_file(const std::map<int, Eigen::VectorXd>& vec_map, const std::string& filename);
void file_2_vec_map(const std::string& filename, std::map<int, Eigen::VectorXd>& vec_map);
void sparse_mat_2_file(const Eigen::SparseMatrix<double>& mat, const std::string& filename);
void file_2_sparse_mat(const std::string& filename, Eigen::SparseMatrix<double>& mat);
void int_vec_2_file(const std::vector<int>& vec, const std::string& filename);
void file_2_int_vec(const std::string& filename, std::vector<int>& vec);
void int_map_2_file(const std::map<int, std::vector<int>>& int_map, const std::string& filename);
void file_2_int_map(const std::string& filename, std::map<int, std::vector<int>>& int_map);
void vec_2_file(const Eigen::VectorXd& vec, const std::string& filename);
void file_2_vec(const std::string& filename, Eigen::VectorXd& vec);
void sparse_map_2_file(const std::map<int, Eigen::SparseMatrix<double>>& sparse_map, const std::string& filename);
void file_2_sparse_map(const std::string& filename, std::map<int, Eigen::SparseMatrix<double>>& sparse_map);


} // end namespace TextIO

#endif // _TEXT_FILE_IO_HPP_