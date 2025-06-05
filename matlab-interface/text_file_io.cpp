#include "text_file_io.hpp"

void TextIO::int_2_txt(int val, std::ofstream& stream)
{
    stream << val << std::endl;
}

void TextIO::txt_2_int(std::ifstream& stream, int& val)
{
    stream >> val;
}

void TextIO::eigen_mat_2_txt(const Eigen::MatrixXd& mat, std::ofstream& stream)
{
    stream << mat.rows() << " " << mat.cols() << std::endl;
    for (int i = 0; i < mat.rows(); i++)
    {
        for (int j = 0; j < mat.cols(); j++)
        {
            stream << mat(i, j) << " ";
        }
        stream << std::endl;
    }
}

void TextIO::txt_2_eigen_mat(std::ifstream& stream, Eigen::MatrixXd& mat)
{
    int rows, cols;
    stream >> rows >> cols;
    mat.resize(rows, cols);
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            stream >> mat(i, j);
        }
    }
}

void TextIO::eigen_sparse_2_txt(const Eigen::SparseMatrix<double>& mat, std::ofstream& stream)
{
    stream << mat.rows() << " " << mat.cols() << " " << mat.nonZeros() << std::endl;
    for (int k = 0; k < mat.outerSize(); ++k)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator it(mat, k); it; ++it)
        {
            stream << it.row() << " " << it.col() << " " << it.value() << std::endl;
        }
    }
}

void TextIO::txt_2_eigen_sparse(std::ifstream& stream, Eigen::SparseMatrix<double>& mat)
{
    int rows, cols, nnz;
    stream >> rows >> cols >> nnz;
    mat.resize(rows, cols);
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(nnz);
    int row, col;
    double value;
    int cnt = 0;
    while ((cnt < nnz) && (stream >> row >> col >> value))
    {
        triplets.push_back(Eigen::Triplet<double>(row, col, value));
        cnt++;
    }
    mat.setFromTriplets(triplets.begin(), triplets.end());
}

void TextIO::eigen_vec_2_txt(const Eigen::VectorXd& vec, std::ofstream& stream)
{
    stream << vec.size() << std::endl;
    for (int i = 0; i < vec.size(); i++)
    {
        stream << vec(i) << " ";
    }
    stream << std::endl;
}

void TextIO::txt_2_eigen_vec(std::ifstream& stream, Eigen::VectorXd& vec)
{
    int size;
    stream >> size;
    vec.resize(size);
    for (int i = 0; i < size; i++)
    {
        stream >> vec(i);
    }
}

void TextIO::sparse_map_2_txt(const std::map<int, Eigen::SparseMatrix<double>>& sparse_map, std::ofstream& stream)
{
    stream << sparse_map.size() << std::endl;
    for (const auto& pair : sparse_map)
    {
        stream << pair.first << std::endl;
        eigen_sparse_2_txt(pair.second, stream);
    }
}

void TextIO::txt_2_sparse_map(std::ifstream& stream, std::map<int, Eigen::SparseMatrix<double>>& sparse_map)
{
    int size, key;
    Eigen::SparseMatrix<double> mat;
    stream >> size;
    for (int i = 0; i < size; i++)
    {
        stream >> key;
        txt_2_eigen_sparse(stream, mat);
        sparse_map[key] = mat;
    }
}

void TextIO::vec_map_2_txt(const std::map<int, Eigen::VectorXd>& vec_map, std::ofstream& stream)
{
    stream << vec_map.size() << std::endl;
    for (const auto& pair : vec_map)
    {
        stream << pair.first << std::endl;
        eigen_vec_2_txt(pair.second, stream);
    }
}

void TextIO::txt_2_vec_map(std::ifstream& stream, std::map<int, Eigen::VectorXd>& vec_map)
{
    int size, key;
    Eigen::VectorXd vec;
    stream >> size;
    for (int i = 0; i < size; i++)
    {
        stream >> key;
        txt_2_eigen_vec(stream, vec);
        vec_map[key] = vec;
    }
}

void TextIO::int_vec_2_txt(const std::vector<int>& vec, std::ofstream& stream)
{
    stream << vec.size() << std::endl;
    for (auto it = vec.begin(); it != vec.end(); it++)
    {
        stream << *it << " ";
    }
    stream << std::endl;
}

void TextIO::txt_2_int_vec(std::ifstream& stream, std::vector<int>& vec)
{
    int size;
    stream >> size;
    vec.clear();
    int val;
    int cnt = 0;
    while ((cnt < size) && (stream >> val))
    {
        vec.push_back(val);
        cnt++;
    }
}

void TextIO::int_map_2_txt(const std::map<int, std::vector<int>>& int_map, std::ofstream& stream)
{
    stream << int_map.size() << std::endl;
    for (const auto& pair : int_map)
    {
        stream << pair.first << std::endl;
        int_vec_2_txt(pair.second, stream);
    }
}

void TextIO::txt_2_int_map(std::ifstream& stream, std::map<int, std::vector<int>>& int_map)
{
    int size, key;
    std::vector<int> vec;
    stream >> size;
    for (int i = 0; i < size; i++)
    {
        stream >> key;
        txt_2_int_vec(stream, vec);
        int_map[key] = vec;
    }
}

void TextIO::qpipmpc_2_file(const MI_QPIPMPC::QPIPMPC_Data& qpipmpc_data, const std::string& filename)
{
    std::ofstream file_stream;
    file_stream.open(filename);

    // write data using methods in text_file_io
    eigen_sparse_2_txt(qpipmpc_data.P_i_nom, file_stream);
    eigen_vec_2_txt(qpipmpc_data.q_i_nom, file_stream);
    sparse_map_2_txt(qpipmpc_data.P_i_vec, file_stream);
    vec_map_2_txt(qpipmpc_data.q_i_vec, file_stream);
    eigen_sparse_2_txt(qpipmpc_data.C_i_nom, file_stream);
    eigen_sparse_2_txt(qpipmpc_data.D_i_nom, file_stream);
    eigen_vec_2_txt(qpipmpc_data.crhs_i_nom, file_stream);
    sparse_map_2_txt(qpipmpc_data.C_i_vec, file_stream);
    sparse_map_2_txt(qpipmpc_data.D_i_vec, file_stream);
    vec_map_2_txt(qpipmpc_data.crhs_i_vec, file_stream);
    eigen_sparse_2_txt(qpipmpc_data.G_i_nom, file_stream);
    eigen_vec_2_txt(qpipmpc_data.w_i_nom, file_stream);
    sparse_map_2_txt(qpipmpc_data.G_i_vec, file_stream);
    vec_map_2_txt(qpipmpc_data.w_i_vec, file_stream);
    int_vec_2_txt(qpipmpc_data.idx_i_nom, file_stream);
    int_map_2_txt(qpipmpc_data.idx_i_vec, file_stream);
    int_map_2_txt(qpipmpc_data.idx_state, file_stream);
    int_map_2_txt(qpipmpc_data.idx_input, file_stream);
    int_map_2_txt(qpipmpc_data.idx_binvar, file_stream);
    int_map_2_txt(qpipmpc_data.idx_y, file_stream);
    int_map_2_txt(qpipmpc_data.idx_eq, file_stream);
    int_map_2_txt(qpipmpc_data.idx_ineq, file_stream);
    eigen_mat_2_txt(qpipmpc_data.Pp, file_stream);
    int_2_txt(qpipmpc_data.n_horizon, file_stream);

    // close file
    file_stream.close();
}

void TextIO::file_2_qpipmpc(const std::string& filename, MI_QPIPMPC::QPIPMPC_Data& qpipmpc_data)
{
    std::ifstream file_stream;
    file_stream.open(filename);

    // read data using methods in text_file_io
    txt_2_eigen_sparse(file_stream, qpipmpc_data.P_i_nom);
    txt_2_eigen_vec(file_stream, qpipmpc_data.q_i_nom);
    txt_2_sparse_map(file_stream, qpipmpc_data.P_i_vec);
    txt_2_vec_map(file_stream, qpipmpc_data.q_i_vec);
    txt_2_eigen_sparse(file_stream, qpipmpc_data.C_i_nom);
    txt_2_eigen_sparse(file_stream, qpipmpc_data.D_i_nom);
    txt_2_eigen_vec(file_stream, qpipmpc_data.crhs_i_nom);
    txt_2_sparse_map(file_stream, qpipmpc_data.C_i_vec);
    txt_2_sparse_map(file_stream, qpipmpc_data.D_i_vec);
    txt_2_vec_map(file_stream, qpipmpc_data.crhs_i_vec);
    txt_2_eigen_sparse(file_stream, qpipmpc_data.G_i_nom);
    txt_2_eigen_vec(file_stream, qpipmpc_data.w_i_nom);
    txt_2_sparse_map(file_stream, qpipmpc_data.G_i_vec);
    txt_2_vec_map(file_stream, qpipmpc_data.w_i_vec);
    txt_2_int_vec(file_stream, qpipmpc_data.idx_i_nom);
    txt_2_int_map(file_stream, qpipmpc_data.idx_i_vec);
    txt_2_int_map(file_stream, qpipmpc_data.idx_state);
    txt_2_int_map(file_stream, qpipmpc_data.idx_input);
    txt_2_int_map(file_stream, qpipmpc_data.idx_binvar);
    txt_2_int_map(file_stream, qpipmpc_data.idx_y);
    txt_2_int_map(file_stream, qpipmpc_data.idx_eq);
    txt_2_int_map(file_stream, qpipmpc_data.idx_ineq);
    txt_2_eigen_mat(file_stream, qpipmpc_data.Pp);
    txt_2_int(file_stream, qpipmpc_data.n_horizon);

    // close file
    file_stream.close();
}

void TextIO::mi_settings_2_file(const MI_QPIPMPC::MI_QPIPMPC_Settings& mi_settings, const std::string& filename)
{
    std::ofstream file_stream;
    file_stream.open(filename);

    // copy settings to file
    file_stream << mi_settings.eps_feas << " ";
    file_stream << mi_settings.max_iter_bb << " ";
    file_stream << mi_settings.verbose << " ";
    file_stream << mi_settings.conv_rel << " ";
    file_stream << mi_settings.conv_abs << " ";
    file_stream << mi_settings.T_max << " ";
    file_stream << mi_settings.n_threads << std::endl;

    // close file
    file_stream.close();
}

void TextIO::file_2_mi_settings(const std::string& filename, MI_QPIPMPC::MI_QPIPMPC_Settings& mi_settings)
{
    std::ifstream file_stream;
    file_stream.open(filename);

    // read settings from file
    file_stream >> mi_settings.eps_feas;
    file_stream >> mi_settings.max_iter_bb;
    file_stream >> mi_settings.verbose;
    file_stream >> mi_settings.conv_rel;
    file_stream >> mi_settings.conv_abs;
    file_stream >> mi_settings.T_max;
    file_stream >> mi_settings.n_threads;

    // close file
    file_stream.close();
}

void TextIO::qp_settings_2_file(const QP_IP_MPC::QP_settings& qp_settings, const std::string& filename)
{
    std::ofstream file_stream;
    file_stream.open(filename);

    // copy settings to file
    file_stream << qp_settings.mu_term << " ";
    file_stream << qp_settings.mu_feas << " ";
    file_stream << qp_settings.mu_max << " ";
    file_stream << qp_settings.mu_init << " ";
    file_stream << qp_settings.iter_max << " ";
    file_stream << qp_settings.gamma << " ";
    file_stream << qp_settings.auto_warm_start << " ";
    file_stream << qp_settings.t_ls << " ";
    file_stream << qp_settings.T_max << std::endl;

    // close file
    file_stream.close();
}

void TextIO::file_2_qp_settings(const std::string& filename, QP_IP_MPC::QP_settings& qp_settings)
{
    std::ifstream file_stream;
    file_stream.open(filename);

    // read settings from file
    file_stream >> qp_settings.mu_term;
    file_stream >> qp_settings.mu_feas;
    file_stream >> qp_settings.mu_max;
    file_stream >> qp_settings.mu_init;
    file_stream >> qp_settings.iter_max;
    file_stream >> qp_settings.gamma;
    file_stream >> qp_settings.auto_warm_start;
    file_stream >> qp_settings.t_ls;
    file_stream >> qp_settings.T_max;

    // close file
    file_stream.close();
}

void TextIO::R_region2region_2_file(const std::map<int, std::map<int, std::vector<int>>>& R_region2region, const std::string& filename)
{
    std::ofstream file_stream;
    file_stream.open(filename);

    // write data using methods in text_file_io
    file_stream << R_region2region.size() << std::endl;
    for (const auto& pair : R_region2region)
    {
        file_stream << pair.first << std::endl;
        int_map_2_txt(pair.second, file_stream);
    }

    // close file
    file_stream.close();
}

void TextIO::file_2_R_region2region(const std::string& filename, std::map<int, std::map<int, std::vector<int>>>& R_region2region)
{
    std::ifstream file_stream;
    file_stream.open(filename);

    // read data using methods in text_file_io
    int size, key;
    std::map<int, std::vector<int>> int_map;
    file_stream >> size;
    for (int i = 0; i < size; i++)
    {
        file_stream >> key;
        txt_2_int_map(file_stream, int_map);
        R_region2region[key] = int_map;
    }

    // close file
    file_stream.close();
}

void TextIO::R_x02region_2_file(const std::map<int, std::vector<int>>& R_x02region, const std::string& filename)
{
    std::ofstream file_stream;
    file_stream.open(filename);

    // write data using methods in text_file_io
    int_map_2_txt(R_x02region, file_stream);

    // close file
    file_stream.close();
}

void TextIO::file_2_R_x02region(const std::string& filename, std::map<int, std::vector<int>>& R_x02region)
{
    std::ifstream file_stream;
    file_stream.open(filename);

    // read data using methods in text_file_io
    txt_2_int_map(file_stream, R_x02region);

    // close file
    file_stream.close();
}

void TextIO::int_2_file(int val, const std::string& filename)
{
    std::ofstream file_stream;
    file_stream.open(filename);

    // write data using methods in text_file_io
    file_stream << val << std::endl;

    // close file
    file_stream.close();
}

void TextIO::file_2_int(const std::string& filename, int& val)
{
    std::ifstream file_stream;
    file_stream.open(filename);

    // read data using methods in text_file_io
    file_stream >> val;

    // close file
    file_stream.close();
}

void TextIO::mat_map_2_file(const std::map<int, Eigen::MatrixXd>& mat_map, const std::string& filename)
{
    std::ofstream file_stream;
    file_stream.open(filename);

    // write data using methods in text_file_io
    file_stream << mat_map.size() << std::endl;
    for (const auto& pair : mat_map)
    {
        file_stream << pair.first << std::endl;
        eigen_mat_2_txt(pair.second, file_stream);
    }

    // close file
    file_stream.close();
}

void TextIO::file_2_mat_map(const std::string& filename, std::map<int, Eigen::MatrixXd>& mat_map)
{
    std::ifstream file_stream;
    file_stream.open(filename);

    // read data using methods in text_file_io
    int size, key;
    Eigen::MatrixXd mat;
    file_stream >> size;
    for (int i = 0; i < size; i++)
    {
        file_stream >> key;
        txt_2_eigen_mat(file_stream, mat);
        mat_map[key] = mat;
    }

    // close file
    file_stream.close();
}

void TextIO::vec_map_2_file(const std::map<int, Eigen::VectorXd>& vec_map, const std::string& filename)
{
    std::ofstream file_stream;
    file_stream.open(filename);

    vec_map_2_txt(vec_map, file_stream);

    // close file
    file_stream.close();
}

void TextIO::file_2_vec_map(const std::string& filename, std::map<int, Eigen::VectorXd>& vec_map)
{
    std::ifstream file_stream;
    file_stream.open(filename);

    txt_2_vec_map(file_stream, vec_map);

    // close file
    file_stream.close();
}

void TextIO::sparse_mat_2_file(const Eigen::SparseMatrix<double>& mat, const std::string& filename)
{
    std::ofstream file_stream;
    file_stream.open(filename);

    eigen_sparse_2_txt(mat, file_stream);

    // close file
    file_stream.close();
}

void TextIO::file_2_sparse_mat(const std::string& filename, Eigen::SparseMatrix<double>& mat)
{
    std::ifstream file_stream;
    file_stream.open(filename);

    txt_2_eigen_sparse(file_stream, mat);

    // close file
    file_stream.close();
}

void TextIO::int_vec_2_file(const std::vector<int>& vec, const std::string& filename)
{
    std::ofstream file_stream;
    file_stream.open(filename);

    int_vec_2_txt(vec, file_stream);

    // close file
    file_stream.close();
}

void TextIO::file_2_int_vec(const std::string& filename, std::vector<int>& vec)
{
    std::ifstream file_stream;
    file_stream.open(filename);

    txt_2_int_vec(file_stream, vec);

    // close file
    file_stream.close();
}

void TextIO::int_map_2_file(const std::map<int, std::vector<int>>& int_map, const std::string& filename)
{
    std::ofstream file_stream;
    file_stream.open(filename);

    int_map_2_txt(int_map, file_stream);

    // close file
    file_stream.close();
}

void TextIO::file_2_int_map(const std::string& filename, std::map<int, std::vector<int>>& int_map)
{
    std::ifstream file_stream;
    file_stream.open(filename);

    txt_2_int_map(file_stream, int_map);

    // close file
    file_stream.close();
}

void TextIO::vec_2_file(const Eigen::VectorXd& vec, const std::string& filename)
{
    std::ofstream file_stream;
    file_stream.open(filename);

    eigen_vec_2_txt(vec, file_stream);

    // close file
    file_stream.close();
}

void TextIO::file_2_vec(const std::string& filename, Eigen::VectorXd& vec)
{
    std::ifstream file_stream;
    file_stream.open(filename);

    txt_2_eigen_vec(file_stream, vec);

    // close file
    file_stream.close();
}

void TextIO::sparse_map_2_file(const std::map<int, Eigen::SparseMatrix<double>>& sparse_map, const std::string& filename)
{
    std::ofstream file_stream;
    file_stream.open(filename);

    sparse_map_2_txt(sparse_map, file_stream);

    // close file
    file_stream.close();
}

void TextIO::file_2_sparse_map(const std::string& filename, std::map<int, Eigen::SparseMatrix<double>>& sparse_map)
{
    std::ifstream file_stream;
    file_stream.open(filename);

    txt_2_sparse_map(file_stream, sparse_map);

    // close file
    file_stream.close();
}