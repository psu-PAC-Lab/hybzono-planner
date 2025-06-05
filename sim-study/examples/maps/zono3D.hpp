#ifndef _ZONO3D_HPP_
#define _ZONO3D_HPP_

#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "ZonoCpp.hpp"
#include <map>
#include <utility>
#include <limits>
#include <array>

using namespace HybZonoMPC;

void zono3D(ZonoCpp::HybZono &X, std::map<int, Eigen::MatrixXd> &A_Hrep, std::map<int, Eigen::VectorXd> &b_Hrep,
    Eigen::VectorXd &pos_0, Eigen::VectorXd &pos_ref, double &T_sim, Eigen::MatrixXd &d_mat,
    std::vector<Eigen::MatrixXd> &V_obs, std::vector<Eigen::MatrixXd> &V_free, Eigen::MatrixXd &V_lim)
{
    // range
    std::pair<double,double> x_range = std::make_pair(-15, 15);
    std::pair<double,double> y_range = std::make_pair(-10, 10);
    std::pair<double,double> z_range = std::make_pair(-5, 5);    

    // initial / final positions
    pos_0.resize(3);
    pos_0 << -15, 0, 0;
    pos_ref.resize(3);
    pos_ref << 15, 0, 0;

    // sim time
    T_sim = 40;

    // grid size
    int m_grid = 7;
    int n_grid = 3;
    int l_grid = 3;

    // obstacles
    std::vector<int> i_obs_grid = {2, 2, 2, 2, 2, 2, 2, 4, 4, 4, 4, 4, 4, 4, 6, 6, 6, 6, 6, 6, 6};
    std::vector<int> j_obs_grid = {1, 1, 1, 2, 2, 3, 3, 1, 1, 1, 2, 3, 3, 3, 1, 1, 1, 2, 2, 2, 3};
    std::vector<int> k_obs_grid = {1, 2, 3, 1, 2, 1, 2, 1, 2, 3, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1};

    // invert wrt x
    for (int i=0; i<i_obs_grid.size(); i++)
    {
        i_obs_grid[i] = m_grid - i_obs_grid[i] + 1;
    }

    // free space grid
    std::vector<std::tuple<int, int, int>> free_space_grid;
    for (int i=1; i<=m_grid; i++)
    {
        for (int j=1; j<=n_grid; j++)
        {
            for (int k=1; k<=l_grid; k++)
            {
                for (int ind_obs=0; ind_obs<i_obs_grid.size(); ind_obs++)
                {
                    if (i_obs_grid[ind_obs] == i && j_obs_grid[ind_obs] == j && k_obs_grid[ind_obs] == k)
                        break;
                    if (ind_obs == i_obs_grid.size()-1)
                    {
                        free_space_grid.push_back({i-1, j-1, k-1});
                    }
                }
            }
        }
    }

    // grid dimensions
    double dx = (x_range.second - x_range.first)/m_grid;
    double dy = (y_range.second - y_range.first)/n_grid;
    double dz = (z_range.second - z_range.first)/l_grid;

    // get tile centers
    std::vector<double> xc_tiles, yc_tiles, zc_tiles;
    for (int i=0; i<free_space_grid.size(); i++)
    {
        xc_tiles.push_back(x_range.first + (((double) std::get<0>(free_space_grid[i])) + 0.5)*dx);
        yc_tiles.push_back(y_range.first + (((double) std::get<1>(free_space_grid[i])) + 0.5)*dy);
        zc_tiles.push_back(z_range.first + (((double) std::get<2>(free_space_grid[i])) + 0.5)*dz);
    }

    // directly construct hybzono

    // generator
    Eigen::Matrix<double, 3, 3> Gc = Eigen::Matrix<double, 3, 3>::Zero();
    Gc.diagonal() << dx, dy, dz;
    Eigen::MatrixXd Gb = Eigen::MatrixXd::Zero(3, xc_tiles.size());
    for (int i=0; i<xc_tiles.size(); i++)
    {
        Gb(0, i) = xc_tiles[i];
        Gb(1, i) = yc_tiles[i];
        Gb(2, i) = zc_tiles[i];
    }
    Eigen::Vector<double, 3> c;
    c << -dx/2, -dy/2, -dz/2;

    // equality constraints
    Eigen::MatrixXd Ac = Eigen::MatrixXd::Zero(1, 3);
    Eigen::MatrixXd Ab = Eigen::MatrixXd::Ones(1, xc_tiles.size());
    Eigen::Vector<double, 1> b = Eigen::Vector<double, 1>::Ones();

    // set hybzono
    X.set(Gc.sparseView(), Gb.sparseView(), c, Ac.sparseView(), Ab.sparseView(), b, true);

    // get H-rep
    A_Hrep.clear();
    b_Hrep.clear();

    // polytope for template cell
    Eigen::Matrix<double, 6, 3> A_temp;
    A_temp << 1, 0, 0,
             -1, 0, 0,
              0, 1, 0,
              0, -1, 0,
              0, 0, 1,
              0, 0, -1;
    Eigen::Vector<double, 6> b_temp;
    b_temp << dx/2, dx/2, dy/2, dy/2, dz/2, dz/2;

    // get H-rep for each tile
    Eigen::Vector<double, 3> c_tile;
    for (int i=0; i<xc_tiles.size(); i++)
    {
        c_tile(0) = xc_tiles[i];
        c_tile(1) = yc_tiles[i];
        c_tile(2) = zc_tiles[i];
        A_Hrep[i+1] = A_temp; // 1-based indexing
        b_Hrep[i+1] = b_temp + A_temp*c_tile; // 1-based indexing
    }

    // calculate distance between regions
    std::array<std::tuple<double, double, double>, 8> delta_corners = {
        std::make_tuple(-dx/2, -dy/2, -dz/2),
        std::make_tuple(-dx/2, dy/2, -dz/2),
        std::make_tuple(dx/2, dy/2, -dz/2),
        std::make_tuple(dx/2, -dy/2, -dz/2),
        std::make_tuple(-dx/2, -dy/2, dz/2),
        std::make_tuple(-dx/2, dy/2, dz/2),
        std::make_tuple(dx/2, dy/2, dz/2),
        std::make_tuple(dx/2, -dy/2, dz/2)
    };

    d_mat = Eigen::MatrixXd::Zero(xc_tiles.size(), xc_tiles.size()); // init
    double dist, min_dist, dx2, dy2, dz2; // declare
    for (int i=0; i<xc_tiles.size()-1; i++)
    {
        for (int j=i+1; j<xc_tiles.size(); j++)
        {   
            min_dist = std::numeric_limits<double>::max(); // init
            for (int ii=0; ii<8; ii++)
            {
                for (int jj=0; jj<8; jj++)
                {
                    dx2 = std::pow(xc_tiles[i] + std::get<0>(delta_corners[ii]) - 
                        (xc_tiles[j] + std::get<0>(delta_corners[jj])), 2);
                    dy2 = std::pow(yc_tiles[i] + std::get<1>(delta_corners[ii]) - 
                        (yc_tiles[j] + std::get<1>(delta_corners[jj])), 2);
                    dz2 = std::pow(zc_tiles[i] + std::get<2>(delta_corners[ii]) - 
                        (zc_tiles[j] + std::get<2>(delta_corners[jj])), 2);

                    dist = std::sqrt(dx2 + dy2 + dz2);
                    if (dist < min_dist)
                        min_dist = dist;
                }
            }
            d_mat(i,j) = min_dist;
            d_mat(j,i) = min_dist; // symmetric
        }
    }

    // get obstacle grid
    std::vector<std::tuple<int, int, int>> obs_grid;
    for (int i=0; i<m_grid; i++)
    {
        for (int j=0; j<n_grid; j++)
        {
            for (int k=0; k<l_grid; k++)
            {
                if (std::find(free_space_grid.begin(), free_space_grid.end(), std::make_tuple(i, j, k)) == free_space_grid.end())
                    obs_grid.push_back({i, j, k});
            }
        }
    }

    // cell vertices relative to center
    Eigen::Matrix<double, 8, 3> dVr, Vtmp;
    dVr << -dx/2, -dy/2, -dz/2,
           dx/2, -dy/2, -dz/2,
           dx/2, dy/2, -dz/2,
           -dx/2, dy/2, -dz/2,
           -dx/2, -dy/2, dz/2,
           dx/2, -dy/2, dz/2,
           dx/2, dy/2, dz/2,
           -dx/2, dy/2, dz/2;

    // center
    Eigen::Matrix<double, 1, 3> Vc;

    // obstacle vertices
    V_obs.clear();
    for (auto it=obs_grid.begin(); it!=obs_grid.end(); ++it)
    {
        Vc(0) = x_range.first + (((double) std::get<0>(*it))+0.5)*dx;
        Vc(1) = y_range.first + (((double) std::get<1>(*it))+0.5)*dy;
        Vc(2) = z_range.first + (((double) std::get<2>(*it))+0.5)*dz;
        Vtmp = Vc.replicate(8, 1) + dVr;
        V_obs.push_back(Vtmp);
    }

    // free space vertices
    V_free.clear();
    for (auto it=free_space_grid.begin(); it!=free_space_grid.end(); ++it)
    {
        Vc(0) = x_range.first + (((double) std::get<0>(*it))+0.5)*dx;
        Vc(1) = y_range.first + (((double) std::get<1>(*it))+0.5)*dy;
        Vc(2) = z_range.first + (((double) std::get<2>(*it))+0.5)*dz;
        Vtmp = Vc.replicate(8, 1) + dVr;
        V_free.push_back(Vtmp);
    }

    // limit vertices
    V_lim.resize(8, 3);
    V_lim << x_range.first, y_range.first, z_range.first,
             x_range.second, y_range.first, z_range.first,
             x_range.second, y_range.second, z_range.first,
             x_range.first, y_range.second, z_range.first,
             x_range.first, y_range.second, z_range.second,
             x_range.second, y_range.second, z_range.second,
             x_range.second, y_range.first, z_range.second,
             x_range.first, y_range.first, z_range.second;
}

#endif