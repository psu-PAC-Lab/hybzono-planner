#ifndef _ZONOPSU_HPP_
#define _ZONOPSU_HPP_

#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "ZonoCpp.hpp"
#include <map>
#include <utility>
#include <limits>
#include <array>

using namespace HybZonoMPC;

void zonopsu(ZonoCpp::HybZono &X, std::map<int, Eigen::MatrixXd> &A_Hrep, std::map<int, Eigen::VectorXd> &b_Hrep,
    Eigen::VectorXd &pos_0, Eigen::VectorXd &pos_ref, double &T_sim, Eigen::MatrixXd &d_mat,
    std::vector<Eigen::MatrixXd> &V_obs, std::vector<Eigen::MatrixXd> &V_free, Eigen::MatrixXd &V_lim)
{
    // range
    std::pair<double,double> x_range = std::make_pair(-37.5, 37.5);
    std::pair<double,double> y_range = std::make_pair(-10, 10);    

    // initial / final positions
    pos_0.resize(2);
    pos_0 << -37.5, 0;
    pos_ref.resize(2);
    pos_ref << 37.5, 0;

    // sim time
    T_sim = 80;

    // free space grid
    std::vector<std::vector<int>> free_space_grid;
    free_space_grid.push_back({2});
    free_space_grid.push_back({2});
    free_space_grid.push_back({2});
    free_space_grid.push_back({2});
    free_space_grid.push_back({2});
    free_space_grid.push_back({1, 2, 3});
    free_space_grid.push_back({1, 3});
    free_space_grid.push_back({1, 3});
    free_space_grid.push_back({1, 3});
    free_space_grid.push_back({1, 3});
    free_space_grid.push_back({2});
    free_space_grid.push_back({2});
    free_space_grid.push_back({2});
    free_space_grid.push_back({2});
    free_space_grid.push_back({1, 2, 3});
    free_space_grid.push_back({3});
    free_space_grid.push_back({1, 2, 3});
    free_space_grid.push_back({1});
    free_space_grid.push_back({1, 2, 3});
    free_space_grid.push_back({2});
    free_space_grid.push_back({2});
    free_space_grid.push_back({2});
    free_space_grid.push_back({2});
    free_space_grid.push_back({1, 2});
    free_space_grid.push_back({1});
    free_space_grid.push_back({1, 2, 3});
    free_space_grid.push_back({1, 3});
    free_space_grid.push_back({1, 2, 3});
    free_space_grid.push_back({2});
    free_space_grid.push_back({2});

    std::vector<int> i_grid, j_grid;
    for (int i=0; i<free_space_grid.size(); i++)
    {
        for (int j : free_space_grid[i])
        {
            i_grid.push_back(i);
            j_grid.push_back(j);
        }
    }

    // grid dimensions
    int m_grid = free_space_grid.size();
    int n_grid = 5;
    double dx = (x_range.second - x_range.first)/m_grid;
    double dy = (y_range.second - y_range.first)/n_grid;

    // get tile centers
    std::vector<double> xc_tiles, yc_tiles;
    for (int i=0; i<i_grid.size(); i++)
    {
        xc_tiles.push_back(x_range.first + (((double) i_grid[i]) + 0.5)*dx);
        yc_tiles.push_back(y_range.first + (((double) j_grid[i]) + 0.5)*dy);
    }

    // directly construct hybzono

    // generator
    Eigen::Matrix<double, 2, 2> Gc = Eigen::Matrix<double, 2, 2>::Zero();
    Gc.diagonal() << dx, dy;
    Eigen::MatrixXd Gb = Eigen::MatrixXd::Zero(2, xc_tiles.size());
    for (int i=0; i<xc_tiles.size(); i++)
    {
        Gb(0, i) = xc_tiles[i];
        Gb(1, i) = yc_tiles[i];
    }
    Eigen::Vector<double, 2> c;
    c << -dx/2, -dy/2;

    // equality constraints
    Eigen::MatrixXd Ac = Eigen::MatrixXd::Zero(1, 2);
    Eigen::MatrixXd Ab = Eigen::MatrixXd::Ones(1, xc_tiles.size());
    Eigen::Vector<double, 1> b = Eigen::Vector<double, 1>::Ones();

    // set hybzono
    X.set(Gc.sparseView(), Gb.sparseView(), c, Ac.sparseView(), Ab.sparseView(), b, true);

    // get H-rep
    A_Hrep.clear();
    b_Hrep.clear();

    // polytope for template cell
    Eigen::Matrix<double, 4, 2> A_temp;
    A_temp << 1, 0,
             -1, 0,
              0, 1,
              0, -1;
    Eigen::Vector<double, 4> b_temp;
    b_temp << dx/2, dx/2, dy/2, dy/2;

    // get H-rep for each tile
    Eigen::Vector<double, 2> c_tile;
    for (int i=0; i<xc_tiles.size(); i++)
    {
        c_tile(0) = xc_tiles[i];
        c_tile(1) = yc_tiles[i];
        A_Hrep[i+1] = A_temp; // 1-based indexing
        b_Hrep[i+1] = b_temp + A_temp*c_tile; // 1-based indexing
    }

    // calculate distance between regions
    std::array<std::pair<double, double>, 4> delta_corners = {
        std::make_pair(-dx/2, -dy/2),
        std::make_pair(-dx/2, dy/2),
        std::make_pair(dx/2, dy/2),
        std::make_pair(dx/2, -dy/2)
    };

    d_mat = Eigen::MatrixXd::Zero(xc_tiles.size(), xc_tiles.size()); // init
    double dist, min_dist; // declare
    for (int i=0; i<xc_tiles.size()-1; i++)
    {
        for (int j=i+1; j<xc_tiles.size(); j++)
        {   
            min_dist = std::numeric_limits<double>::max(); // init
            for (int ii=0; ii<4; ii++)
            {
                for (int jj=0; jj<4; jj++)
                {
                    dist = std::sqrt(std::pow(xc_tiles[i] + delta_corners[ii].first - (xc_tiles[j] + delta_corners[jj].first), 2) + 
                        std::pow(yc_tiles[i] + delta_corners[ii].second - (yc_tiles[j] + delta_corners[jj].second), 2));
                    if (dist < min_dist)
                        min_dist = dist;
                }
            }
            d_mat(i,j) = min_dist;
            d_mat(j,i) = min_dist; // symmetric
        }
    }

    // switch free space grid format
    std::vector<std::pair<int,int>> free_space;
    for (int i=0; i<m_grid; i++)
    {
        for (auto it = free_space_grid[i].begin(); it != free_space_grid[i].end(); ++it)
        {
            free_space.push_back({i+1, *it+1}); // one-based indexing
        }
    }

    // get obstacle grid
    std::vector<std::pair<int,int>> obs_space;
    for (int i=0; i<m_grid; i++)
    {
        for (int j=0; j<n_grid; j++)
        {
            if (std::find(free_space.begin(), free_space.end(), std::make_pair(i+1, j+1)) == free_space.end())
                obs_space.push_back({i+1, j+1}); // one-based indexing
        }
    }

    // cell vertices relative to center
    Eigen::Matrix<double, 4, 2> dVr, Vtmp;
    dVr << -dx/2, -dy/2,
           dx/2, -dy/2,
           dx/2, dy/2,
           -dx/2, dy/2;

    // center
    Eigen::Matrix<double, 1, 2> Vc;

    // obstacle vertices
    V_obs.clear();
    for (auto it=obs_space.begin(); it!=obs_space.end(); ++it)
    {   
        Vc(0) = x_range.first + (((double) it->first)-0.5)*dx;
        Vc(1) = y_range.first + (((double) it->second)-0.5)*dy;
        Vtmp = Vc.replicate(4, 1) + dVr;
        V_obs.push_back(Vtmp);
    }

    // free space vertices
    V_free.clear();
    for (auto it=free_space.begin(); it!=free_space.end(); ++it)
    {
        Vc(0) = x_range.first + (((double) it->first)-0.5)*dx;
        Vc(1) = y_range.first + (((double) it->second)-0.5)*dy;
        Vtmp = Vc.replicate(4, 1) + dVr;
        V_free.push_back(Vtmp);
    }

    // limits
    V_lim.resize(4, 2);
    V_lim << x_range.first, y_range.first,
             x_range.second, y_range.first,
             x_range.second, y_range.second,
             x_range.first, y_range.second;
}

#endif