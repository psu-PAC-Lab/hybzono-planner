#ifndef _ZONOCOST1_HPP_
#define _ZONOCOST1_HPP_

#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "ZonoCpp.hpp"
#include <map>
#include <utility>
#include <limits>
#include <array>

using namespace HybZonoMPC;

void zonocost1(ZonoCpp::HybZono &X, std::map<int, Eigen::MatrixXd> &A_Hrep, std::map<int, Eigen::VectorXd> &b_Hrep,
    Eigen::VectorXd &pos_0, Eigen::VectorXd &pos_ref, double &T_sim, Eigen::VectorXd &cost_vec, Eigen::MatrixXd &d_mat,
    std::vector<Eigen::MatrixXd> &V_obs, std::vector<Eigen::MatrixXd> &V_free, Eigen::MatrixXd &V_lim)
{
    // range
    std::pair<double,double> x_range = std::make_pair(0, 25);
    std::pair<double,double> y_range = std::make_pair(-5, 5);    

    // initial / final positions
    pos_0.resize(2);
    pos_0 << 0, 2.5;
    pos_ref.resize(2);
    pos_ref << 25, 2.5;

    // sim time
    T_sim = 40;

    // positions of obstacles
    #define n_obs 4
    std::array<double, n_obs> xc_obs_rel = {0.25, 0.4, 0.65, 0.85};
    std::array<double, n_obs> yc_obs_rel = {0.2, 0.8, 0.8, 0.2};
    std::array<double, n_obs> xc_obs, yc_obs;
    for (int i=0; i<n_obs; i++)
    {
        xc_obs[i] = (x_range.second-x_range.first)*xc_obs_rel[i] + x_range.first;
        yc_obs[i] = (y_range.second-y_range.first)*yc_obs_rel[i] + y_range.first;
    }


    // occupied cells
    int nx = 25;
    int ny = 10;

    std::vector<std::pair<int,int>> occupied_cells;
    for (int i : {5,6,7})
    {
        for (int j : {1,2})
            occupied_cells.push_back(std::make_pair(i, j));
    }
    for (int i : {8,9,10,11})
    {
        for (int j : {7,8})
            occupied_cells.push_back(std::make_pair(i, j));
    }
    for (int i : {15,16,17})
    {
        for (int j : {7,8})
            occupied_cells.push_back(std::make_pair(i, j));
    }
    for (int i : {20,21,22})
    {
        for (int j : {1,2})
            occupied_cells.push_back(std::make_pair(i, j));
    }

    // free cells
    std::vector<std::pair<int, int>> free_space_cells;
    for (int i=0; i<nx; i++)
    {
        for (int j=0; j<ny; j++)
        {
            if (std::find(occupied_cells.begin(), occupied_cells.end(), std::make_pair(i, j)) == occupied_cells.end())
            {
                free_space_cells.push_back(std::make_pair(i, j));
            }
        }
    }

    // grid dimensions
    double dx = (x_range.second - x_range.first)/nx;
    double dy = (y_range.second - y_range.first)/ny;

    // get tile centers
    std::vector<double> xc_tiles, yc_tiles;
    for (auto it = free_space_cells.begin(); it != free_space_cells.end(); it++)
    {
        xc_tiles.push_back(x_range.first + (((double) it->first) + 0.5)*dx);
        yc_tiles.push_back(y_range.first + (((double) it->second) + 0.5)*dy);
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

    // define cost associated with grid
    double max_cost = 25;
    double sigma_gauss = 1.25*dx;
    cost_vec.resize(xc_tiles.size());
    for (int i=0; i<xc_tiles.size(); i++)
    {
        double x = xc_tiles[i];
        double y = yc_tiles[i];
        double cost = 0; // init
        for (int j=0; j<n_obs; j++)
        {
            double dist_from_obs = std::sqrt(std::pow(x - xc_obs[j], 2) + std::pow(y - yc_obs[j], 2));
            cost += max_cost*std::exp(-std::pow(dist_from_obs/sigma_gauss, 2) / (2*std::pow(sigma_gauss, 2)));
        }
        cost_vec(i) = cost>max_cost ? max_cost : cost;
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
    for (auto it=occupied_cells.begin(); it!=occupied_cells.end(); ++it)
    {
        Vc(0) = x_range.first + (((double) it->first)+0.5)*dx;
        Vc(1) = y_range.first + (((double) it->second)+0.5)*dy;
        Vtmp = Vc.replicate(4, 1) + dVr;
        V_obs.push_back(Vtmp);
    }

    // free space vertices
    V_free.clear();
    for (auto it=free_space_cells.begin(); it!=free_space_cells.end(); ++it)
    {
        Vc(0) = x_range.first + (((double) it->first)+0.5)*dx;
        Vc(1) = y_range.first + (((double) it->second)+0.5)*dy;
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