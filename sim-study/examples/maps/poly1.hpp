#ifndef _POLY1_HPP_
#define _POLY1_HPP_

#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "ZonoCpp.hpp"
#include "vrep_calcs.hpp"
#include <map>
#include <utility>
#include "polypartition.h"

using namespace HybZonoMPC;

void poly1(ZonoCpp::HybZono &X, std::map<int, Eigen::MatrixXd> &A_Hrep, std::map<int, Eigen::VectorXd> &b_Hrep,
    Eigen::VectorXd &pos_0, Eigen::VectorXd &pos_ref, double &T_sim,
    std::vector<Eigen::MatrixXd> &V_obs, std::vector<Eigen::MatrixXd> &V_free, Eigen::MatrixXd &V_lim)
{
    // range
    std::pair<double,double> x_range = std::make_pair(-20, 20);
    std::pair<double,double> y_range = std::make_pair(-10, 10);    

    // initial / final positions
    pos_0.resize(2);
    pos_0 << -20, 0;
    pos_ref.resize(2);
    pos_ref << 20, 0;

    // sim time
    T_sim = 50;

    // max partition area
    const double min_area = 1e-4;

    // adjustment to get inside map
    double eps = 1e-6;

    // obstacle vertices
    std::vector<std::pair<double,double>> obs1 = {
            {-8.4420,   -0.0208},
            {-10.8413,   -3.6000},
            {-5.7312,   -6.5159},
            {0.1647,    4.9766}};

    std::vector<std::pair<double,double>> obs2 = {
        {5.4836,   -2.8722},
        {14.6216,   -2.2904},
        {-0.1272,  -10.0000+eps},
        {7.3237,  -10.0000+eps},
        {-3.4280,   -7.2480},
        {-0.0700,   -3.4771}};

    std::vector<std::pair<double,double>> obs3 = {
            {14.8001,   -8.6078},
            {14.5745,   -5.5935},
            {12.1184,  -10.0000+eps},
            {12.9135,  -10.0000+eps}};

    std::vector<std::pair<double,double>> obs4 = {
            {14.7605,    2.5948},
            {14.1417,   -1.1628},
            {9.4639,   -2.4678}};

    std::vector<std::vector<std::pair<double,double>>> obs_vec = {obs1, obs2, obs3, obs4};
    for (auto it=obs_vec.begin(); it!=obs_vec.end(); ++it)
        order_pts(*it);

    // free space
    std::vector<std::pair<double,double>> free_space = {
        {x_range.first, y_range.first},
        {x_range.second, y_range.first},
        {x_range.second, y_range.second},
        {x_range.first, y_range.second}};
    order_pts(free_space);

    // partition the free space using polypartition
    std::vector<std::vector<std::pair<double,double>>> poly_part_vec = part_HM(free_space, obs_vec);

    // remove polytopes below min area
    remove_polys_below_area(poly_part_vec, min_area);

    // convert to H-rep
    A_Hrep.clear(); // reset
    b_Hrep.clear(); // reset
    std::pair<Eigen::MatrixXd, Eigen::VectorXd> poly_H;
    int region = 1; // 1-based indexing
    for (auto it=poly_part_vec.begin(); it!=poly_part_vec.end(); ++it)
    {
        vrep2hrep_2D(*it, poly_H);
        A_Hrep[region] = poly_H.first;
        b_Hrep[region] = poly_H.second;
        region++;
    }

    // get hybzono

    // vertex and incidence matrices
    Eigen::MatrixXd V, M;
    poly2mat(poly_part_vec, V, M);

    // construct hybzono
    X.set(V, M);

    // obstacle vertices
    V_obs.clear();
    for (auto it=obs_vec.begin(); it!=obs_vec.end(); ++it)
    {
        Eigen::MatrixXd V_obs_temp(it->size(), 2);
        for (int i=0; i<it->size(); i++)
        {
            V_obs_temp(i, 0) = it->at(i).first;
            V_obs_temp(i, 1) = it->at(i).second;
        }
        V_obs.push_back(V_obs_temp);
    }

    // free space vertices
    V_free.clear();
    for (auto it=poly_part_vec.begin(); it!=poly_part_vec.end(); ++it)
    {
        Eigen::MatrixXd V_free_temp(it->size(), 2);
        for (int i=0; i<it->size(); i++)
        {
            V_free_temp(i, 0) = it->at(i).first;
            V_free_temp(i, 1) = it->at(i).second;
        }
        V_free.push_back(V_free_temp);
    }

    // limit vertices
    V_lim.resize(4, 2);
    V_lim << x_range.first, y_range.first,
             x_range.second, y_range.first,
             x_range.second, y_range.second,
             x_range.first, y_range.second;
}

#endif