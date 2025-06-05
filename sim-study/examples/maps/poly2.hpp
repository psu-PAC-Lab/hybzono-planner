#ifndef _POLY2_HPP_
#define _POLY2_HPP_

#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "ZonoCpp.hpp"
#include "vrep_calcs.hpp"
#include <map>
#include <utility>
#include "polypartition.h"

using namespace HybZonoMPC;

void poly2(ZonoCpp::HybZono &X, std::map<int, Eigen::MatrixXd> &A_Hrep, std::map<int, Eigen::VectorXd> &b_Hrep,
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

    // obstacle vertices
    std::vector<std::pair<double,double>> obs1 = {
            {-3.984556013353694, -1.301164139544404},
            {-5.523791724984051, -6.800061947461296},
            {4.658729278331522, -9.999999000000017},
            {-2.345384754522588, -9.999999000000017},
            {8.620220475239497, -6.547317200323886},
            {0.830500242792606, 4.387700237243461},
            {-2.261328356812852, 1.667311522938448}};

    std::vector<std::pair<double,double>> obs2 = {
        {11.596365301924079, 5.674900893365844},
        {4.743078372628562, 1.647444129829116},
        {9.713509794311136, 9.999999000000006},
        {3.846801197079755, 9.999999000000006},
        {1.744052709329245, 5.672310632448288}};

    std::vector<std::pair<double,double>> obs3 = {
            {-12.128910351505571, -0.468356236963613},
            {-10.378286279099317, 4.978073629913993},
            {-11.503284893166358, 4.067017294893390},
            {-8.039259675853062, 5.864085420571117},
            {-2.463259564905091, 4.892243072726099}};

    std::vector<std::vector<std::pair<double,double>>> obs_vec = {obs1, obs2, obs3};
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