#ifndef _VREP_CALCS_HPP_
#define _VREP_CALCS_HPP_

#include <vector>
#include <utility>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <string>
#include "polypartition.h"
#include "Eigen/Dense"

struct PtAngle
{
    // fields
    std::pair<double,double> pt;
    double angle;

    // constructor
    PtAngle(std::pair<double,double> pt, double angle) : pt(pt), angle(angle) {}
    
    // comparison operators
    bool operator<(const PtAngle &other) const
    {
        return angle < other.angle;
    }

    bool operator>(const PtAngle &other) const
    {
        return angle > other.angle;
    }
};

void order_pts(std::vector<std::pair<double,double>> &pts, bool cw = true)
{
    // find the centroid
    double sum_x = 0;
    double sum_y = 0;
    for (auto it=pts.begin(); it!=pts.end(); ++it)
    {
        sum_x += it->first;
        sum_y += it->second;
    }
    double cx = sum_x / pts.size();
    double cy = sum_y / pts.size();

    // create a vector of points and their angles
    std::vector<PtAngle> pt_angles;
    for (auto it=pts.begin(); it!=pts.end(); ++it)
    {
        double dx = it->first - cx;
        double dy = it->second - cy;
        double angle = std::atan2(dy, dx);
        pt_angles.push_back(PtAngle(*it, angle));
    }

    // sort the points by angle
    if (cw)
        std::sort(pt_angles.begin(), pt_angles.end(), std::greater<PtAngle>());
    else
        std::sort(pt_angles.begin(), pt_angles.end());

    // put into pts variable
    pts.clear();
    for (auto it=pt_angles.begin(); it!=pt_angles.end(); ++it)
        pts.push_back(it->pt);
}

// make polygon
void make_polygon(const std::vector<std::pair<double, double>> &pairs, 
    bool is_hole, TPPLPoly &poly)
{
    // number of vertices
    int num_vertices = pairs.size();
    poly.Init(num_vertices);

    // add vertices to polygon
    for (int i=0; i<num_vertices; i++)
    {
        poly[i].x = pairs[i].first;
        poly[i].y = pairs[i].second;
    }

    // handle hole logic
    poly.SetHole(is_hole);
    if (is_hole)
        poly.SetOrientation(TPPLOrientation::TPPL_ORIENTATION_CW);
    else
        poly.SetOrientation(TPPLOrientation::TPPL_ORIENTATION_CCW);
}

// partition
std::vector<std::vector<std::pair<double,double>>> part_HM(const std::vector<std::pair<double,double>> &poly_free,
    const std::vector<std::vector<std::pair<double,double>>> &poly_obs_vec)
{
    // declare
    TPPLPolyList polys_w_hole;
    TPPLPoly poly;

    // outer polygon
    make_polygon(poly_free, false, poly);
    polys_w_hole.push_back(poly);

    // holes
    for (auto it=poly_obs_vec.begin(); it!=poly_obs_vec.end(); ++it)
    {
        make_polygon(*it, true, poly);
        polys_w_hole.push_back(poly);
    }

    // declare partition object
    TPPLPartition partition;

    // remove holes
    TPPLPolyList polys_wo_hole;
    int out_holes = partition.RemoveHoles(&polys_w_hole, &polys_wo_hole);
    if (out_holes == 0)
        throw std::runtime_error("Hole removal failed");

    TPPLPoly poly_nohole = *(polys_wo_hole.begin());

    // partition
    TPPLPolyList polys_partitioned;
    int out_part = partition.ConvexPartition_HM(&poly_nohole, &polys_partitioned);
    if (out_part == 0)
        throw std::runtime_error("Partitioning failed");

    // put into vector of vector of pairs
    std::vector<std::vector<std::pair<double,double>>> poly_part_vec;
    for (auto it=polys_partitioned.begin(); it!=polys_partitioned.end(); ++it)
    {
        std::vector<std::pair<double,double>> poly_part;
        for (int i=0; i<it->GetNumPoints(); i++)
            poly_part.push_back(std::make_pair(it->GetPoint(i).x, it->GetPoint(i).y));
        poly_part_vec.push_back(poly_part);
    }

    return poly_part_vec;
}

// sign operation
double sgn_dbl(double val) 
{
    if (val < 0)
        return -1;
    else if (val > 0)
        return 1;
    else
        return 0;
}

// V-rep to H-rep (only for 2D)
void vrep2hrep_2D(std::vector<std::pair<double,double>> poly_V, 
    std::pair<Eigen::MatrixXd, Eigen::VectorXd> &poly_H)
{
    // sort the vertices to be counter-clockwise
    order_pts(poly_V, false);

    // number of vertices
    int num_vertices = poly_V.size();

    // create A and b matrices
    Eigen::MatrixXd A(num_vertices, 2);
    Eigen::VectorXd b(num_vertices);

    // declare
    double xi, yi, xip1, yip1, dx, dy, b_raw;

    // fill in A and b
    for (int i=0; i<num_vertices; i++)
    {
        xi = poly_V[i].first;
        yi = poly_V[i].second;
        xip1 = poly_V[(i+1) % num_vertices].first;
        yip1 = poly_V[(i+1) % num_vertices].second;
        dx = xip1 - xi;
        dy = yip1 - yi;

        b_raw = dy*xi - dx*yi;  
        if (std::abs(b_raw) > 1e-10) // normalize
        {
            A(i, 0) = dy/std::abs(b_raw);
            A(i, 1) = -dx/std::abs(b_raw);
            b(i) = sgn_dbl(b_raw) * 1.0;
        }
        else // don't normalize
        {
            A(i, 0) = dy;
            A(i, 1) = -dx;
            b(i) = b_raw;
        }
    }

    // put into pair
    poly_H = std::make_pair(A, b);
}

// compute vertex and incidence matrices for a vector polygons
void poly2mat(const std::vector<std::vector<std::pair<double,double>>> &Vpoly_vec,
    Eigen::MatrixXd &V, Eigen::MatrixXd &M)
{
    // create vector of all vertices
    std::vector<std::pair<double,double>> V_all;
    for (auto it_poly = Vpoly_vec.begin(); it_poly != Vpoly_vec.end(); ++it_poly)
    {
        for (auto it_pt = it_poly->begin(); it_pt != it_poly->end(); ++it_pt)
            V_all.push_back(*it_pt);
    }

    // sort and remove duplicates
    std::sort(V_all.begin(), V_all.end());
    auto it_end = std::unique(V_all.begin(), V_all.end());
    V_all.resize(std::distance(V_all.begin(), it_end));

    // put into matrix
    Eigen::MatrixXd VT (V_all.size(), 2);
    for (int i=0; i<V_all.size(); i++)
    {
        VT(i, 0) = V_all[i].first;
        VT(i, 1) = V_all[i].second;
    }
    V = VT.transpose();

    // get incidence matrix
    int n_verts = V_all.size();
    int n_partitions = Vpoly_vec.size();
    M = Eigen::MatrixXd::Zero(n_verts, n_partitions);
    for (int i=0; i<n_verts; i++)
    {
        for (int j=0; j<n_partitions; j++)
        {
            if (std::find(Vpoly_vec[j].begin(), Vpoly_vec[j].end(), V_all[i]) != Vpoly_vec[j].end())
                M(i, j) = 1;
        }
    }
}

// compute area
double poly_area(const std::vector<std::pair<double,double>> &poly)
{
    // number of vertices
    int num_vertices = poly.size();

    // declare
    double xi, yi, xip1, yip1, area = 0;

    // compute area
    for (int i=0; i<num_vertices; i++)
    {
        xi = poly[i].first;
        yi = poly[i].second;
        xip1 = poly[(i+1) % num_vertices].first;
        yip1 = poly[(i+1) % num_vertices].second;
        area += xi*yip1 - xip1*yi;
    }

    return 0.5 * std::abs(area);
}

// remove polytopes below given area
void remove_polys_below_area(std::vector<std::vector<std::pair<double,double>>> &Vpoly_vec, double min_area)
{
    // iterate through vector of polygons
    for (auto it=Vpoly_vec.begin(); it!=Vpoly_vec.end();)
    {
        // compute area
        double area = poly_area(*it);

        // remove if below max area
        if (area < min_area)
            it = Vpoly_vec.erase(it);
        else
            ++it;
    }
}

// vrep to file
void vrep2file(const std::vector<Eigen::MatrixXd> &Vobs, const std::vector<Eigen::MatrixXd> &Vfree,
    const Eigen::MatrixXd &V_lim, const std::string &map_name)
{
    // create filename
    std::string filename = "../examples/data/" + map_name + "_vrep.txt";

    // open file
    std::ofstream file;

    // write to file
    file.open(filename);
    file << "Limits" << std::endl;
    for (int i=0; i<V_lim.rows(); i++)
        file << V_lim.row(i) << std::endl;
    file << "---" << std::endl;

    file << "Obstacles" << std::endl;
    for (auto it=Vobs.begin(); it!=Vobs.end(); ++it)
    {
        for (int i=0; i<it->rows(); i++)
            file << it->row(i) << std::endl;
        file << "---" << std::endl;
    }

    file << "Free Space" << std::endl;
    for (auto it=Vfree.begin(); it!=Vfree.end(); ++it)
    {
        for (int i=0; i<it->rows(); i++)
            file << it->row(i) << std::endl;
        file << "---" << std::endl;
    }

    // close file
    file.close();
}


void vrep2file(const std::vector<Eigen::MatrixXd> &Vobs, const std::vector<Eigen::MatrixXd> &Vfree,
    const Eigen::MatrixXd &V_lim, const Eigen::VectorXd &cost_vec, const std::string &map_name)
{
    // create filename
    std::string filename = "../examples/data/" + map_name + "_vrep.txt";

    // open file
    std::ofstream file;

    // write to file
    file.open(filename);
    file << "Limits" << std::endl;
    for (int i=0; i<V_lim.rows(); i++)
        file << V_lim.row(i) << std::endl;
    file << "---" << std::endl;

    file << "Obstacles" << std::endl;
    for (auto it=Vobs.begin(); it!=Vobs.end(); ++it)
    {
        for (int i=0; i<it->rows(); i++)
            file << it->row(i) << std::endl;
        file << "---" << std::endl;
    }

    file << "Free Space" << std::endl;
    for (auto it=Vfree.begin(); it!=Vfree.end(); ++it)
    {
        for (int i=0; i<it->rows(); i++)
            file << it->row(i) << " " << cost_vec(i) << std::endl;
        file << "---" << std::endl;
    }

    file << "Costs" << std::endl;
    for (int i=0; i<cost_vec.size(); i++)
        file << cost_vec(i) << std::endl;

    // close file
    file.close();
}


#endif