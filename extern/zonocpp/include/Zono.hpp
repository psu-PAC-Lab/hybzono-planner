#ifndef __ZONO_HPP__
#define __ZONO_HPP__

#include "Eigen/Dense"
#include "Eigen/Sparse"
#include <sstream>

namespace ZonoCpp
{

class Zono
{
    public:

    // fields
    Eigen::SparseMatrix<double> G;
    Eigen::VectorXd c;
    int n, nG;
    bool zero_one_generators = false;

    // constructor
    Zono()
    {
        nG = 0;
        n = 0;
        G = Eigen::SparseMatrix<double>(0, 0);
        c = Eigen::VectorXd(0);
    }

    Zono(const Eigen::SparseMatrix<double> &G, const Eigen::VectorXd &c,
        bool zero_one_generators=false)
    {
        set(G, c, zero_one_generators);
    }

    // set method
    void set(const Eigen::SparseMatrix<double> &G, const Eigen::VectorXd &c,
        bool zero_one_generators=false)
    {
        this->G = G;
        this->c = c;
        nG = G.cols();
        n = G.rows();
        this->zero_one_generators = zero_one_generators;
    }

    // generator conversion between [-1,1] and [0,1]
    void convert_generator_range()
    {
        Eigen::VectorXd c;
        Eigen::SparseMatrix<double> G;

        if (!zero_one_generators) // convert to [0,1] generators
        {
            c = this->c - this->G*Eigen::VectorXd::Ones(this->nG);
            G = 2*this->G;

            set(G, c, true);
        }
        else // convert to [-1,1] generators
        {
            c = this->c + 0.5*this->G*Eigen::VectorXd::Ones(this->nG);
            G = 0.5*this->G;

            set(G, c, false);
        }
    }

    // display methods
    friend std::ostream& operator<<(std::ostream &os, const Zono &Z)
    {
        os << "Zono: " << std::endl;
        os << "n: " << Z.n << std::endl;
        os << "nG: " << Z.nG << std::endl;
        os << "G: " << Eigen::MatrixXd(Z.G) << std::endl;
        os << "c: " << Z.c << std::endl;
        os << "zero_one_generators: " << Z.zero_one_generators << std::endl;
        return os;
    }

    std::string print()
    {
        std::stringstream ss;
        ss << "Zono: " << std::endl;
        ss << "n: " << this->n << std::endl;
        ss << "nG: " << this->nG << std::endl;
        ss << "G: " << Eigen::MatrixXd(this->G) << std::endl;
        ss << "c: " << this->c << std::endl;
        ss << "zero_one_generators: " << this->zero_one_generators << std::endl;
        return ss.str();
    }
};

} // namespace ZonoCpp

#endif