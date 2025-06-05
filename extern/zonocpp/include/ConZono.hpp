#ifndef _CONZONO_HPP_
#define _CONZONO_HPP_

#include "Eigen/Dense"
#include "Eigen/Sparse"
#include <sstream>

namespace ZonoCpp
{

class ConZono
{
    public:

    // fields
    Eigen::SparseMatrix<double> G, A;
    Eigen::VectorXd c, b;
    int n, nC, nG;
    bool zero_one_generators = false;

    // constructor
    ConZono()
    {
        nG = 0;
        nC = 0;
        n = 0;
        G = Eigen::SparseMatrix<double>(0, 0);
        A = Eigen::SparseMatrix<double>(0, 0);
        c = Eigen::VectorXd(0);
        b = Eigen::VectorXd(0);
    }

    ConZono(const Eigen::SparseMatrix<double> &G, const Eigen::VectorXd &c,
        const Eigen::SparseMatrix<double> &A, const Eigen::VectorXd &b, 
        bool zero_one_generators=false)
    {
        set(G, c, A, b, zero_one_generators);
    }

    // set method
    void set(const Eigen::SparseMatrix<double> &G, const Eigen::VectorXd &c,
        const Eigen::SparseMatrix<double> &A, const Eigen::VectorXd &b,
        bool zero_one_generators=false)
    {
        this->G = G;
        this->A = A;
        this->c = c;
        this->b = b;
        nG = G.cols();
        nC = A.rows();
        n = G.rows();
        this->zero_one_generators = zero_one_generators;
    }

    // generator conversion between [-1,1] and [0,1]
    void convert_generator_range()
    {
        Eigen::VectorXd c, b;
        Eigen::SparseMatrix<double> G, A;

        if (!zero_one_generators) // convert to [0,1] generators
        {
            c = this->c - this->G*Eigen::VectorXd::Ones(this->nG);
            b = this->b + this->A*Eigen::VectorXd::Ones(this->nG);
            G = 2*this->G;
            A = 2*this->A;

            set(G, c, A, b, true);
        }
        else // convert to [-1,1] generators
        {
            c = this->c + 0.5*this->G*Eigen::VectorXd::Ones(this->nG);
            b = this->b - 0.5*this->A*Eigen::VectorXd::Ones(this->nG);
            G = 0.5*this->G;
            A = 0.5*this->A;

            set(G, c, A, b, false);
        }
    }

    // display methods
    friend std::ostream& operator<<(std::ostream &os, const ConZono &Zc)
    {
        os << "ConZono: " << std::endl;
        os << "n: " << Zc.n << std::endl;
        os << "nG: " << Zc.nG << std::endl;
        os << "nC: " << Zc.nC << std::endl;
        os << "G: " << Eigen::MatrixXd(Zc.G) << std::endl;
        os << "c: " << Zc.c << std::endl;
        os << "A: " << Eigen::MatrixXd(Zc.A) << std::endl;
        os << "b: " << Zc.b << std::endl;
        os << "zero_one_generators: " << Zc.zero_one_generators << std::endl;
        return os;
    }

    std::string print()
    {
        std::stringstream ss;
        ss << "ConZono: " << std::endl;
        ss << "n: " << this->n << std::endl;
        ss << "nG: " << this->nG << std::endl;
        ss << "nC: " << this->nC << std::endl;
        ss << "G: " << Eigen::MatrixXd(this->G) << std::endl;
        ss << "c: " << this->c << std::endl;
        ss << "A: " << Eigen::MatrixXd(this->A) << std::endl;
        ss << "b: " << this->b << std::endl;
        ss << "zero_one_generators: " << this->zero_one_generators << std::endl;
        return ss.str();
    }
};

} // namespace ZonoCpp


#endif