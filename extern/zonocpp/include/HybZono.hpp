#ifndef _HYBZONO_HPP_
#define _HYBZONO_HPP_

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <sstream>
#include <vector>
#include "ConZono.hpp"

namespace ZonoCpp
{

class HybZono
{
    public:

        // fields
        Eigen::SparseMatrix<double> Gc, Gb, Ac, Ab, G, A;
        Eigen::VectorXd c, b;
        int n, nGc, nGb, nC, nG;
        bool zero_one_generators = false;
        
        // constructors
        HybZono() = default;
        HybZono(const Eigen::SparseMatrix<double> &Gc, const Eigen::SparseMatrix<double> &Gb, const Eigen::VectorXd &c,
            const Eigen::SparseMatrix<double> &Ac, const Eigen::SparseMatrix<double> &Ab,  const Eigen::VectorXd &b,
            bool zero_one_generators=false)
        {
            set(Gc, Gb, c, Ac, Ab, b, zero_one_generators);
        } 
        HybZono(const Eigen::Ref<const Eigen::MatrixXd> V, const Eigen::Ref<const Eigen::MatrixXd> M)
        {
            set(V, M);
        }

        // set methods
        void set(const Eigen::SparseMatrix<double> &Gc, const Eigen::SparseMatrix<double> &Gb, const Eigen::VectorXd &c,
            const Eigen::SparseMatrix<double> &Ac, const Eigen::SparseMatrix<double> &Ab,  const Eigen::VectorXd &b,
            bool zero_one_generators=false)
        {
            this->Gc = Gc;
            this->Gb = Gb;
            this->Ac = Ac;
            this->Ab = Ab;
            this->c = c;
            this->b = b;
            nGc = Gc.cols();
            nGb = Gb.cols();
            nC = Ac.rows();
            n = Gc.rows();
            this->zero_one_generators = zero_one_generators;

            make_G_A();
        }

        // set using V-rep to hybzono
        void set(const Eigen::Ref<const Eigen::MatrixXd> V, const Eigen::Ref<const Eigen::MatrixXd> M)
        {
            // dimensions
            int n = V.rows();
            int nV = V.cols();
            int N = M.cols();

            // Q hybzono
            Eigen::MatrixXd Gc (nV+N, nV);
            Gc << 0.5*Eigen::MatrixXd::Identity(nV, nV), Eigen::MatrixXd::Zero(N, nV);

            Eigen::MatrixXd Gb (nV+N, N);
            Gb << Eigen::MatrixXd::Zero(nV, N), 0.5*Eigen::MatrixXd::Identity(N, N);

            Eigen::VectorXd c = 0.5*Eigen::VectorXd::Ones(nV+N);

            Eigen::MatrixXd Ac (2, nV);
            Ac << Eigen::MatrixXd::Ones(1, nV), Eigen::MatrixXd::Zero(1, nV);

            Eigen::MatrixXd Ab (2, N);
            Ab << Eigen::MatrixXd::Zero(1, N), Eigen::MatrixXd::Ones(1, N);

            Eigen::VectorXd b (2);
            b << 2-nV, 2-N;

            HybZono Q (Gc.sparseView(), Gb.sparseView(), c, Ac.sparseView(), Ab.sparseView(), b);

            // halfspace intersection
            Eigen::MatrixXd R (nV, nV + N);
            R << Eigen::MatrixXd::Identity(nV, nV), -1*M;
            HybZono D = Q.halfspace_intersection(Eigen::MatrixXd::Identity(nV, nV), Eigen::MatrixXd::Zero(nV, 1), R);

            // premultiply
            Eigen::MatrixXd P (n, nV + N);
            P << V, Eigen::MatrixXd::Zero(n, N);
            *this = D.premultiply(P);

            // make G, A
            make_G_A();
        }

        // set using zono union
        void set_from_zono_union(const Eigen::Ref<const Eigen::MatrixXd> S, const Eigen::Ref<const Eigen::MatrixXd> M, 
            const Eigen::Ref<const Eigen::MatrixXd> C)
        {
            // n is dimension of space
            // n_s is number of line segments
            // n_r is number of zonotopes (i.e., convex regions)
            // S is (n x n_s) matrix of line segments: S = [s0, s1, ..., sN-1]
            // M is (n_s x n_r) incidence matrix indicating which line segments appear in each zonotope
            // C is (n x n_r) matrix of zonotope offsets    

            // dimensions
            this->n = S.rows();
            int n_s = S.cols();
            int n_r = C.cols();
            this->nGc = 2*n_s;
            this->nGb = n_r;
            this->nC = n_s + 1;

            // generator matrices
            this->Gc = S.sparseView();
            this->Gc.conservativeResize(this->n, this->nGc); // fill the rest with zeros

            this->Gb = C.sparseView();

            this->c = Eigen::VectorXd::Zero(this->n);

            // constraint matrices
            std::vector<Eigen::Triplet<double>> triplets;
            triplets.reserve(2*n_s);
            for (int i=0; i<n_s; i++)
            {
                // count the number of appearances of the i-th line segment
                int num_seg_appearances = M.row(i).sum();
                triplets.push_back(Eigen::Triplet<double>(i, i, 1));
                triplets.push_back(Eigen::Triplet<double>(i, n_s+i, num_seg_appearances));
            }
            this->Ac.resize(this->nC, this->nGc);
            this->Ac.setFromTriplets(triplets.begin(), triplets.end());

            this->Ab = -1*M.sparseView();
            this->Ab.conservativeResize(this->nC, this->nGb); // fill the extra row with zeros
            for (int i=0; i<this->nGb; i++) // now fill that row with ones
            {
                this->Ab.insert(this->nC-1, i) = 1;
            }

            this->b = Eigen::VectorXd::Zero(this->nC);
            this->b(this->nC-1) = 1;

            // this is using zero one generators
            this->zero_one_generators = true;

            // make G, A
            make_G_A();
        }


        // generator conversion between [-1,1] and [0,1]
        void convert_generator_range()
        {
            Eigen::VectorXd c, b;
            Eigen::SparseMatrix<double> Gb, Ab, Ac, Gc;

            if (!zero_one_generators) // convert to [0,1] generators
            {
                c = this->c - this->G*Eigen::VectorXd::Ones(this->nG);
                b = this->b + this->A*Eigen::VectorXd::Ones(this->nG);
                Gb = 2*this->Gb;
                Ab = 2*this->Ab;
                Gc = 2*this->Gc;
                Ac = 2*this->Ac;

                set(Gc, Gb, c, Ac, Ab, b, true);
            }
            else // convert to [-1,1] generators
            {
                c = this->c + 0.5*this->G*Eigen::VectorXd::Ones(this->nG);
                b = this->b - 0.5*this->A*Eigen::VectorXd::Ones(this->nG);
                Gb = 0.5*this->Gb;
                Ab = 0.5*this->Ab;
                Gc = 0.5*this->Gc;
                Ac = 0.5*this->Ac;

                set(Gc, Gb, c, Ac, Ab, b, false);
            }
        }

        // convex relaxation
        ConZono convex_relaxation()
        {
            return ConZono(this->G, this->c, this->A, this->b, this->zero_one_generators);
        }

        // halfspace intersection
        HybZono halfspace_intersection(const Eigen::Ref<const Eigen::MatrixXd> H, 
            const Eigen::Ref<const Eigen::VectorXd> f, const Eigen::Ref<const Eigen::MatrixXd> R)
        {
            // H - nH x m real matrix defining the nH halfspace normal vectors
            // f - nH x 1 real vector defining the nH offsets
            // R - m x n real matrix

            // dim
            int nH = H.rows();

            // convert to {-1,1} binaries if necessary
            if (zero_one_generators)
                convert_generator_range();

            // init working variables
            Eigen::MatrixXd Ac = Eigen::MatrixXd(this->Ac);
            Eigen::MatrixXd Ab = Eigen::MatrixXd(this->Ab);
            Eigen::VectorXd b = this->b;
            Eigen::MatrixXd Gc = Eigen::MatrixXd(this->Gc);
            Eigen::MatrixXd Gb = Eigen::MatrixXd(this->Gb);
            Eigen::VectorXd c = this->c;
            double d_max;
            Eigen::MatrixXd Gc_Gb;
            Eigen::MatrixXd Ac_kp1, Ab_kp1, b_kp1, Gc_kp1;

            // compute intersection
            for (int i=0; i<nH; i++)
            {
                Gc_Gb = Eigen::MatrixXd(Gc.rows(), Gc.cols() + Gb.cols());
                Gc_Gb << Gc, Gb;

                // Computes the distance of the hyperplane from the farthest vertex of the zonotope
                d_max = f(i) - (H.row(i)*R)*c + ((H.row(i)*R*Gc_Gb).cwiseAbs()).sum();

                Ac_kp1 = Eigen::MatrixXd(Ac.rows() + 1, Ac.cols() + 1);
                Ac_kp1 << Ac, Eigen::MatrixXd::Zero(Ac.rows(), 1), H.row(i)*R*Gc, d_max/2;

                Ab_kp1 = Eigen::MatrixXd(Ab.rows() + 1, Ab.cols());
                Ab_kp1 << Ab, H.row(i)*R*Gb;

                b_kp1 = Eigen::VectorXd(b.rows() + 1);
                b_kp1 << b, f(i) - H.row(i)*R*c - d_max/2;

                Gc_kp1 = Eigen::MatrixXd(Gc.rows(), Gc.cols() + 1);
                Gc_kp1 << Gc, Eigen::MatrixXd::Zero(Gc.rows(), 1);

                // increment
                Ac = Ac_kp1;
                Ab = Ab_kp1;
                b = b_kp1;
                Gc = Gc_kp1;
            }

            return HybZono(Gc.sparseView(), Gb.sparseView(), c, 
                Ac.sparseView(), Ab.sparseView(), b, false);
        }

        // premultiply by matrix
        HybZono premultiply(const Eigen::Ref<const Eigen::MatrixXd> M)
        {
            Eigen::SparseMatrix<double> M_sp = M.sparseView();
            return premultiply(M_sp);
        }
        
        HybZono premultiply(Eigen::SparseMatrix<double> &M)
        {
            Eigen::SparseMatrix<double> Gc, Gb;
            Eigen::VectorXd c;

            // premultiply
            Gc = M*this->Gc;
            Gb = M*this->Gb;
            c = M*this->c;

            return HybZono(Gc, Gb, c, this->Ac, this->Ab, this->b, this->zero_one_generators);
        }

        // display methods
        friend std::ostream& operator<<(std::ostream &os, const HybZono &Zh)
        {
            os << "HybZono: " << std::endl;
            os << "n: " << Zh.n << std::endl;
            os << "nGc: " << Zh.nGc << std::endl;
            os << "nGb: " << Zh.nGb << std::endl;
            os << "nC: " << Zh.nC << std::endl;
            os << "Gc: " << Eigen::MatrixXd(Zh.Gc) << std::endl;
            os << "Gb: " << Eigen::MatrixXd(Zh.Gb) << std::endl;
            os << "c: " << Zh.c << std::endl;
            os << "Ac: " << Eigen::MatrixXd(Zh.Ac) << std::endl;
            os << "Ab: " << Eigen::MatrixXd(Zh.Ab) << std::endl;
            os << "b: " << Zh.b << std::endl;
            os << "zero_one_generators: " << Zh.zero_one_generators;
            return os;
        }

        std::string print()
        {
            std::stringstream ss;
            ss << "HybZono: " << std::endl;
            ss << "n: " << this->n << std::endl;
            ss << "nGc: " << this->nGc << std::endl;
            ss << "nGb: " << this->nGb << std::endl;
            ss << "nC: " << this->nC << std::endl;
            ss << "Gc: " << Eigen::MatrixXd(this->Gc) << std::endl;
            ss << "Gb: " << Eigen::MatrixXd(this->Gb) << std::endl;
            ss << "c: " << this->c << std::endl;
            ss << "Ac: " << Eigen::MatrixXd(this->Ac) << std::endl;
            ss << "Ab: " << Eigen::MatrixXd(this->Ab) << std::endl;
            ss << "b: " << this->b << std::endl;
            ss << "zero_one_generators: " << this->zero_one_generators;
            return ss.str();
        }

    protected:

        // make G, A
        void make_G_A()
        {
            // implemenent get_triplets_offset as a lambda function
            auto get_triplets_offset = [] (const Eigen::SparseMatrix<double> &mat, std::vector<Eigen::Triplet<double>> &triplets, 
                    int i_offset, int j_offset) -> void
            {
                for (int k=0; k<mat.outerSize(); ++k)
                {
                    for (Eigen::SparseMatrix<double>::InnerIterator it(mat, k); it; ++it)
                    {
                        triplets.push_back(Eigen::Triplet<double>(it.row() + i_offset, it.col() + j_offset, it.value()));
                    }
                }
            };

            std::vector<Eigen::Triplet<double>> G_triplets;
            G_triplets.reserve(Gc.nonZeros() + Gb.nonZeros());
            get_triplets_offset(Gc, G_triplets, 0, 0);
            get_triplets_offset(Gb, G_triplets, 0, nGc);
            this->G.resize(n, nGc + nGb);
            this->G.setFromTriplets(G_triplets.begin(), G_triplets.end());

            std::vector<Eigen::Triplet<double>> A_triplets;
            A_triplets.reserve(Ac.nonZeros() + Ab.nonZeros());
            get_triplets_offset(Ac, A_triplets, 0, 0);
            get_triplets_offset(Ab, A_triplets, 0, nGc);
            this->A.resize(nC, nGc + nGb);
            this->A.setFromTriplets(A_triplets.begin(), A_triplets.end());

            this->nG = nGc + nGb;
        }
};

} // end namespace ZonoCpp

#endif