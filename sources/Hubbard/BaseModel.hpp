#pragma once
#include <Eigen/Dense>
#include <iostream>
#include "ModelParameters.hpp"

namespace Hubbard {
    // Defines the working precision of the entire project
	// Change to float, double or long double - so far double produces the best results
	typedef double double_prec;
	constexpr double_prec L_PI = 3.141592653589793238462643383279502884L; /* pi */
	typedef Eigen::Matrix<double_prec, Eigen::Dynamic, Eigen::Dynamic> Matrix_L;
	typedef Eigen::Vector<double_prec, Eigen::Dynamic> Vector_L;

	typedef std::complex<double_prec> complex_prec;
	typedef Eigen::Matrix<complex_prec, Eigen::Dynamic, Eigen::Dynamic> MatrixCL;
	typedef Eigen::Vector<complex_prec, Eigen::Dynamic> VectorCL;

	template <const int vector_size>
	void printAsRow(Eigen::Vector<double_prec, vector_size>& printer) {
		for (size_t i = 0; i < printer.size(); i++)
		{
			std::cout << "\t" << printer(i);
			if ((i + 1) % 8 == 0) {
				std::cout << "\n\t    ";
			}
		}
		std::cout << std::endl;
	}

	template <const int vector_size>
	void printAsRow(Eigen::Vector<complex_prec, vector_size>& printer) {
		for (size_t i = 0; i < printer.size(); i++)
		{
			std::cout << " \t" << printer(i);
			if ((i + 1) % 4 == 0) {
				std::cout << "\n\t    ";
			}
		}
		std::cout << std::endl;
	}

    class BaseModel {
    protected:
        typedef Eigen::Matrix<complex_prec, Eigen::Dynamic, Eigen::Dynamic> SpinorMatrix;
        const complex_prec I = { 0, 1 };
		SpinorMatrix hamilton;

        double_prec temperature;
		double_prec U;
		double_prec V;
		double_prec U_OVER_N;
		double_prec V_OVER_N;

        size_t TOTAL_BASIS;
		double_prec delta_sc, delta_cdw, delta_afm, delta_eta;
        double_prec chemical_potential;

		inline virtual void computeChemicalPotential() {
			this->chemical_potential = 0.5 * U + 4 * V;
		};

        template<typename... Args>
		inline double_prec gamma(Args... ks) const {
			return (cos(ks) + ...);
		}
        template<typename... Args>
		inline double_prec unperturbed_energy(Args... ks) const {
			return -2 * gamma(ks...);
		};

        inline double_prec fermi_dirac(double_prec energy) const {
			if (temperature > 1e-8) {
				return (1. / (1 + exp(energy / temperature)));
			}
			else {
				if (std::abs(energy) < 1e-12) {
					return 0.5;
				}
				return ((energy > 0) ? 0 : 1);
			}
		};
		inline void fillRho(SpinorMatrix& rho, const Eigen::SelfAdjointEigenSolver<SpinorMatrix>& solvedHamilton) const {
			rho.fill(0);
			for (int i = 0; i < rho.rows(); i++)
			{
				rho(i, i) = 1 - fermi_dirac(solvedHamilton.eigenvalues()(i));
			}
			rho = solvedHamilton.eigenvectors() * rho * solvedHamilton.eigenvectors().adjoint();
		};

	public:
		BaseModel(const ModelParameters& _params);
    };
}