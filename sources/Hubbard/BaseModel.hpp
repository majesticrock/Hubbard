#pragma once
#include <Eigen/Dense>
#include <iostream>
#include <cmath>
#include <cstdarg>
#include <type_traits>
#include "Constants.hpp"
#include "ModelParameters.hpp"
#include "../Utility/Resolvent.hpp"
#include "BaseModelAttributes.hpp"

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
	using SpinorMatrix = MatrixCL;
	typedef Utility::Resolvent<double_prec> Resolvent_L;

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

	typedef Eigen::Vector<complex_prec, Eigen::Dynamic> ComplexParameterVector;

	template <typename DataType>
	class BaseModel :
		public std::conditional_t<std::is_same_v<DataType, double_prec>, BaseModelRealAttributes, BaseModelComplexAttributes>
	{
	private:
		inline void init()
		{
			this->U_OVER_N = U / Constants::BASIS_SIZE;
			this->V_OVER_N = V / Constants::BASIS_SIZE;

			this->SPINOR_SIZE = 4;
			computeChemicalPotential();
		};

	protected:
		using BaseAttributes = std::conditional_t<std::is_same_v<DataType, double_prec>, BaseModelRealAttributes, BaseModelComplexAttributes>;
		typedef Eigen::Vector<DataType, Eigen::Dynamic> ParameterVector;

		const complex_prec I = { 0, 1 };
		SpinorMatrix hamilton;

		double_prec temperature;
		double_prec U;
		double_prec V;
		double_prec U_OVER_N;
		double_prec V_OVER_N;
		double_prec chemical_potential;

		size_t TOTAL_BASIS;
		size_t SPINOR_SIZE;

		// Stores the coefficients for the parameters (e.g. V/N) with the appropriate index
		std::vector<double_prec> parameterCoefficients;

		inline virtual void computeChemicalPotential() {
			this->chemical_potential = 0.5 * U + 4 * V;
		};

		inline double_prec fermi_dirac(double_prec energy) const {
			if (temperature > 1e-8) {
				return (1. / (1. + exp(energy / temperature)));
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
		void setParameters(ParameterVector& F) {
			for (size_t i = 0; i < F.size(); i++)
			{
				F(i) *= parameterCoefficients[i];
			}

			constexpr double_prec new_weight = 0.5;
			for (size_t i = 0; i < F.size(); i++)
			{
				*(this->parameterMapper[i]) = new_weight * F(i) + (1 - new_weight) * (*(this->parameterMapper[i]));
			}
		};

		virtual void fillHamiltonianHelper(va_list args) = 0;
		void fillHamiltonian(int variadic_count, ...) {
			va_list args;
			va_start(args, variadic_count);
			fillHamiltonianHelper(args);
			va_end(args);
		};

		virtual void addToParameterSetHelper(const SpinorMatrix& rho, ComplexParameterVector& F, va_list args) = 0;
		void addToParameterSet(const SpinorMatrix& rho, ComplexParameterVector& F, int variadic_count, ...) {
			va_list args;
			va_start(args, variadic_count);
			addToParameterSetHelper(rho, F, args);
			va_end(args);
		};

	public:
		BaseModel(const ModelParameters& _params)
			: BaseAttributes(_params), temperature(_params.temperature), U(_params.U), V(_params.V)
		{
			init();
		};

		BaseModel(const ModelParameters& _params, const BaseAttributes& startingValues)
			: BaseAttributes(startingValues), temperature(_params.temperature), U(_params.U), V(_params.V)
		{
			init();
		};

		virtual BaseModelRealAttributes computePhases(const bool print = false) = 0;

		inline virtual double_prec entropyPerSite() = 0;
		inline virtual double_prec internalEnergyPerSite() = 0;
		inline virtual double_prec freeEnergyPerSite() {
			return this->internalEnergyPerSite() - temperature * this->entropyPerSite();
		};
	};
}