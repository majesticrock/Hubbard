#pragma once
#include "../Utility/Resolvent.hpp"
#include "Constants.hpp"
#include "ModelAttributes.hpp"

namespace Hubbard {
	struct PhaseDebuggingPolicy{
		bool printAll{false};
		bool convergenceWarning{true};

		PhaseDebuggingPolicy() = default;
		PhaseDebuggingPolicy(bool _printAll, bool _convergenceWarning) 
			: printAll{_printAll}, convergenceWarning(_convergenceWarning) {} ;
	};

	constexpr double L_PI = 3.141592653589793238462643383279502884L; /* pi */
	typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Matrix_L;
	typedef Eigen::Vector<double, Eigen::Dynamic> Vector_L;

	typedef std::complex<double> complex_prec;
	typedef Eigen::Matrix<complex_prec, Eigen::Dynamic, Eigen::Dynamic> MatrixCL;
	typedef Eigen::Vector<complex_prec, Eigen::Dynamic> VectorCL;
	using SpinorMatrix = MatrixCL;
	typedef Utility::Resolvent<double> Resolvent_L;

	template <const int vector_size = Eigen::Dynamic>
	void printAsRow(Eigen::Vector<double, vector_size>& printer) {
		for (size_t i = 0U; i < printer.size(); ++i)
		{
			std::cout << " \t" << printer(i);
			if ((i + 1U) % 8U == 0U) {
				std::cout << "\n\t    ";
			}
		}
		std::cout << std::endl;
	}

	template <const int vector_size = Eigen::Dynamic>
	void printAsRow(Eigen::Vector<complex_prec, vector_size>& printer) {
		for (size_t i = 0U; i < printer.size(); ++i)
		{
			std::cout << " \t" << printer(i);
			if ((i + 1U) % 4U == 0U) {
				std::cout << "\n\t    ";
			}
		}
		std::cout << std::endl;
	}
	const complex_prec I = { 0, 1 };

	using ComplexParameterVector = Eigen::Vector<complex_prec, Eigen::Dynamic>;

	inline void complexParametersToReal(const ComplexParameterVector& c, Eigen::VectorXd& r) {
		r(0) = c(0).real(); // CDW
		r(1) = c(1).real(); // AFM
		r(2) = c(2).real(); // SC
		r(3) = c(3).real(); // Gamma SC
		r(4) = c(4).imag(); // Xi SC
		r(5) = c(5).imag(); // Eta
		r(6) = c(6).real(); // Gamma Occupation Up
		r(7) = c(7).real(); // Gamma Occupation Down
	};

	template <typename DataType>
	class BaseModel
	{
	private:
		inline void init()
		{
			computeChemicalPotential();
		};

	protected:
		using ParameterVector = Eigen::Vector<DataType, Eigen::Dynamic>;
		using BaseAttributes = ModelAttributes<DataType>;

		ModelAttributes<DataType> model_attributes;
		// Stores the coefficients for the parameters (e.g. V/N) with the appropriate index
		std::vector<double> parameterCoefficients;

		SpinorMatrix hamilton;

		double temperature{};
		double U{};
		double V{};
		double U_OVER_N{ U / Constants::BASIS_SIZE };
		double V_OVER_N{ V / Constants::BASIS_SIZE };
		double chemical_potential{};

		size_t TOTAL_BASIS{};
		size_t SPINOR_SIZE{ 4 };

		inline virtual void computeChemicalPotential() {
			this->chemical_potential = 0.5 * U + 4 * V;
		};

		inline double fermi_dirac(double energy) const {
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
		inline void multiplyParametersByCoefficients(ParameterVector& F) const {
			for (size_t i = 0U; i < F.size(); ++i)
			{
				F(i) *= parameterCoefficients[i];
			}
		};
		inline void setParameters(ParameterVector& F) {
			constexpr double new_weight = 0.5;
			for (size_t i = 0U; i < F.size(); ++i)
			{
				this->model_attributes[i] = new_weight * F(i) + (1 - new_weight) * this->model_attributes[i];
			}
		};

	public:
		explicit BaseModel(const ModelParameters& _params)
			: model_attributes(_params), temperature(_params.temperature), U(_params.U), V(_params.V)
		{
			init();
		};

		template<typename StartingValuesDataType>
		BaseModel(const ModelParameters& _params, const ModelAttributes< StartingValuesDataType >& startingValues)
			: model_attributes(startingValues), temperature(_params.temperature), U(_params.U), V(_params.V)
		{
			init();
		};
		virtual ~BaseModel() = default;

		virtual ModelAttributes<double> computePhases(const PhaseDebuggingPolicy debugPolicy=PhaseDebuggingPolicy{}) = 0;

		inline double getTotalGapValue() {
			return this->model_attributes.getTotalGapValue();
		}
		inline virtual double entropyPerSite() = 0;
		inline virtual double internalEnergyPerSite() = 0;
		inline double freeEnergyPerSite() {
			return this->internalEnergyPerSite() - temperature * this->entropyPerSite();
		};
	};
}