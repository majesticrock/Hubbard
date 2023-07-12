#pragma once
#include "../MomentumBasedModel.hpp"
#include <numeric>

#define DELTA_CDW this->model_attributes[0]
#define DELTA_AFM this->model_attributes[1]
#define DELTA_SC this->model_attributes[2]
#define GAMMA_SC this->model_attributes[3]
#define XI_SC this->model_attributes[4]
#define DELTA_ETA this->model_attributes[5]
#define GAMMA_OCCUPATION_UP this->model_attributes[6]
#define GAMMA_OCCUPATION_DOWN this->model_attributes[7]

namespace Hubbard::SquareLattice
{
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
	class Model2D : public MomentumBasedModel<DataType, 2>
	{
	private:
		void init() {
			this->parameterCoefficients = {
				0.5 * this->U_OVER_N - 4. * this->V_OVER_N, // CDW
				0.5 * this->U_OVER_N, // AFM
				this->U_OVER_N, // SC
				this->V_OVER_N, // Gamma SC
				this->V_OVER_N, // Xi SC
				this->U_OVER_N, // Eta
				this->V_OVER_N, // Occupation Up
				this->V_OVER_N, // Occupation Down
			};
		};
	protected:
		using ParameterVector = typename BaseModel<DataType>::ParameterVector;

		inline void iterationStep(const ParameterVector& x, ParameterVector& F) {
			F.fill(0);
			std::conditional_t<std::is_same_v<DataType, complex_prec>,
				ComplexParameterVector&, ComplexParameterVector> complex_F = F;

			SpinorMatrix rho{ SpinorMatrix::Zero(this->SPINOR_SIZE, this->SPINOR_SIZE) };
			Eigen::SelfAdjointEigenSolver<SpinorMatrix> solver;

			std::copy(x.begin(), x.end(), this->model_attributes.selfconsistency_values.begin());

			for (int k = -Constants::K_DISCRETIZATION; k < Constants::K_DISCRETIZATION; ++k)
			{
				double k_x = (k * L_PI) / Constants::K_DISCRETIZATION;
				for (int l = -Constants::K_DISCRETIZATION; l < 0; ++l)
				{
					double k_y = (l * L_PI) / Constants::K_DISCRETIZATION;
					this->fillHamiltonian(2, k_x, k_y);
					solver.compute(this->hamilton);
					this->fillRho(rho, solver);

					this->addToParameterSet(rho, complex_F, 2, k_x, k_y);
				}
			}

			/* { // Checks for numerical accurarcy
				const double ERROR_MARGIN = 1e-10 * Constants::BASIS_SIZE;
				if (std::abs(complex_F(2).imag()) > ERROR_MARGIN) {
					std::cout << "sc: " << complex_F(2) << "\t Params: " << this->temperature << ", " << this->U << ", " << this->V << std::endl;
				}
				if (std::abs(complex_F(3).imag()) > ERROR_MARGIN) {
					std::cout << "gamma sc: " << complex_F(3) << "\t Params: " << this->temperature << ", " << this->U << ", " << this->V << std::endl;
				}
				if (std::abs(complex_F(4).real()) > ERROR_MARGIN) {
					std::cout << "xi sc: " << complex_F(4) << "\t Params: " << this->temperature << ", " << this->U << ", " << this->V << std::endl;
				}
				if (std::abs(complex_F(5).real()) > ERROR_MARGIN) {
					std::cout << "eta: " << complex_F(5) << "\t Params: " << this->temperature << ", " << this->U << ", " << this->V << std::endl;
				}
			}*/

			if constexpr (!std::is_same_v<DataType, complex_prec>) {
				complexParametersToReal(complex_F, F);
			}
			this->multiplyParametersByCoefficients(F);
			for (auto& value : F) {
				if (std::abs(value) < 1e-12) value = 0.;
			}
			this->setParameters(F);

			F -= x;
		};
	public:
		Model2D(const ModelParameters& _params)
			: MomentumBasedModel<DataType, 2>(_params)
		{
			this->init();
		};

		template<typename StartingValuesDataType>
		Model2D(const ModelParameters& _params, const ModelAttributes<StartingValuesDataType>& startingValues)
			: MomentumBasedModel<DataType, 2>(_params, startingValues)
		{
			this->init();
		};

		// saves all one particle energies to reciever
		inline void getAllEnergies(std::vector<std::vector<double>>& reciever)
		{
			reciever.resize(4 * Constants::K_DISCRETIZATION, std::vector<double>(2 * Constants::K_DISCRETIZATION));
			Eigen::SelfAdjointEigenSolver<SpinorMatrix> solver;
			double k_val = 0;
			double l_val = 0;
			for (int k = -Constants::K_DISCRETIZATION; k < Constants::K_DISCRETIZATION; k++)
			{
				k_val = k * L_PI / Constants::K_DISCRETIZATION;
				for (int l = -Constants::K_DISCRETIZATION; l < Constants::K_DISCRETIZATION; l++)
				{
					l_val = l * L_PI / Constants::K_DISCRETIZATION;
					this->fillHamiltonian(2, k_val, l_val);
					solver.compute(this->hamilton, false);
					reciever[k + Constants::K_DISCRETIZATION][l + Constants::K_DISCRETIZATION] = solver.eigenvalues()(0);

					for (int i = 1; i < 4; i++)
					{
						if (std::abs(solver.eigenvalues()(0) - solver.eigenvalues()(i)) > 1e-8) {
							reciever[k + 3 * Constants::K_DISCRETIZATION][l + Constants::K_DISCRETIZATION] = solver.eigenvalues()(i);
							break;
						}
					}
				}
			}
		};

		inline virtual double entropyPerSite() override {
			double entropy = 0;
			Eigen::SelfAdjointEigenSolver<SpinorMatrix> solver;
			for (int k = -Constants::K_DISCRETIZATION; k < Constants::K_DISCRETIZATION; ++k)
			{
				double k_x = (k * L_PI) / Constants::K_DISCRETIZATION;
				for (int l = -Constants::K_DISCRETIZATION; l < 0; ++l)
				{
					double k_y = (l * L_PI) / Constants::K_DISCRETIZATION;
					this->fillHamiltonian(2, k_x, k_y);
					solver.compute(this->hamilton, false);

					entropy += std::accumulate(solver.eigenvalues().begin(), solver.eigenvalues().end(), double{},
						[this](double current, double toAdd) {
							auto occ = BaseModel<DataType>::fermi_dirac(toAdd);
							// Let's just not take the ln of 0. Negative numbers cannot be reached (because math...)
							return (occ > 1e-12 ? current - occ * std::log(occ) : current);
						});
				}
			}
			return entropy / Constants::BASIS_SIZE;
		};

		inline virtual double internalEnergyPerSite() override {
			double energy = 0;
			Eigen::SelfAdjointEigenSolver<SpinorMatrix> solver;
			for (int k = -Constants::K_DISCRETIZATION; k < Constants::K_DISCRETIZATION; ++k)
			{
				double k_x = (k * L_PI) / Constants::K_DISCRETIZATION;
				for (int l = -Constants::K_DISCRETIZATION; l < 0; ++l)
				{
					double k_y = (l * L_PI) / Constants::K_DISCRETIZATION;
					this->fillHamiltonian(2, k_x, k_y);
					solver.compute(this->hamilton, false);

					energy += std::accumulate(solver.eigenvalues().begin(), solver.eigenvalues().end(), double{},
						[this](double current, double toAdd) {
							return current + toAdd * BaseModel<DataType>::fermi_dirac(toAdd);
						});
				}
			}
			return energy / Constants::BASIS_SIZE;
		};

		inline void computeExpectationValues(std::vector<MatrixCL>& expecs, std::vector<complex_prec>& sum_of_all) {
			SpinorMatrix rho = SpinorMatrix::Zero(this->SPINOR_SIZE, this->SPINOR_SIZE);
			Eigen::SelfAdjointEigenSolver<SpinorMatrix> solver;

			expecs = std::vector<MatrixCL>(8, Matrix_L::Zero(2 * Constants::K_DISCRETIZATION, 2 * Constants::K_DISCRETIZATION));
			sum_of_all = std::vector<std::complex<double>>(8, 0.0);

			for (int k = -Constants::K_DISCRETIZATION; k < Constants::K_DISCRETIZATION; k++)
			{
				for (int l = -Constants::K_DISCRETIZATION; l < Constants::K_DISCRETIZATION; l++)
				{
					this->fillHamiltonian(2, (k * L_PI) / Constants::K_DISCRETIZATION, (l * L_PI) / Constants::K_DISCRETIZATION);
					solver.compute(this->hamilton);
					this->fillRho(rho, solver);

					// n_up
					expecs[0](k + Constants::K_DISCRETIZATION, l + Constants::K_DISCRETIZATION) = 1 - rho(0, 0).real();
					// g_up
					expecs[1](k + Constants::K_DISCRETIZATION, l + Constants::K_DISCRETIZATION) = -rho(1, 0);
					// f
					expecs[2](k + Constants::K_DISCRETIZATION, l + Constants::K_DISCRETIZATION) = -rho(0, 2);
					// eta
					expecs[3](k + Constants::K_DISCRETIZATION, l + Constants::K_DISCRETIZATION) = -rho(0, 3);
					// n_down
					expecs[4](k + Constants::K_DISCRETIZATION, l + Constants::K_DISCRETIZATION) = rho(2, 2).real();
					// g_down
					expecs[5](k + Constants::K_DISCRETIZATION, l + Constants::K_DISCRETIZATION) = rho(2, 3);
					// n_up + n_down
					expecs[6](k + Constants::K_DISCRETIZATION, l + Constants::K_DISCRETIZATION) = 1 - (rho(0, 0) - rho(2, 2)).real();
					// g_up + g_down
					expecs[7](k + Constants::K_DISCRETIZATION, l + Constants::K_DISCRETIZATION) = rho(2, 3) - rho(1, 0);
					for (size_t idx = 0U; idx < 8U; ++idx)
					{
						sum_of_all[idx] += expecs[idx](k + Constants::K_DISCRETIZATION, l + Constants::K_DISCRETIZATION);
					}

					if (std::abs(rho(3, 0)) > 1e-10) {
						std::cerr << "Warning: <eta> does not vanish! " << rho(3, 0) << std::endl;
					}
				}
			}
		};
	};
}