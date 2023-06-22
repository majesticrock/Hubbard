#pragma once
#include "../MomentumBasedModel.hpp"

namespace Hubbard::SquareLattice {
	template <typename DataType>
	class Model : public MomentumBasedModel<DataType, 2>
	{
	protected:
		using ParameterVector = typename BaseModel<DataType>::ParameterVector;

		inline void iterationStep(const ParameterVector& x, ParameterVector F) {
			F.fill(0);
			SpinorMatrix rho = SpinorMatrix::Zero(this->SPINOR_SIZE, this->SPINOR_SIZE);
			Eigen::SelfAdjointEigenSolver<SpinorMatrix> solver;

			for (size_t i = 0; i < this->parameterMapper.size(); i++)
			{
				*(this->parameterMapper[i]) = x(i);
			}

			for (int k = -Constants::K_DISCRETIZATION; k < Constants::K_DISCRETIZATION; k++)
			{
				double_prec k_x = (k * L_PI) / Constants::K_DISCRETIZATION;
				for (int l = -Constants::K_DISCRETIZATION; l < 0; l++)
				{
					double_prec k_y = (l * L_PI) / Constants::K_DISCRETIZATION;
					this->fillHamiltonian(k_x, k_y);
					solver.compute(this->hamilton);
					this->fillRho(rho, solver);

					this->addToParameterSet(rho, F, k_x, k_y);
				}
			}

			this->setParameters(F);
			F -= x;
		};
	public:
		Model(const ModelParameters& _params) : MomentumBasedModel<DataType, 2>(_params) { };

		// saves all one particle energies to reciever
		inline void getAllEnergies(std::vector<std::vector<double>>& reciever)
		{
			reciever.resize(4 * Constants::K_DISCRETIZATION, std::vector<double>(2 * Constants::K_DISCRETIZATION));
			Eigen::SelfAdjointEigenSolver<SpinorMatrix> solver;
			double_prec k_val = 0;
			double_prec l_val = 0;
			for (int k = -Constants::K_DISCRETIZATION; k < Constants::K_DISCRETIZATION; k++)
			{
				k_val = k * L_PI / Constants::K_DISCRETIZATION;
				for (int l = -Constants::K_DISCRETIZATION; l < Constants::K_DISCRETIZATION; l++)
				{
					l_val = l * L_PI / Constants::K_DISCRETIZATION;
					this->fillHamiltonian(k_val, l_val);
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

		inline void computeExpectationValues(std::vector<MatrixCL>& expecs, std::vector<complex_prec>& sum_of_all) {
			SpinorMatrix rho = SpinorMatrix::Zero(4, 4);
			Eigen::SelfAdjointEigenSolver<SpinorMatrix> solver;

			expecs = std::vector<MatrixCL>(8, Matrix_L::Zero(2 * Constants::K_DISCRETIZATION, 2 * Constants::K_DISCRETIZATION));
			sum_of_all = std::vector<std::complex<double>>(8, 0.0);

			for (int k = -Constants::K_DISCRETIZATION; k < Constants::K_DISCRETIZATION; k++)
			{
				for (int l = -Constants::K_DISCRETIZATION; l < Constants::K_DISCRETIZATION; l++)
				{
					this->fillHamiltonian((k * L_PI) / Constants::K_DISCRETIZATION, (l * L_PI) / Constants::K_DISCRETIZATION);
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
					for (int idx = 0; idx < 8; idx++)
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