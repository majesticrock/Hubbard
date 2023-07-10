#include "HubbardCDW.hpp"

constexpr size_t NUMBER_OF_PARAMETERS = 16;

namespace Hubbard::SquareLattice {
	void HubbardCDW::init()
	{
		DELTA_CDW = -0. * (0.5 * U - 4. * V);
		DELTA_AFM = 0.1 * 0.5 * U;
		DELTA_SC = -0.3 * U;

		if (V > 0) {
			DELTA_SC *= 0;
		}
		//else if (V < 0) { -1.24153679
		//	delta_cdw *= 0;
		//}
		DELTA_ETA = 0.4 * I * U;
		GAMMA_OCCUPATION_UP = V * 0.2;
		GAMMA_OCCUPATION_DOWN = V * 0.2;
		GAMMA_SC = V * 0.1;
		XI_SC = I * std::abs(V) * 0.15;

		model_attributes.selfconsistency_values.reserve(16);
		model_attributes.push_back(V * 0.1);

		model_attributes.push_back(I * V * 0.);
		model_attributes.push_back(I * V * 0.);
		model_attributes.push_back(I * V * 0.);
		model_attributes.push_back(I * V * 0.);
		model_attributes.push_back(V * 0.1);
		model_attributes.push_back(V * 0.1);
		model_attributes.push_back(I * V * 0.);
		model_attributes.push_back(V * 0.);

		this->hamilton = SpinorMatrix::Zero(4, 4);

		parameterCoefficients = {
			0.5 * U_OVER_N - 4. * V_OVER_N, // CDW - 0
			0.5 * U_OVER_N, // AFM - 1
			U_OVER_N, // SC - 2
			V_OVER_N, // Gamma SC - 3
			V_OVER_N, // Xi SC - 4
			U_OVER_N, // Eta -
			V_OVER_N, // Gamma Occupation up   - 6
			V_OVER_N, // Gamma Occupation down - 7
			0.5 * V_OVER_N, // Gamma CDW - 8
			0.5 * V_OVER_N, // Xi CDW - 9
			0.5 * V_OVER_N, // Gamma AFM - 10
			0.5 * V_OVER_N, // Xi AFM - 11
			V_OVER_N, // Xi Occupation Up   - 12
			V_OVER_N, // Xi Occupation Down - 13
			V_OVER_N, // Gamma SC - 14
			V_OVER_N // Xi SC - 15
		};
	}
	void HubbardCDW::fillHamiltonianHelper(va_list args)
	{
		UNPACK_2D;
		hamilton.fill(0.0);
		const double GAMMA = gamma(k_x, k_y);
		const double XI = xi(k_x, k_y);

		hamilton(0, 1) = DELTA_CDW - DELTA_AFM;// + ((gamma_cdw - this->gamma_afm) * GAMMA + (xi_cdw - this->xi_afm) * XI);

		hamilton(0, 2) = DELTA_SC + (GAMMA_SC * GAMMA + XI_SC * XI);
		hamilton(0, 3) = DELTA_ETA;// +(this->gamma_eta * GAMMA + this->xi_eta * XI);

		hamilton(1, 2) = DELTA_ETA;// - (this->gamma_eta * GAMMA + this->xi_eta * XI);
		hamilton(1, 3) = DELTA_SC - (GAMMA_SC * GAMMA + XI_SC * XI);

		hamilton(2, 3) = -DELTA_CDW - DELTA_AFM;// - ((gamma_cdw + this->gamma_afm) * GAMMA + (xi_cdw + this->xi_afm) * XI);

		SpinorMatrix buffer = hamilton.adjoint();
		hamilton += buffer;
		double eps = model_attributes.renormalizedEnergy_up(GAMMA) - 2. * model_attributes[12].real() * XI;
		hamilton(0, 0) = eps;
		hamilton(1, 1) = -eps;
		eps = model_attributes.renormalizedEnergy_down(GAMMA) - 2. * model_attributes[13].real() * XI;
		hamilton(2, 2) = -eps;
		hamilton(3, 3) = eps;
	}
	HubbardCDW::HubbardCDW(const ModelParameters& _params)
		: Model2D(_params)
	{
		init();
	}
	ModelAttributes<double> HubbardCDW::computePhases(const bool print)
	{
		computeChemicalPotential();

		constexpr double EPSILON = 1e-12;
		double error = 100;
		constexpr int MAX_STEPS = 1500;

		ParameterVector f0 = ParameterVector(NUMBER_OF_PARAMETERS);
		for (size_t i = 0; i < NUMBER_OF_PARAMETERS; i++)
		{
			f0(i) = model_attributes[i];
		}

		ParameterVector x0 = f0;

		if (print) {
			std::cout << "-1:\t" << std::fixed << std::setprecision(8);
			printAsRow<-1>(x0);
		}
		model_attributes.converged = true;
		for (size_t i = 0U; i < MAX_STEPS && error > EPSILON; ++i)
		{
			iterationStep(x0, f0);
			error = f0.norm();
			for (size_t j = 0U; j < NUMBER_OF_PARAMETERS; ++j)
			{
				x0(j) = model_attributes[j];
			}
			if (print) {
				std::cout << i << ":\t" << std::fixed << std::setprecision(8);
				printAsRow<-1>(x0);
			}
			if (i == MAX_STEPS - 1) {
				std::cerr << "[T, U, V] = [" << this->temperature << ", " << this->U << "," << this->V
					<< "]\tConvergence at " << error << std::endl;
				for (auto& value : model_attributes.selfconsistency_values) {
					value = 0.;
				}
				model_attributes.converged = false;
			}
		}

		if (std::abs(DELTA_SC.imag()) > 1e-8) {
			std::cout << "[T, U, V] = [" << this->temperature << ", " << this->U << "," << this->V
				<< "]" << std::endl;
		}
		if (std::abs(GAMMA_SC.real()) > 1e-8) {
			std::cout << "[T, U, V] = [" << this->temperature << ", " << this->U << "," << this->V
				<< "]" << std::endl;
		}
		if (std::abs(XI_SC.real()) > 1e-8) {
			std::cout << "[T, U, V] = [" << this->temperature << ", " << this->U << "," << this->V
				<< "]" << std::endl;
		}
		if (std::abs(DELTA_ETA) > 1e-8) {
			std::cout << "[T, U, V] = [" << this->temperature << ", " << this->U << "," << this->V
				<< "]" << std::endl;
		}

		return ModelAttributes<double>(this->model_attributes);
	}
}