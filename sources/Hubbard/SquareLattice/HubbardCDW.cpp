#include "HubbardCDW.hpp"

#define GAMMA_CDW this->model_attributes[8]
#define XI_CDW this->model_attributes[9]
#define GAMMA_AFM this->model_attributes[10]
#define XI_AFM this->model_attributes[11]
#define XI_OCCUPATiON_UP this->model_attributes[12]
#define XI_OCCUPATiON_DOWN this->model_attributes[13]
#define GAMMA_ETA this->model_attributes[14]
#define XI_ETA this->model_attributes[15]

namespace Hubbard::SquareLattice {
	void HubbardCDW::init()
	{
		DELTA_ETA *= I;
		XI_SC *= I;
		GAMMA_SC = V * 0.2;

		model_attributes.selfconsistency_values.reserve(16);

		model_attributes.push_back(I * V * 0.);
		model_attributes.push_back(I * V * 0.);
		model_attributes.push_back(I * V * 0.);
		model_attributes.push_back(I * V * 0.);
		model_attributes.push_back(V * 0.);
		model_attributes.push_back(V * 0.);
		model_attributes.push_back(I * V * 0.);
		model_attributes.push_back(V * 0.);

		this->hamilton = SpinorMatrix::Zero(4, 4);

		parameterCoefficients = {
			0.5 * U_OVER_N - 4. * V_OVER_N, // CDW - 0
			0.5 * U_OVER_N, // AFM - 1
			U_OVER_N, // SC - 2
			V_OVER_N, // Gamma SC - 3
			V_OVER_N, // Xi SC - 4
			U_OVER_N, // Eta - 5
			V_OVER_N, // Gamma Occupation up   - 6
			V_OVER_N, // Gamma Occupation down - 7
			0.5 * V_OVER_N, // Gamma CDW - 8
			0.5 * V_OVER_N, // Xi CDW - 9
			0.5 * V_OVER_N, // Gamma AFM - 10
			0.5 * V_OVER_N, // Xi AFM - 11
			V_OVER_N, // Xi Occupation Up   - 12
			V_OVER_N, // Xi Occupation Down - 13
			V_OVER_N, // Gamma ETA - 14
			V_OVER_N // Xi ETA - 15
		};
	}
	void HubbardCDW::fillHamiltonian(const NumericalMomentum<2>& k_values)
	{
		hamilton.fill(0.0);
		const double GAMMA = k_values.gamma();
		const double XI = xi(k_values);

		hamilton(0, 1) = DELTA_CDW - DELTA_AFM + ((GAMMA_CDW - GAMMA_AFM) * GAMMA + (XI_CDW - XI_AFM) * XI);

		hamilton(0, 2) = DELTA_SC + (GAMMA_SC * GAMMA + XI_SC * XI);
		hamilton(0, 3) = DELTA_ETA + (GAMMA_ETA * GAMMA + XI_ETA * XI);

		hamilton(1, 2) = DELTA_ETA - (GAMMA_ETA * GAMMA + XI_ETA * XI);
		hamilton(1, 3) = DELTA_SC - (GAMMA_SC * GAMMA + XI_SC * XI);

		hamilton(2, 3) = -DELTA_CDW - DELTA_AFM - ((GAMMA_CDW + GAMMA_AFM) * GAMMA + (XI_CDW + XI_AFM) * XI);

		SpinorMatrix buffer{ hamilton.adjoint() };
		hamilton += buffer;
		double eps = model_attributes.renormalizedEnergy_up(GAMMA) - model_attributes[12].real() * XI;
		hamilton(0, 0) = eps;
		hamilton(1, 1) = -eps;
		eps = model_attributes.renormalizedEnergy_down(GAMMA) - model_attributes[13].real() * XI;
		hamilton(2, 2) = -eps;
		hamilton(3, 3) = eps;
	}

	void HubbardCDW::addToParameterSet(const SpinorMatrix& rho, ParameterVector& F, const NumericalMomentum<2>& k_values)
	{
			const double GAMMA = k_values.gamma();
			const double XI = xi(k_values);

			F(0) -= (rho(0, 1) + rho(1, 0) - rho(2, 3) - rho(3, 2)).real(); // CDW
			F(1) -= (rho(0, 1) + rho(1, 0) + rho(2, 3) + rho(3, 2)).real(); // AFM
			F(2) -= (rho(0, 2) + rho(1, 3)); // SC
			F(3) -= GAMMA * (rho(0, 2) - rho(1, 3)); // Gamma SC
			F(4) -= XI * (rho(0, 2) - rho(1, 3)); // Xi SC
			F(5) -= (rho(0, 3) + rho(1, 2)); // Eta
			F(6) -= GAMMA * (rho(0, 0) - rho(1, 1)).real(); // Gamma Occupation Up
			F(7) += GAMMA * (rho(2, 2) - rho(3, 3)).real(); // Gamma Occupation Down
			F(8) += I * GAMMA * (rho(0, 1) - rho(1, 0) + rho(2, 3) - rho(3, 2)).imag(); // Gamma CDW
			F(9) += I * XI * (rho(0, 1) - rho(1, 0) + rho(2, 3) - rho(3, 2)).imag(); // Xi CDW
			F(10) += I * GAMMA * (rho(0, 1) - rho(1, 0) - rho(2, 3) + rho(3, 2)).imag(); // Gamma AFM
			F(11) += I * XI * (rho(0, 1) - rho(1, 0) - rho(2, 3) + rho(3, 2)).imag(); // Xi AFM
			F(12) -= XI * (rho(0, 0) - rho(1, 1)).real(); // Xi Occupation Up
			F(13) += XI * (rho(2, 2) - rho(3, 3)).real(); // Xi Occupation Down
			F(14) -= GAMMA * (rho(0, 3) - rho(1, 2)); // Gamma eta
			F(15) -= XI * (rho(0, 3) - rho(1, 2)); // Xi eta
		}

	HubbardCDW::HubbardCDW(const ModelParameters& _params)
		: Model2D(_params)
	{
		init();
	}
	ModelAttributes<double> HubbardCDW::computePhases(const PhaseDebuggingPolicy debugPolicy/*=PhaseDebuggingPolicy{}*/)
	{
		constexpr double EPSILON = 1e-12;
		double error = 100;
		constexpr size_t MAX_STEPS = 1000;
		const size_t NUMBER_OF_PARAMETERS = this->model_attributes.size();

		ParameterVector f0{ ParameterVector::Zero(NUMBER_OF_PARAMETERS) };
		std::copy(model_attributes.selfconsistency_values.begin(), model_attributes.selfconsistency_values.end(), f0.begin());

		ParameterVector x0 = f0;

		if (debugPolicy.printAll) {
			std::cout << "-1:\t" << std::fixed << std::setprecision(8);
			printAsRow<-1>(x0);
		}
		model_attributes.converged = true;
		for (size_t i = 0U; i < MAX_STEPS && error > EPSILON; ++i)
		{
			iterationStep(x0, f0);
			error = f0.norm();
			std::copy(model_attributes.selfconsistency_values.begin(), model_attributes.selfconsistency_values.end(), x0.begin());

			if (debugPolicy.printAll) {
				std::cout << i << ":\t" << std::fixed << std::setprecision(8);
				printAsRow<-1>(x0);
			}
			if (i == MAX_STEPS - 1) {
				if (debugPolicy.convergenceWarning){
					std::cerr << "No convergence for [T U V] = [" << std::fixed << std::setprecision(8)
					<< this->temperature << " " << this->U << " " << this->V << "]" << std::endl;
				}

				std::fill(model_attributes.selfconsistency_values.begin(), model_attributes.selfconsistency_values.end(), 0.);
				model_attributes.converged = false;
			}
		}

		return ModelAttributes<double>(this->model_attributes);
	}
}