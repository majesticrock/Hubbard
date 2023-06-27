#include "HubbardCDW.hpp"

constexpr size_t NUMBER_OF_PARAMETERS = 16;

namespace Hubbard::SquareLattice {
	void HubbardCDW::init()
	{
		this->delta_cdw = -0. * (0.5 * U - 4. * V);
		this->delta_afm = 0.1 * 0.5 * U;
		this->delta_sc = -0.3 * U;

		if (V > 0) {
			this->delta_sc *= 0;
		}
		//else if (V < 0) { -1.24153679
		//	this->delta_cdw *= 0;
		//}
		this->delta_eta = 0.4 * I * U;
		this->gamma_occupation_up = V * 0.2;
		this->gamma_occupation_down = V * 0.2;
		this->xi_occupation_up = V * 0.1;
		this->xi_occupation_down = V * 0.1;
		this->gamma_sc = V * 0.1;
		this->xi_sc = I * std::abs(V) * 0.15;

		this->gamma_cdw = I * V * 0.;
		this->xi_cdw = I * V * 0.;
		this->gamma_afm = I * V * 0.;
		this->xi_afm = I * V * 0.;
		this->gamma_eta = I * V * 0.;
		this->xi_eta = V * 0.;

		this->hamilton = SpinorMatrix::Zero(4, 4);

		parameterMapper = {
			&(this->delta_cdw),
			&(this->delta_afm),
			&(this->delta_sc),
			&(this->gamma_sc),
			&(this->xi_sc),
			&(this->delta_eta),
			&(this->gamma_occupation_up),
			&(this->gamma_occupation_down),
			&(this->gamma_cdw),
			&(this->xi_cdw),
			&(this->gamma_afm),
			&(this->xi_afm),
			&(this->xi_occupation_up),
			&(this->xi_occupation_down),
			&(this->gamma_eta),
			&(this->xi_eta)
		};

		parameterCoefficients = {
			0.5 * U_OVER_N - 4. * V_OVER_N, // CDW
			0.5 * U_OVER_N, // AFM
			U_OVER_N, // SC
			V_OVER_N, // Gamma SC
			V_OVER_N, // Xi SC
			U_OVER_N, // Eta
			V_OVER_N, // Gamma Occupation up
			V_OVER_N, // Gamma Occupation down
			0.5 * V_OVER_N, // Gamma CDW
			0.5 * V_OVER_N, // Xi CDW
			0.5 * V_OVER_N, // Gamma AFM
			0.5 * V_OVER_N, // Xi AFM
			V_OVER_N, // Xi Occupation Up
			V_OVER_N, // Xi Occupation Down
			V_OVER_N, // Gamma SC
			V_OVER_N // Xi SC
		};
	}
	void HubbardCDW::fillHamiltonianHelper(va_list args)
	{
		UNPACK_2D;
		hamilton.fill(0.0);
		const double_prec GAMMA = gamma(k_x, k_y);
		const double_prec XI = xi(k_x, k_y);

		hamilton(0, 1) = this->delta_cdw - this->delta_afm;// + ((gamma_cdw - this->gamma_afm) * GAMMA + (xi_cdw - this->xi_afm) * XI);

		hamilton(0, 2) = this->delta_sc + (this->gamma_sc * GAMMA + this->xi_sc * XI);
		hamilton(0, 3) = this->delta_eta;// +(this->gamma_eta * GAMMA + this->xi_eta * XI);

		hamilton(1, 2) = this->delta_eta;// - (this->gamma_eta * GAMMA + this->xi_eta * XI);
		hamilton(1, 3) = this->delta_sc - (this->gamma_sc * GAMMA + this->xi_sc * XI);

		hamilton(2, 3) = -this->delta_cdw - this->delta_afm;// - ((gamma_cdw + this->gamma_afm) * GAMMA + (xi_cdw + this->xi_afm) * XI);

		SpinorMatrix buffer = hamilton.adjoint();
		hamilton += buffer;
		double_prec eps = renormalizedEnergy_up(GAMMA) - 2. * this->xi_occupation_up.real() * XI;
		hamilton(0, 0) = eps;
		hamilton(1, 1) = -eps;
		eps = renormalizedEnergy_down(GAMMA) - 2. * this->xi_occupation_down.real() * XI;
		hamilton(2, 2) = -eps;
		hamilton(3, 3) = eps;
	}
	HubbardCDW::HubbardCDW(const ModelParameters& _params)
		: Model(_params)
	{
		init();
	}
	HubbardCDW::HubbardCDW(const ModelParameters& _params, const BaseAttributes& startingValues)
		: Model(_params, startingValues)
	{
		init();
	}
	BaseModelRealAttributes HubbardCDW::computePhases(const bool print)
	{
		computeChemicalPotential();

		constexpr double_prec EPSILON = 1e-12;
		double_prec error = 100;
		constexpr int MAX_STEPS = 1500;

		ParameterVector f0 = ParameterVector(NUMBER_OF_PARAMETERS);
		for (size_t i = 0; i < NUMBER_OF_PARAMETERS; i++)
		{
			f0(i) = *(parameterMapper[i]);
		}

		ParameterVector x0 = f0;

		if (print) {
			std::cout << "-1:\t" << std::fixed << std::setprecision(8);
			printAsRow(x0);
		}
		for (size_t i = 0; i < MAX_STEPS && error > EPSILON; i++)
		{
			iterationStep(x0, f0);
			error = f0.norm();
			for (size_t i = 0; i < NUMBER_OF_PARAMETERS; i++)
			{
				x0(i) = *(parameterMapper[i]);
			}
			if (print) {
				std::cout << i << ":\t" << std::fixed << std::setprecision(8);
				printAsRow(x0);
			}
			if (i == MAX_STEPS - 1) {
				std::cerr << "[T, U, V] = [" << this->temperature << ", " << this->U << "," << this->V
					<< "]\tConvergence at " << error << std::endl;
				delta_cdw = 0;
				delta_afm = 0;
				delta_sc = 0;
				delta_eta = 0;
			}
		}

		if (std::abs(delta_sc.imag()) > 1e-8) {
			std::cout << "[T, U, V] = [" << this->temperature << ", " << this->U << "," << this->V
				<< "]" << std::endl;
		}
		if (std::abs(gamma_sc.real()) > 1e-8) {
			std::cout << "[T, U, V] = [" << this->temperature << ", " << this->U << "," << this->V
				<< "]" << std::endl;
		}
		if (std::abs(xi_sc.real()) > 1e-8) {
			std::cout << "[T, U, V] = [" << this->temperature << ", " << this->U << "," << this->V
				<< "]" << std::endl;
		}
		if (std::abs(delta_eta) > 1e-8) {
			std::cout << "[T, U, V] = [" << this->temperature << ", " << this->U << "," << this->V
				<< "]" << std::endl;
		}
		if (std::abs(gamma_eta) > 1e-8) {
			std::cout << "[T, U, V] = [" << this->temperature << ", " << this->U << "," << this->V
				<< "]" << std::endl;
		}
		if (std::abs(xi_eta) > 1e-8) {
			std::cout << "[T, U, V] = [" << this->temperature << ", " << this->U << "," << this->V
				<< "]" << std::endl;
		}

		return BaseModelRealAttributes(*this);
	}
}