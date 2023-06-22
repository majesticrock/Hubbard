#include "HubbardCDW.hpp"

constexpr size_t NUMBER_OF_PARAMETERS = 16;

namespace Hubbard::SquareLattice {
	void HubbardCDW::fillHamiltonianHelper(va_list args)
	{
		UNPACK_2D;
		hamilton.fill(0.0);
		const double_prec GAMMA = gamma(k_x, k_y);
		const double_prec XI = xi(k_x, k_y);

		hamilton(0, 1) = this->delta_cdw - this->delta_afm;// +((gamma_cdw - this->gamma_afm) * GAMMA + (xi_cdw - this->xi_afm) * XI);

		hamilton(0, 2) = this->delta_sc + (this->gamma_sc * GAMMA + this->xi_sc * XI);
		hamilton(0, 3) = this->delta_eta + (this->gamma_eta * GAMMA + this->xi_eta * XI);

		hamilton(1, 2) = this->delta_eta - (this->gamma_eta * GAMMA + this->xi_eta * XI);
		hamilton(1, 3) = this->delta_sc - (this->gamma_sc * GAMMA + this->xi_sc * XI);

		hamilton(2, 3) = -this->delta_cdw - this->delta_afm;// - ((gamma_cdw + this->gamma_afm) * GAMMA + (xi_cdw + this->xi_afm) * XI);

		SpinorMatrix buffer = hamilton.adjoint();
		hamilton += buffer;
		double_prec eps = renormalizedEnergy_up(k_x, k_y);
		hamilton(0, 0) = eps;
		hamilton(1, 1) = -eps;
		eps = renormalizedEnergy_down(k_x, k_y);
		hamilton(2, 2) = -eps;
		hamilton(3, 3) = eps;
	}
	HubbardCDW::HubbardCDW(const ModelParameters& _params)
		: Model(_params), BaseModelComplexAttributes(_params)
	{
		this->delta_cdw = (std::abs(U) + V) * 0.5 + 0.1;
		this->delta_sc = (std::abs(U + std::abs(V)) * 0.3 + 0.05);
		if (V > 0) {
			this->delta_sc *= 0;
		}
		else if (V < 0) {
			this->delta_cdw *= 0;
		}
		if (U > 0) {
			this->delta_afm = std::abs(U - std::abs(V)) * 0.5 + 0.1;
		}
		else {
			this->delta_afm = 0;
		}
		this->delta_eta = 0;//U * 0.1;
		this->delta_occupation_up = -V * 0.2;
		this->delta_occupation_down = -V * 0.2;
		this->delta_occupation_up_y = V * 0.2;
		this->delta_occupation_down_y = V * 0.2;
		this->gamma_sc = V * 0.05;
		this->xi_sc = I * std::abs(V) * 0.1;
		this->gamma_cdw = 0;//I * V * 0.15;
		this->xi_cdw = 0;//I * V * 0.2;
		this->gamma_afm = 0;//I * V * 0.05;
		this->xi_afm = 0;//I * V * 0.04;
		this->gamma_eta = 0;//V * 0.01;
		this->xi_eta = 0;// V * 0.01;

		this->hamilton = SpinorMatrix::Zero(4, 4);

		parameterMapper = {
			&(this->delta_cdw),
			&(this->delta_afm),
			&(this->delta_sc),
			&(this->gamma_sc),
			&(this->xi_sc),
			&(this->delta_eta),
			&(this->delta_occupation_up),
			&(this->delta_occupation_down),
			&(this->gamma_cdw),
			&(this->xi_cdw),
			&(this->gamma_afm),
			&(this->xi_afm),
			&(this->delta_occupation_up_y),
			&(this->delta_occupation_down_y),
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
			V_OVER_N, // Occupation Up
			V_OVER_N, // Occupation Down
			0.5 * V_OVER_N, // Gamma CDW
			0.5 * V_OVER_N, // Xi CDW
			0.5 * V_OVER_N, // Gamma AFM
			0.5 * V_OVER_N, // Xi AFM
			V_OVER_N, // Occupation Up y
			V_OVER_N, // Occupation Down y
			V_OVER_N, // Gamma SC
			V_OVER_N // Xi SC
		};
	}
	PhaseDataSet HubbardCDW::computePhases(const bool print)
	{
		computeChemicalPotential();

		constexpr double_prec EPSILON = 1e-12;
		double_prec error = 100;
		constexpr int MAX_STEPS = 2500;

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

		PhaseDataSet ret;
		ret.delta_cdw = this->delta_cdw.real();
		ret.delta_afm = this->delta_afm.real();
		ret.delta_sc = this->delta_sc.real();
		ret.gamma_sc = this->gamma_sc.imag();
		ret.xi_sc = this->xi_sc.imag();
		ret.delta_eta = this->delta_eta.imag();

		return ret;
	}
}