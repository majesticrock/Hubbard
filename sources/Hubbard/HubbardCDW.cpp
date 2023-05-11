#define _USE_MATH_DEFINES

#include "HubbardCDW.hpp"

namespace Hubbard {
	void HubbardCDW::computeChemicalPotential()
	{
		Model::computeChemicalPotential();
		chemical_potential += 0.5 * V;
	}
	void HubbardCDW::fillHamiltonian(double_prec k_x, double_prec k_y)
	{
		complex_h.fill(0);

		complex_h(0, 1) = delta_cdw;
		complex_h(0, 2) = delta_sc - I * (2 * xi_sc * (cos(k_x) + cos(k_y)));
		complex_h(0, 3) = I * (delta_eta - 2 * xi_sc * (cos(k_x) + cos(k_y)));

		complex_h(1, 2) = I * (delta_eta + 2 * xi_sc * (cos(k_x) + cos(k_y)));
		complex_h(1, 3) = delta_sc + I * (2 * xi_sc * (cos(k_x) + cos(k_y)));
		complex_h(2, 3) = -delta_cdw;

		Matrix_4cL buffer = complex_h.adjoint();
		complex_h += buffer;
		const double_prec eps = unperturbed_energy(k_x, k_y) - (2 * (cos(k_x) + cos(k_y)) * delta_occupation);
		complex_h(0, 0) = eps;
		complex_h(1, 1) = -eps;
		complex_h(2, 2) = -eps;
		complex_h(3, 3) = eps;
	}
	HubbardCDW::HubbardCDW(ModelParameters& _params, int _number_of_basis_terms, int _start_basis_at)
		: Model(_params, _number_of_basis_terms, _start_basis_at), V(_params.V)
	{
		this->delta_cdw = std::abs(U - V) * 0.5;
		this->delta_sc = std::abs(U + V) * 0.5;
		if (V > 0) {
			this->delta_sc *= 0.25;
		}
		else if (V < 0) {
			this->delta_cdw *= 0.25;
		}
		this->delta_occupation = V * 0.1;
		this->xi_sc = V * 0.1;
		this->xi_eta = V * 0.1;

		this->hamilton = Matrix_L::Zero(4, 4);
	}
	Model::data_set HubbardCDW::computePhases(const bool print)
	{
		Matrix_4cL rho = Matrix_4cL::Zero(4, 4);
		Eigen::SelfAdjointEigenSolver<Matrix_4cL> solver;
		constexpr double_prec EPSILON = 1e-8;
		double_prec sc = 0, cdw = 0, eta = 0, cos_n = 0;
		double_prec old_parameters[4] = { 100, 100, 100, 100 };
		double_prec error = 100;
		double_prec error_cdw = 100, error_cdw_osc = 100;
		double_prec error_sc = 100, error_sc_osc = 100;
		double_prec error_eta = 100, error_eta_osc = 100;
		double_prec error_cos_n = 100;
		constexpr int MAX_STEPS = 1000;
		for (size_t i = 0; i < MAX_STEPS && error > EPSILON; i++)
		{
			sc = 0;
			cdw = 0;
			eta = 0;
			cos_n = 0;

			complex_prec c_cdw = { 0, 0 }, c_sc = { 0, 0 }, c_eta = { 0, 0 };
			complex_prec c_xi_sc = { 0,0 }, c_xi_eta = { 0,0 };

			for (int k = -Constants::K_DISCRETIZATION; k < Constants::K_DISCRETIZATION; k++)
			{
				double_prec k_x = (k * M_PI) / Constants::K_DISCRETIZATION;
				for (int l = -Constants::K_DISCRETIZATION; l < Constants::K_DISCRETIZATION; l++)
				{
					fillHamiltonian(k_x, (l * M_PI) / Constants::K_DISCRETIZATION);
					solver.compute(complex_h);

					rho.fill(0);
					for (int i = 0; i < 4; i++)
					{
						rho(i, i) = fermi_dirac(solver.eigenvalues()(i));
					}
					rho = solver.eigenvectors() * rho * solver.eigenvectors().adjoint();

					c_cdw += 0.5 * (rho(0, 1) - rho(2, 3));
					c_sc += rho(0, 2);
					c_eta += rho(0, 3);
					cos_n += 0.5 * cos(k_x) * (rho(0, 0).real() + 1 - rho(2, 2).real());

					c_xi_sc += cos(k_x) * rho(0, 2);
					c_xi_eta += cos(k_x) * rho(0, 3);
				}
			}

			if (std::abs(c_cdw.imag()) > 1e-8) {
				std::cerr << "cdw: " << c_cdw << std::endl;
			}
			if (std::abs(c_sc.imag()) > 1e-8) {
				std::cerr << "sc: " << c_sc << std::endl;
			}
			if (std::abs(c_eta.real()) > 1e-8) {
				std::cerr << "eta: " << c_eta << std::endl;
			}
			if (std::abs(c_xi_sc.real()) > 1e-8) {
				std::cerr << "xi sc: " << c_eta << std::endl;
			}
			if (std::abs(c_xi_eta.real()) > 1e-8) {
				std::cerr << "xi eta: " << c_eta << std::endl;
			}

			cdw = c_cdw.real();
			sc  = c_sc.real();
			eta = c_eta.imag();

			old_parameters[0] = delta_cdw;
			old_parameters[1] = delta_sc;
			old_parameters[2] = delta_eta;
			old_parameters[3] = delta_occupation;

			setParameters(cdw, sc, eta, cos_n, c_xi_sc.imag(), c_xi_eta.imag());

			error_cdw = std::abs(delta_cdw - old_parameters[0]);
			error_sc = std::abs(delta_sc - old_parameters[1]);
			error_eta = std::abs(delta_eta - old_parameters[2]);
			error_cos_n = std::abs(delta_occupation - old_parameters[3]);
			error = error_cdw + error_sc + error_eta + error_cos_n;

			delta_cdw = 0.5 * old_parameters[0] + 0.5 * delta_cdw;
			delta_sc = 0.5 * old_parameters[1] + 0.5 * delta_sc;
			delta_eta = 0.5 * old_parameters[2] + 0.5 * delta_eta;
			delta_occupation = 0.5 * old_parameters[3] + 0.5 * delta_occupation;

			error_cdw_osc = std::abs(delta_cdw + old_parameters[0]);
			error_sc_osc = std::abs(delta_sc + old_parameters[1]);
			error_eta_osc = std::abs(delta_eta + old_parameters[2]);

			delta_cdw = ((error_cdw_osc) < EPSILON) ? 0 : delta_cdw;
			delta_sc = ((error_sc_osc) < EPSILON) ? 0 : delta_sc;
			delta_eta = ((error_eta_osc) < EPSILON) ? 0 : delta_eta;

			if (print) {
				double_prec total = 0;
				for (size_t i = 0; i < 3; i++)
				{
					total += old_parameters[i] * old_parameters[i];
				}
				std::cout << i << ":\t" << std::fixed << std::setprecision(8)
					<< delta_cdw << "\t" << delta_sc << "\t" << delta_eta << "\t" << delta_occupation 
					<< "\t" << xi_sc << "\t" << xi_eta << std::endl;
			}
			if (i == MAX_STEPS - 1) {
				std::cerr << "[T, U] = [" << this->temperature << ", " << this->U << "]\tConvergence at " << error << std::endl;
				delta_cdw = 0;
				delta_sc = 0;
				delta_eta = 0;
			}
		}
		data_set ret;
		ret.delta_cdw = delta_cdw;
		ret.delta_sc = delta_sc;
		ret.delta_eta = delta_eta;

		return ret;
	}
}