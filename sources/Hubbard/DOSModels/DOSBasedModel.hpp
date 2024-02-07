#pragma once
#include "../BaseModel.hpp"
#include "../../../../FermionCommute/sources/Coefficient.hpp"
#include "../DensityOfStates/BaseDOS.hpp"
#include "../DensityOfStates/Square.hpp"
#include "../DensityOfStates/SimpleCubic.hpp"
#include "../DensityOfStates/DOSIntegrator.hpp"
#include <algorithm>
#include <mutex>
#include <numeric>

#define DELTA_CDW this->model_attributes[0]
#define DELTA_AFM this->model_attributes[1]
#define DELTA_SC this->model_attributes[2]
#define GAMMA_SC this->model_attributes[3]
// 4 is unused
#define DELTA_ETA this->model_attributes[5]
#define GAMMA_OCCUPATION_UP this->model_attributes[6]
#define GAMMA_OCCUPATION_DOWN this->model_attributes[7]
#define DELTA_PS this->model_attributes[8]

namespace Hubbard {
	template <typename DataType, class DOS>
	class DOSBasedModel : public BaseModel<DataType>
	{
	private:
		static std::mutex dos_mutex;
		constexpr static int NUMBER_OF_PARAMETERS = 9;

		void init() {
			this->hamilton = SpinorMatrix::Zero(Constants::SPINOR_SIZE, Constants::SPINOR_SIZE);
			this->rho = SpinorMatrix::Zero(Constants::SPINOR_SIZE, Constants::SPINOR_SIZE);

			this->model_attributes[4] = 0.;
			// The "1/N"-part is handled in the integration method
			this->V_OVER_N = this->V;
			this->U_OVER_N = this->U;
			this->parameterCoefficients = {
				(0.5 * this->U_OVER_N - DOS::COORDINATION_NUMBER * this->V_OVER_N), // CDW
				0.5 * this->U_OVER_N, // AFM
				this->U_OVER_N, // SC
				this->V_OVER_N, // Gamma SC
				0, // unused
				this->U_OVER_N, // Eta
				this->V_OVER_N, // Occupation Up
				this->V_OVER_N, // Occupation Down
				(0.5 * this->U_OVER_N + DOS::COORDINATION_NUMBER * this->V_OVER_N) // Phase seperation
			};
			this->computeChemicalPotential();

			if (!DOS::computed) {
				std::lock_guard<std::mutex> guard(dos_mutex);
				// Might have been changed by another thread
				if (!DOS::computed) {
					DOS dos;
					dos.computeValues();
					std::cout << "1 - DOS-Norm = " << std::scientific << 1. - DensityOfStates::computeNorm<DOS>() << std::endl;
#ifdef _EXACT_DOS
					//Constants::setBasis(2U * DOS::size());
#else
					Constants::BASIS_SIZE += 1;
#endif
				}
			}
		};

		inline void setParameterSet(ComplexParameterVector& F, const global_floating_type gamma) const {
			F(0) = -(this->rho(0, 1) + this->rho(1, 0) - this->rho(2, 3) - this->rho(3, 2)).real(); // CDW
			F(1) = -(this->rho(0, 1) + this->rho(1, 0) + this->rho(2, 3) + this->rho(3, 2)).real(); // AFM
			F(2) = -(this->rho(0, 2) + this->rho(1, 3)); // SC
			F(3) = -(2.0 / DOS::DIMENSION) * gamma * (this->rho(0, 2) - this->rho(1, 3)); // Gamma SC
			F(4) = 0; // unused
			F(5) = -(this->rho(0, 3) + this->rho(1, 2)); // Eta
			F(6) = -(2.0 / DOS::DIMENSION) * gamma * (this->rho(0, 0) - this->rho(1, 1)).real(); // Gamma Occupation Up
			F(7) = +(2.0 / DOS::DIMENSION) * gamma * (this->rho(2, 2) - this->rho(3, 3)).real(); // Gamma Occupation Down
			//F(8) = -this->rho(0, 4).real() - this->rho(1, 5).real() + this->rho(2, 6).real() + this->rho(3, 7).real(); // PS
		};

	protected:
		using _scalar_integrator = typename DOS::Integrator<global_floating_type>;
		using ParameterVector = typename BaseModel<DataType>::ParameterVector;
		typename DOS::Integrator<ComplexParameterVector> _self_consistency_integrator;

		virtual void computeChemicalPotential() override {
			this->chemical_potential = 0.5 * this->U + DOS::COORDINATION_NUMBER * this->V;
		};

		inline void fillHamiltonian(const global_floating_type& gamma) {
			this->hamilton.fill(global_floating_type{});

			this->hamilton(0, 1) = DELTA_CDW - DELTA_AFM;
			this->hamilton(0, 2) = DELTA_SC + GAMMA_SC * gamma;
			this->hamilton(0, 3) = I * DELTA_ETA;

			this->hamilton(1, 2) = I * DELTA_ETA;
			this->hamilton(1, 3) = DELTA_SC - GAMMA_SC * gamma;
			this->hamilton(2, 3) = -DELTA_CDW - DELTA_AFM;

			SpinorMatrix buffer{ this->hamilton.adjoint() };
			this->hamilton += buffer;
			global_floating_type eps = this->model_attributes.renormalizedEnergy_up(gamma);
			this->hamilton(0, 0) = eps;
			this->hamilton(1, 1) = -eps;
			eps = this->model_attributes.renormalizedEnergy_down(gamma);
			this->hamilton(2, 2) = -eps;
			this->hamilton(3, 3) = eps;
		};

		virtual void iterationStep(const ParameterVector& x, ParameterVector& F) override {
			F.fill(global_floating_type{});
			std::conditional_t<Utility::is_complex<DataType>(), ComplexParameterVector&, ComplexParameterVector> complex_F = F;

			std::copy(x.begin(), x.end(), this->model_attributes.begin());

			auto expectationValues = [this](global_floating_type gamma, ComplexParameterVector& result) {
				this->fillHamiltonian(gamma);
				this->fillRho();
				this->setParameterSet(result, gamma);
				};
	
			complex_F = _self_consistency_integrator.integrate_by_reference_symmetric(expectationValues);

			if constexpr (!std::is_same_v<DataType, complex_prec>) {
				complexParametersToReal(complex_F, F);
			}
			this->applyIteration(F);
			F -= x;
		};
	public:
		DOSBasedModel(const ModelParameters& _params) : BaseModel<DataType>(_params, DOS::DIMENSION),
			_self_consistency_integrator(ComplexParameterVector::Zero(NUMBER_OF_PARAMETERS))
		{
			init();
		};

		template<typename StartingValuesDataType>
		DOSBasedModel(const ModelParameters& _params, const ModelAttributes<StartingValuesDataType>& startingValues)
			: BaseModel<DataType>(_params, startingValues),
			_self_consistency_integrator(ComplexParameterVector::Zero(NUMBER_OF_PARAMETERS))
		{
			init();
		};

		inline virtual global_floating_type entropyPerSite() override {
			using std::log;
			Eigen::SelfAdjointEigenSolver<SpinorMatrix> solver;
			auto procedure = [this, &solver](global_floating_type gamma) {
				this->fillHamiltonian(gamma);
				solver.compute(this->hamilton, false);

				return std::accumulate(solver.eigenvalues().begin(), solver.eigenvalues().end(), global_floating_type{},
					[this](global_floating_type current, global_floating_type toAdd) {
						auto occ = BaseModel<DataType>::fermi_dirac(toAdd);
						// Let's just not take the ln of 0. Negative numbers cannot be reached (because math...)
						return (occ > 1e-12 ? current - occ * log(occ) : current);
					});
				};
			// Devide by two because the matrix representation already includes gamma and -gamma.
			return _scalar_integrator().integrate_by_value(procedure) / 2;
		};

		inline virtual global_floating_type internalEnergyPerSite() override {
			Eigen::SelfAdjointEigenSolver<SpinorMatrix> solver;
			auto procedure = [this, &solver](global_floating_type gamma) {
				this->fillHamiltonian(gamma);
				solver.compute(this->hamilton, false);

				return std::accumulate(solver.eigenvalues().begin(), solver.eigenvalues().end(), global_floating_type{},
					[this](global_floating_type current, global_floating_type toAdd) {
						return current + toAdd * BaseModel<DataType>::fermi_dirac(toAdd);
					});
				};
			// Devide by two because the matrix representation already includes gamma and -gamma.
			return _scalar_integrator().integrate_by_value(procedure) / 2;
		};

		inline global_floating_type computeCoefficient(const SymbolicOperators::Coefficient& coeff, const global_floating_type& gamma) const {
			if (coeff.name == "\\epsilon_0") {
				if (!coeff.dependsOn('k')) throw std::runtime_error("Epsilon should always be k-dependent.");
				if (coeff.momentum.momentum_list.size() > 1U) throw std::runtime_error("Epsilon depends on more than k.");
				return ((coeff.momentum.add_Q ? 2 * gamma : -2 * gamma) - this->chemical_potential);
			}
			if (coeff.name == "\\frac{U}{N}") {
				return this->U_OVER_N;
			}
			if (coeff.name == "\\tilde{V}") {
				if (coeff.momentum.add_Q) throw std::runtime_error("V includes Q, this should not occur.");
				if (coeff.dependsOnMomentum()) {
					return  this->V_OVER_N * (coeff.dependsOnTwoMomenta() ? (gamma / DOS::DIMENSION) : gamma);
				}
				else {
					return DOS::DIMENSION * this->V_OVER_N;
				}
			}
			throw(std::invalid_argument("Could not find the coefficient: " + coeff.name));
		};

		inline global_floating_type getGammaFromIndex(int gamma_idx) const {
#ifdef _EXACT_DOS
			return (Constants::HALF_BASIS - gamma_idx - 0.5) * DOS::LOWER_BORDER / Constants::BASIS_SIZE;
			//return DOS::abscissa_v(gamma_idx);
#else
			return 2 * (Constants::HALF_BASIS - gamma_idx) * DOS::LOWER_BORDER / (Constants::BASIS_SIZE - 1);
#endif
		};

		inline int shiftByQ(int k) const {
#ifdef _EXACT_DOS
			return (-k + 2 * Constants::HALF_BASIS - 1);
			//return k + (k < Constants::HALF_BASIS ? Constants::HALF_BASIS : -Constants::HALF_BASIS);
#else
			return (-k + 2 * Constants::HALF_BASIS);
#endif
		};

		virtual void computeExpectationValues(std::vector<ValueArray>& expecs, ValueArray& sum_of_all) override {
			expecs = std::vector<ValueArray>(8U, ValueArray::Zero(Constants::BASIS_SIZE, 1));
			sum_of_all = ValueArray::Zero(8, 2);

			for (int i = 0; i < Constants::BASIS_SIZE; ++i)
			{
				const global_floating_type gamma = this->getGammaFromIndex(i);

				this->fillHamiltonian(gamma);
				this->fillRho();
				// n_up
				expecs[0](i, 0) = this->get_n_up();
				// g_up
				expecs[1](i, 0) = this->get_g_up();
				// f
				expecs[2](i, 0) = this->get_f();
				// eta
				expecs[3](i, 0) = this->get_eta();
				// n_down
				expecs[4](i, 0) = this->get_n_down();
				// g_down
				expecs[5](i, 0) = this->get_g_down();
				// n_up + n_down
				expecs[6](i, 0) = this->get_n_up_plus_down();
				// g_up + g_down
				expecs[7](i, 0) = this->get_g_up_plus_down();
				if (abs(this->rho(3, 0)) > 1e-10) {
					std::cerr << "Warning: <eta> does not vanish! " << this->rho(3, 0) << std::endl;
				}
			}

			auto func_sum = [this](const global_floating_type& gamma) {
				this->fillHamiltonian(gamma);
				this->fillRho();
				ValueArray ret = ValueArray::Zero(8, 2);
				ret(0, 0) = this->get_n_up();
				ret(1, 0) = this->get_g_up();
				ret(2, 0) = this->get_f();
				ret(3, 0) = this->get_eta();
				ret(4, 0) = this->get_n_down();
				ret(5, 0) = this->get_g_down();
				ret(6, 0) = this->get_n_up_plus_down();
				ret(7, 0) = this->get_g_up_plus_down();

				ret(0, 1) = gamma * this->get_n_up();
				ret(1, 1) = gamma * this->get_g_up();
				ret(2, 1) = gamma * this->get_f();
				ret(3, 1) = gamma * this->get_eta();
				ret(4, 1) = gamma * this->get_n_down();
				ret(5, 1) = gamma * this->get_g_down();
				ret(6, 1) = gamma * this->get_n_up_plus_down();
				ret(7, 1) = gamma * this->get_g_up_plus_down();

				return ret;
				};

			typename DOS::Integrator<ValueArray> integrator;
			sum_of_all = integrator.integrate_by_value(func_sum);
		};

		// saves all one particle energies to reciever
		virtual void getAllEnergies(std::vector<global_floating_type>& reciever) override
		{
			reciever.reserve(Constants::SPINOR_SIZE * (Constants::BASIS_SIZE + 1));
			Eigen::SelfAdjointEigenSolver<SpinorMatrix> solver;
			const global_floating_type step{ DOS::LOWER_BORDER / Constants::BASIS_SIZE };

			for (int g = 0; g <= Constants::BASIS_SIZE; ++g)
			{
				this->fillHamiltonian(g * step);
				solver.compute(this->hamilton, false);

				for (int i = 0; i < Constants::SPINOR_SIZE; ++i)
				{
					reciever.push_back(solver.eigenvalues()(i));
				}
			}

			std::sort(reciever.begin(), reciever.end());
		};

		global_floating_type higgs_sum_rule() {
			auto compute_values = [this](global_floating_type gamma) {
				Eigen::Vector4d values = Eigen::Vector4d::Zero();

				this->fillHamiltonian(gamma);
				this->fillRho();

				values(0) = this->get_f().real(); // f
				values(1) = this->get_n_down(); // n
				values(2) = this->get_n_down() * this->get_n_down(); // n^2
				values(3) = this->get_f().real() * this->get_f().real(); // f^2

				return values;
				};

			typename DOS::Integrator<Eigen::Vector4d> integ;
			Eigen::Vector4d res = integ.integrate_by_value(compute_values);
			res(0) *= 4 * res(0);

			std::cout << res(0) << "   " << res(1) << "   " << res(2) << "   " << res(3) << std::endl;
			return 1 - 2 * (res(1) - res(2) + res(3));
		};

		global_floating_type cdw_in_sc_sum_rule() {
			auto compute_values = [this](global_floating_type gamma) {
				Eigen::Vector2d values = Eigen::Vector2d::Zero();

				this->fillHamiltonian(gamma);
				this->fillRho();

				values(0) = this->get_f().real() * this->get_f_Q().real(); // f_k * f_k+q
				values(1) = this->get_n_down() * this->get_n_down_Q(); // n_k * n_k+q

				return values;
				};

			typename DOS::Integrator<Eigen::Vector2d> integ;
			Eigen::Vector2d res = integ.integrate_by_value(compute_values);

			std::cout << res(0) << "   " << res(1) << std::endl;
			return 1 + 2 * (res(0) - res(1));
		};
	};

	template <typename DataType, class DOS>
	std::mutex DOSBasedModel<DataType, DOS>::dos_mutex;
}