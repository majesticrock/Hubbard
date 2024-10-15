#include "EMCoupling.hpp"
#include <Utility/Selfconsistency/IterativeSolver.hpp>

#define DELTA_SC(k) this->model_attributes[this->get_sc_index((k))]
#define DELTA_N(k) this->model_attributes[this->get_n_index((k))]
#define DELTA_AFM(k) this->model_attributes[this->get_afm_index((k))]
#define DELTA_CDW(k) this->model_attributes[this->get_cdw_index((k))]
#define DELTA_ETA(k) this->model_attributes[this->get_eta_index((k))]

namespace Hubbard::Models {
	void EMCoupling::init()
	{
		this->hamilton = SpinorMatrix::Zero(this->SPINOR_SIZE, this->SPINOR_SIZE);
		this->rho = SpinorMatrix::Zero(this->SPINOR_SIZE, this->SPINOR_SIZE);

		this->computeChemicalPotential();
		this->parameterCoefficients = {
			0.5 * this->U_OVER_N - 4. * this->V_OVER_N, // CDW
			0.5 * this->U_OVER_N, // AFM
			this->U_OVER_N, // SC
			this->V_OVER_N, // Gamma SC
			this->V_OVER_N, // Xi SC
			this->U_OVER_N, // Eta
			this->V_OVER_N, // Occupation Up
			this->V_OVER_N, // Occupation Down
			(0.5 * this->U_OVER_N + 4. * this->V_OVER_N) // Phase seperation
		};
	}
	void EMCoupling::computeChemicalPotential()
	{
		this->chemical_potential = 0.5 * this->U + 4 * this->V;
	}
	void EMCoupling::fillHamiltonian(const NumericalMomentum<2>& k)
	{
		hamilton.fill(global_floating_type{});

		hamilton(0, 1) = DELTA_CDW(k) - DELTA_AFM(k);

		hamilton(0, 2) = DELTA_SC(k);
		hamilton(0, 3) = DELTA_ETA(k);

		hamilton(1, 2) = DELTA_ETA(k);
		hamilton(1, 3) = DELTA_SC(k);

		hamilton(2, 3) = -DELTA_CDW(k) - DELTA_AFM(k);

		SpinorMatrix buffer{ hamilton.adjoint() };
		hamilton += buffer;
		global_floating_type eps = -2.0 * k.gamma() + DELTA_N(k);
		hamilton(0, 0) = eps;
		hamilton(1, 1) = -eps;
		hamilton(2, 2) = -eps;
		hamilton(3, 3) = eps;
	}

	void EMCoupling::addToParameterSet(ComplexParameterVector& F, const NumericalMomentum<2>& k)
	{
		NumericalMomentum<2> q{};
		NumericalMomentum<2> q_plus_k, q_minus_k;
		do {
			q_plus_k = q + k;
			q_minus_k = q - k;

			for (int i = 0; i < N_PARAMETERS; ++i) {
				F(q.getIndex() + static_cast<size_t>(i * Constants::BASIS_SIZE)) +=
					(this->phis[q.getIndex()] + parameterCoefficients[i] + 0.5 * V_OVER_N * q.gamma())
					* (this->expectation_values(q_plus_k.getIndex(), i) + this->expectation_values(q_minus_k.getIndex(), i));
			}
		} while (q.iterateFullBZ());
		q = NumericalMomentum<2>::GammaPoint();
	}

	EMCoupling::EMCoupling(const ModelParameters& _params)
		: MomentumBasedModel(_params, 4U, static_cast<size_t>(2 * Constants::BASIS_SIZE))
	{
		auto guess = [&]() -> global_floating_type {
			if (std::abs(this->U) > 1e-12) {
				return std::abs(this->U) * 4. * exp(-2 * 3.1415926 / sqrt(std::abs(this->U)));
			}
			return 0.0;
			};
		// TODO: aendern zu fill
		this->model_attributes[get_cdw_index(NumericalMomentum<2>::GammaPoint())] = guess() / sqrt(2.0);
		this->model_attributes[get_sc_index(NumericalMomentum<2>::GammaPoint())] = guess() / sqrt(2.0);
		//NumericalMomentum<Dimension> q;
		//do {
		//	auto mod = pow(cos(q.squared_norm() / (2 * BASE_PI)), 2);
		//	this->model_attributes[get_sc_index(q)] = this->U * mod;
		//	this->model_attributes[get_cdw_index(q)] = this->U * mod;
		//} while (q.iterateFullBZ());
	}

	void EMCoupling::iterationStep(const ParameterVector& x, ParameterVector& F)
	{
		std::cerr << "TODO!!!" << std::endl;
		F.fill(global_floating_type{});
		std::copy(x.begin(), x.end(), this->model_attributes.begin());

		//this->fillHamiltonian();
		this->fillRho();
		//this->setParameterSet(F);

		for (int i = 0; i < Constants::BASIS_SIZE; ++i) {
			F(i) *= 0.5 * this->U_OVER_N; // SC
			F(i + Constants::BASIS_SIZE) *= 0.5 * this->U_OVER_N; // CDW
		}

		this->setParameters(F);
		F -= x;
	}

	void EMCoupling::getAllEnergies(std::vector<global_floating_type>& reciever) {
		std::cerr << "To be implemented" << std::endl;
	}

	global_floating_type EMCoupling::entropyPerSite() {
		std::cerr << "To be implemented" << std::endl;
		return 0;
	}

	global_floating_type EMCoupling::internalEnergyPerSite() {
		std::cerr << "To be implemented" << std::endl;
		return 0;
	}

	ModelAttributes<global_floating_type> EMCoupling::computePhases()
	{
		auto solver = Utility::Selfconsistency::make_iterative<global_floating_type>(this, &model_attributes);
		return solver.compute();
	}

	void EMCoupling::computeExpectationValues(std::vector<ValueArray>& expecs, ValueArray& sum_of_all)
	{
		std::cerr << "To be implemented" << std::endl;
	}
}