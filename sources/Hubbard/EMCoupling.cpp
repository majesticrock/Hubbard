#include "EMCoupling.hpp"
#include "Selfconsistency/IterativeSolver.hpp"

namespace Hubbard {
	void EMCoupling::fillHamiltonian()
	{
		hamilton.fill(global_floating_type{});

		NumericalMomentum<Dimension> k; // base momentum
		NumericalMomentum<Dimension> q; // transfer momentum
		size_t i, j;

		do {
			i = k.getIndex();
			hamilton(i, i) = -2. * k.gamma();
			hamilton(i + Constants::BASIS_SIZE, i + Constants::BASIS_SIZE) = 2. * k.gamma();
			do {
				j = (k + q).getIndex();
				hamilton(i, j + Constants::BASIS_SIZE) = this->model_attributes[this->get_sc_index(q)];
				hamilton(j + Constants::BASIS_SIZE, i) = conj(this->model_attributes[this->get_sc_index(q)]);

				if (not q.isZero()) // The q = 0 case is treated on the diagonals
				{
					hamilton(i, j) = this->model_attributes[this->get_cdw_index(-q)];
					hamilton(i + Constants::BASIS_SIZE, j + Constants::BASIS_SIZE) = -conj(this->model_attributes[this->get_cdw_index(-q)]);

					hamilton(j, i) = conj(hamilton(i, j));
					hamilton(j + Constants::BASIS_SIZE, i + Constants::BASIS_SIZE) = conj(hamilton(i + Constants::BASIS_SIZE, j + Constants::BASIS_SIZE));
				}
			} while (q.iterateFullBZ());
			q.reset();
		} while (k.iterateFullBZ());
	}

	void EMCoupling::setParameterSet(ComplexParameterVector& F)
	{
		NumericalMomentum<Dimension> k; // base momentum
		NumericalMomentum<Dimension> q; // transfer momentum
		size_t i, j;

		do {
			i = k.getIndex();
			do {
				j = (k + q).getIndex();
				F(this->get_sc_index(q)) -= this->rho(i, j + Constants::BASIS_SIZE);

				if (not q.isZero()) {
					F(this->get_cdw_index(q)) -= this->rho(i, j) - this->rho(i + Constants::BASIS_SIZE, j + Constants::BASIS_SIZE);
				}
			} while (q.iterateFullBZ());
			q.reset();
		} while (k.iterateFullBZ());
	}

	EMCoupling::EMCoupling(const ModelParameters& _params)
		: BaseModel(_params, static_cast<size_t>(2 * Constants::BASIS_SIZE), static_cast<size_t>(2 * Constants::BASIS_SIZE))
	{
		auto guess = [&]() -> double {
			if (std::abs(this->U) > 1e-12) {
				return std::abs(this->U) * 4. * exp(-2 * 3.1415926 / sqrt(std::abs(this->U)));
			}
			return 0.0;
			};

		this->model_attributes[get_cdw_index(NumericalMomentum<Dimension>())] = guess() / sqrt(2.0);
		this->model_attributes[get_sc_index(NumericalMomentum<Dimension>(0))] = guess() / sqrt(2.0);
		//NumericalMomentum<Dimension> q;
		//do {
		//	auto mod = pow(cos(q.squared_norm() / (2 * BASE_PI)), 2);
		//	this->model_attributes[get_sc_index(q)] = this->U * mod;
		//	this->model_attributes[get_cdw_index(q)] = this->U * mod;
		//} while (q.iterateFullBZ());
	}

	void EMCoupling::iterationStep(const ParameterVector& x, ParameterVector& F)
	{
		F.fill(complex_prec{});
		std::copy(x.begin(), x.end(), this->model_attributes.begin());

		this->fillHamiltonian();
		this->fillRho();
		this->setParameterSet(F);

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

	ModelAttributes<global_floating_type> EMCoupling::computePhases(const PhaseDebuggingPolicy debugPolicy/*=NoWarning*/)
	{
		Selfconsistency::IterativeSolver<complex_prec> solver(this, &model_attributes);

		return ModelAttributes<global_floating_type>(solver.computePhases(debugPolicy), SeperateRealAndImaginary);
	}

	void EMCoupling::computeExpectationValues(std::vector<ValueArray>& expecs, ValueArray& sum_of_all)
	{
		std::cerr << "To be implemented" << std::endl;
	}
}