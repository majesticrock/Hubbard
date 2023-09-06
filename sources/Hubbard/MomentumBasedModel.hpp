#pragma once
#include "BaseModel.hpp"
#include "../../../FermionCommute/sources/Coefficient.hpp"
#include "NumericalMomentum.hpp"

namespace Hubbard {
	template <typename DataType, size_t Dimension>
	class MomentumBasedModel : public BaseModel<DataType>
	{
	protected:
		using ParameterVector = typename BaseModel<DataType>::ParameterVector;
		virtual void fillHamiltonian(const NumericalMomentum<Dimension>& k_values) = 0;
		virtual void addToParameterSet(ComplexParameterVector& F, const NumericalMomentum<Dimension>& k_values) = 0;
	public:
		virtual void iterationStep(const ParameterVector& x, ParameterVector& F) override {
			F.fill(global_floating_type{});
			std::conditional_t<Utility::is_complex<DataType>(),
				ComplexParameterVector&, ComplexParameterVector> complex_F = F;

			std::copy(x.begin(), x.end(), this->model_attributes.begin());

			NumericalMomentum<Dimension> ks;
			while (ks.iterateHalfBZ()) {
				this->fillHamiltonian(ks);
				this->fillRho();
				this->addToParameterSet(complex_F, ks);
			}

			if constexpr (!Utility::is_complex<DataType>()) {
				complexParametersToReal(complex_F, F);
			}
			this->applyIteration(F);

			F -= x;
		};
		
		MomentumBasedModel(const ModelParameters& _params)
			: BaseModel<DataType>(_params) {};

		template<typename StartingValuesDataType>
		MomentumBasedModel(const ModelParameters& _params, const ModelAttributes<StartingValuesDataType>& startingValues)
			: BaseModel<DataType>(_params, startingValues) {};

		inline global_floating_type computeCoefficient(const SymbolicOperators::Coefficient& coeff, const Eigen::Vector<int, Dimension>& momentum) const {
			if (coeff.name == "\\epsilon_0") {
				NumericalMomentum<Dimension> temp{ (momentum.array() - Constants::K_DISCRETIZATION).eval() };
				return temp.unperturbed_energy() - this->chemical_potential;
			}
			if (coeff.name == "\\frac{U}{N}") {
				return this->U_OVER_N;
			}
			if (coeff.name == "\\tilde{V}") {
				NumericalMomentum<Dimension> temp{ (momentum.array() - Constants::K_DISCRETIZATION).eval() };
				return this->V_OVER_N * temp.gamma();
			}
			throw(std::invalid_argument("Could not find the coefficient: " + coeff.name));
		};

		// saves all one particle energies to reciever
		virtual void getAllEnergies(std::vector<global_floating_type>& reciever) override
		{
			reciever.reserve(2 * Constants::BASIS_SIZE);
			Eigen::SelfAdjointEigenSolver<SpinorMatrix> solver;

			NumericalMomentum<Dimension> ks;
			while (ks.iterateHalfBZ()) {
				this->fillHamiltonian(ks);
				solver.compute(this->hamilton, false);
				for (int i = 0; i < 4; i++)
				{
					reciever.push_back(solver.eigenvalues()(i));
				}
			}
		};

		inline virtual global_floating_type entropyPerSite() override {
			using std::log;
			global_floating_type entropy{};

			NumericalMomentum<Dimension> ks;
			while (ks.iterateHalfBZ()) {
				this->fillHamiltonian(ks);
				this->hamilton_solver.compute(this->hamilton, false);

				entropy += std::accumulate(this->hamilton_solver.eigenvalues().begin(), this->hamilton_solver.eigenvalues().end(), global_floating_type{},
					[this](global_floating_type current, global_floating_type toAdd) {
						auto occ = BaseModel<DataType>::fermi_dirac(toAdd);
						// Let's just not take the ln of 0. Negative numbers cannot be reached (because math...)
						return (occ > 1e-12 ? current - occ * log(occ) : current);
					});
			}
			return entropy / Constants::BASIS_SIZE;
		};

		inline virtual global_floating_type internalEnergyPerSite() override {
			global_floating_type energy{};

			NumericalMomentum<Dimension> ks;
			while (ks.iterateHalfBZ()) {
				this->fillHamiltonian(ks);
				this->hamilton_solver.compute(this->hamilton, false);

				energy += std::accumulate(this->hamilton_solver.eigenvalues().begin(), this->hamilton_solver.eigenvalues().end(), global_floating_type{},
					[this](global_floating_type current, global_floating_type toAdd) {
						return current + toAdd * BaseModel<DataType>::fermi_dirac(toAdd);
					});
			}
			return energy / Constants::BASIS_SIZE;
		};
	};
}