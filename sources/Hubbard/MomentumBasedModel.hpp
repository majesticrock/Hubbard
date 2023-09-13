#pragma once
#include "BaseModel.hpp"
#include "../../../FermionCommute/sources/Coefficient.hpp"
#include "NumericalMomentum.hpp"

namespace Hubbard {
	inline int addQTo(int k){
		Eigen::Vector2i l_buf_vec = { x(k), y(k) };
		l_buf_vec(0) += Constants::K_DISCRETIZATION;
		l_buf_vec(1) += Constants::K_DISCRETIZATION;
		clean_factor_2pi(l_buf_vec);
		return l_buf_vec(0) * 2 * Constants::K_DISCRETIZATION + l_buf_vec(1);
	};

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
			do {
				this->fillHamiltonian(ks);
				this->fillRho();
				this->addToParameterSet(complex_F, ks);
			} while (ks.iterateHalfBZ());

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
				NumericalMomentum<Dimension> temp{ index_vector_to_k_vector(momentum) };
				return temp.unperturbed_energy() - this->chemical_potential;
			}
			if (coeff.name == "\\frac{U}{N}") {
				return this->U_OVER_N;
			}
			if (coeff.name == "\\tilde{V}") {
				NumericalMomentum<Dimension> temp{ index_vector_to_k_vector(momentum) };
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
			do {
				this->fillHamiltonian(ks);
				solver.compute(this->hamilton, false);
				for (int i = 0; i < 4; i++)
				{
					reciever.push_back(solver.eigenvalues()(i));
				}
			} while (ks.iterateHalfBZ());
		};

		inline virtual global_floating_type entropyPerSite() override {
			using std::log;
			global_floating_type entropy{};

			NumericalMomentum<Dimension> ks;
			do {
				this->fillHamiltonian(ks);
				this->hamilton_solver.compute(this->hamilton, false);

				entropy += std::accumulate(this->hamilton_solver.eigenvalues().begin(), this->hamilton_solver.eigenvalues().end(), global_floating_type{},
					[this](global_floating_type current, global_floating_type toAdd) {
						auto occ = BaseModel<DataType>::fermi_dirac(toAdd);
						// Let's just not take the ln of 0. Negative numbers cannot be reached (because math...)
						return (occ > 1e-12 ? current - occ * log(occ) : current);
					});
			} while (ks.iterateHalfBZ());
			return entropy / Constants::BASIS_SIZE;
		};

		inline virtual global_floating_type internalEnergyPerSite() override {
			global_floating_type energy{};

			NumericalMomentum<Dimension> ks;
			do {
				this->fillHamiltonian(ks);
				this->hamilton_solver.compute(this->hamilton, false);

				energy += std::accumulate(this->hamilton_solver.eigenvalues().begin(), this->hamilton_solver.eigenvalues().end(), global_floating_type{},
					[this](global_floating_type current, global_floating_type toAdd) {
						return current + toAdd * BaseModel<DataType>::fermi_dirac(toAdd);
					});
			} while (ks.iterateHalfBZ());
			return energy / Constants::BASIS_SIZE;
		};
	};
}