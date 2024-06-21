#pragma once
#include "MomentumBasedModel.hpp"

namespace Hubbard {
	class EMCoupling : public MomentumBasedModel<global_floating_type, 2>
	{
	private:
		std::vector<global_floating_type> phis;
		ValueArray expectation_values;

		static constexpr int N_PARAMETERS = 5;
		/* We structure the selfconsistency parameters as follows:
		*  ---All Delta_SC for different momenta   (N entries)
		*  ---All Delta_n for different momenta  (N entries)
		*  ---All Delta_CDW for different momenta  (N entries)
		*  ---All Delta_AFM for different momenta  (N entries)
		*  ---All Delta_eta for different momenta  (N entries)
		*/
		void init() override;
		void computeChemicalPotential() override;

		inline size_t get_sc_index(const NumericalMomentum<2>& q) const {
			return q.getIndex();
		};
		inline size_t get_n_index(const NumericalMomentum<2>& q) const {
			return q.getIndex() + static_cast<size_t>(Constants::BASIS_SIZE);
		};
		inline size_t get_cdw_index(const NumericalMomentum<2>& q) const {
			return q.getIndex() + static_cast<size_t>(2 * Constants::BASIS_SIZE);
		};
		inline size_t get_afm_index(const NumericalMomentum<2>& q) const {
			return q.getIndex() + static_cast<size_t>(3 * Constants::BASIS_SIZE);
		};
		inline size_t get_eta_index(const NumericalMomentum<2>& q) const {
			return q.getIndex() + static_cast<size_t>(4 * Constants::BASIS_SIZE);
		};

	protected:
		virtual void fillHamiltonian(const NumericalMomentum<2>& k) override;
		void addToParameterSet(ComplexParameterVector& F, const NumericalMomentum<2>& k) override;
	public:
		EMCoupling(const ModelParameters& _params);

		EMCoupling(const ModelParameters& _params, const ModelAttributes<global_floating_type>& startingValues)
			: MomentumBasedModel(_params, startingValues) {};

		virtual void iterationStep(const ParameterVector& x, ParameterVector& F) override;

		// saves all one particle energies to reciever
		virtual void getAllEnergies(std::vector<global_floating_type>& reciever) override;

		virtual global_floating_type entropyPerSite() override;

		virtual global_floating_type internalEnergyPerSite() override;

		virtual ModelAttributes<global_floating_type> computePhases() override;

		virtual void computeExpectationValues(std::vector<ValueArray>& expecs, ValueArray& sum_of_all) override;
	};
}