#pragma once
#include "BaseModel.hpp"
#include "NumericalMomentum.hpp"

namespace Hubbard {
	class EMCoupling : public BaseModel<complex_prec>
	{
	private:
		/* We structure the selfconsistency parameters as follows:
		*  ---All Delta_SC for different momenta   (N entries)
		*  ---All Delta_CDW for different momenta  (N entries)
		*/
		static constexpr int Dimension = 2;
		
		inline size_t get_sc_index(const NumericalMomentum<Dimension>& q) const {
			return q.getIndex();
		};
		inline size_t get_cdw_index(const NumericalMomentum<Dimension>& q) const {
			return q.getIndex() + static_cast<size_t>(Constants::BASIS_SIZE);
		};
	protected:
		virtual void fillHamiltonian();
		virtual void setParameterSet(ComplexParameterVector& F);
	public:
		EMCoupling(const ModelParameters& _params)
			: BaseModel(_params, static_cast<size_t>(2 * Constants::BASIS_SIZE), static_cast<size_t>(2 * Constants::BASIS_SIZE)) 
		{};

		EMCoupling(const ModelParameters& _params, const ModelAttributes<complex_prec>& startingValues)
			: BaseModel(_params, startingValues) {};

		virtual void iterationStep(const ParameterVector& x, ParameterVector& F) override;

		// saves all one particle energies to reciever
		virtual void getAllEnergies(std::vector<global_floating_type>& reciever) override;

		virtual global_floating_type entropyPerSite() override;

		virtual global_floating_type internalEnergyPerSite() override;

		virtual ModelAttributes<global_floating_type> computePhases(const PhaseDebuggingPolicy debugPolicy = NoWarning) override;

		virtual void computeExpectationValues(std::vector<ValueArray>& expecs, ValueArray& sum_of_all) override;
	};
}