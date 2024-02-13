#pragma once
#include "BaseModel.hpp"
#include "NumericalMomentum.hpp"

namespace Hubbard {
	class EMCoupling : public BaseModel<global_floating_type>
	{
	private:
		/* We structure the selfconsistency parameters as follows:
		*  ---All Delta_SC for different momenta   (N entries)
		*  ---All Delta_CDW for different momenta  (N entries)
		*/
		static constexpr int Dimension = 2;
	protected:
		virtual void fillHamiltonian();
		virtual void addToParameterSet(ComplexParameterVector& F);
	public:
		EMCoupling(const ModelParameters& _params)
			: BaseModel(_params, static_cast<size_t>(2 * Constants::BASIS_SIZE), static_cast<size_t>(2 * Constants::BASIS_SIZE)) 
		{};

		EMCoupling(const ModelParameters& _params, const ModelAttributes<global_floating_type>& startingValues)
			: BaseModel(_params, startingValues) {};

		virtual void iterationStep(const ParameterVector& x, ParameterVector& F) override;

		// saves all one particle energies to reciever
		virtual void getAllEnergies(std::vector<global_floating_type>& reciever) override;

		inline virtual global_floating_type entropyPerSite() override;

		inline virtual global_floating_type internalEnergyPerSite();
	};
}