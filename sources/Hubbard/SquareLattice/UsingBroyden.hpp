#pragma once
#include "Model.hpp"
#include "../BaseModelAttributes.hpp"

namespace Hubbard::SquareLattice {
	class UsingBroyden : public Model<double_prec>, public BaseModelRealAttributes
	{
	protected:
		virtual void fillHamiltonianHelper(va_list args) override;

		inline void addToParameterSetHelper(const SpinorMatrix& rho, ParameterVector& F, va_list args) override {
			UNPACK_2D;
			fillHamiltonian(k_x, k_y);
		};
	public:
		UsingBroyden(const ModelParameters& _params);

		PhaseDataSet computePhases(const bool print = false) override;
	};
}