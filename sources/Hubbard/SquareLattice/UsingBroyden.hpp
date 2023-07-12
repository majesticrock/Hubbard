#pragma once
#include "Model2D.hpp"

namespace Hubbard::SquareLattice {
	class UsingBroyden : public Model2D<double>
	{
	private:
		const size_t MaxPreBroydenIterations;
		void init();
	protected:
		virtual void fillHamiltonianHelper(va_list args) override;

		inline void addToParameterSetHelper(const SpinorMatrix& rho, ComplexParameterVector& F, va_list args) override {
			UNPACK_2D;

			const double GAMMA = gamma(k_x, k_y);
			const double XI = xi(k_x, k_y);

			F(0) -= (rho(0, 1) + rho(1, 0) - rho(2, 3) - rho(3, 2)).real(); // CDW
			F(1) -= (rho(0, 1) + rho(1, 0) + rho(2, 3) + rho(3, 2)).real(); // AFM
			F(2) -= (rho(0, 2) + rho(1, 3)); // SC
			F(3) -= GAMMA * (rho(0, 2) - rho(1, 3)); // Gamma SC
			F(4) -= XI * (rho(0, 2) - rho(1, 3)); // Xi SC
			F(5) -= (rho(0, 3) + rho(1, 2)); // Eta
			F(6) -= GAMMA * (rho(0, 0) - rho(1, 1)).real(); // Gamma Occupation Up
			F(7) += GAMMA * (rho(2, 2) - rho(3, 3)).real(); // Gamma Occupation Down
		};
	public:
		explicit UsingBroyden(const ModelParameters& _params, size_t _MaxPreBroydenIterations = 300U);
		UsingBroyden(const ModelParameters& _params, const BaseAttributes& startingValues, size_t _MaxPreBroydenIterations = 300U);

		ModelAttributes<double> computePhases(const PhaseDebuggingPolicy debugPolicy=PhaseDebuggingPolicy{}) override;
	};
}