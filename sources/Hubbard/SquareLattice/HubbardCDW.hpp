#pragma once
#include "Model2D.hpp"

namespace Hubbard::SquareLattice {
	class HubbardCDW : public Model2D<complex_prec>
	{
	private:
		void init();
	protected:
		virtual void fillHamiltonianHelper(va_list args) override;

		inline void addToParameterSetHelper(const SpinorMatrix& rho, ParameterVector& F, va_list args) override {
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
			F(8) += I * GAMMA * (rho(0, 1) - rho(1, 0) + rho(2, 3) - rho(3, 2)).imag(); // Gamma CDW
			F(9) += I * XI * (rho(0, 1) - rho(1, 0) + rho(2, 3) - rho(3, 2)).imag(); // Xi CDW
			F(10) += I * GAMMA * (rho(0, 1) - rho(1, 0) - rho(2, 3) + rho(3, 2)).imag(); // Gamma AFM
			F(11) += I * XI * (rho(0, 1) - rho(1, 0) - rho(2, 3) + rho(3, 2)).imag(); // Xi AFM
			F(12) -= XI * (rho(0, 0) - rho(1, 1)).real(); // Xi Occupation Up
			F(13) += XI * (rho(2, 2) - rho(3, 3)).real(); // Xi Occupation Down
			F(14) -= GAMMA * (rho(0, 3) - rho(1, 2)); // Gamma eta
			F(15) -= XI * (rho(0, 3) - rho(1, 2)); // Xi eta
		};
	public:
		HubbardCDW(const ModelParameters& _params);

		template<typename StartingValuesDataType>
		HubbardCDW(const ModelParameters& _params, const ModelAttributes< StartingValuesDataType >& startingValues)
			: Model2D(_params, startingValues) {};

		virtual ModelAttributes<double> computePhases(const PhaseDebuggingPolicy debugPolicy=PhaseDebuggingPolicy{}) override;
	};
}