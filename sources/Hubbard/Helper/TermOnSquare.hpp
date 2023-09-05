#pragma once
#include "../GlobalDefinitions.hpp"
#include "../Constants.hpp"
#include "../../../../FermionCommute/sources/WickTerm.hpp"
#include "../SquareLattice/UsingBroyden.hpp"
#include "DetailModelConstructor.hpp"

namespace Hubbard::Helper {
	class TermOnSquare : protected DetailModelConstructor<Hubbard::SquareLattice::UsingBroyden>
	{
	protected:
		complex_prec getExpectationValue(const SymbolicOperators::WickOperator& op, const Eigen::Vector2i& momentum_value) const;
		///////////////////////
		// Utility functions //
		///////////////////////
		// Computes the respective x or y component from a given input index
		inline int x(int idx) const {
			return idx / (2 * Constants::K_DISCRETIZATION) - Constants::K_DISCRETIZATION;
		};
		inline int y(int idx) const {
			return idx % (2 * Constants::K_DISCRETIZATION) - Constants::K_DISCRETIZATION;
		};
		inline int equal_up_to_Q(const Eigen::Vector2i& l, const Eigen::Vector2i& r) const {
			if (l == r) return 0;
			if (l(0) == r(0) + Constants::K_DISCRETIZATION || l(0) == r(0) - Constants::K_DISCRETIZATION) {
				if (l(1) == r(1) + Constants::K_DISCRETIZATION || l(1) == r(1) - Constants::K_DISCRETIZATION) {
					return 1;
				}
			}
			return -1;
		};
		// returns a value in [0, N_K), note that N_K = 2*constants::k_disc
		inline void clean_factor_2pi(Eigen::Vector2i& toClean) const {
			// + Q is required for the modulo operation later
			// as well as referencing, which works on indizes from 0 to [2pi] and not from [-pi] to [pi]
			for (int i = 0; i < 2; i++)
			{
				toClean(i) += Constants::K_DISCRETIZATION;
				if (toClean(i) < 0) {
					toClean(i) = ((2 * Constants::K_DISCRETIZATION)
						- abs(toClean(i) % (2 * Constants::K_DISCRETIZATION))) % (2 * Constants::K_DISCRETIZATION);
				}
				else {
					toClean(i) = (toClean(i) % (2 * Constants::K_DISCRETIZATION)) % (2 * Constants::K_DISCRETIZATION);
				}
			}
		};
		
		Eigen::Vector2i computeMomentum(const SymbolicOperators::Momentum& momentum, const std::vector<Eigen::Vector2i>& indizes, const std::vector<char>& momenta) const;

		complex_prec computeTerm(const SymbolicOperators::WickTerm& term, int l, int k) const;

	public:
		TermOnSquare(Utility::InputFileReader& input) : DetailModelConstructor(input) {};
	};
}