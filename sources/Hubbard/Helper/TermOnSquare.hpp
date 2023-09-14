#pragma once
#include "../GlobalDefinitions.hpp"
#include "../Constants.hpp"
#include "../SquareLattice/UsingBroyden.hpp"
#include "DetailModelConstructor.hpp"
#include "../MomentumIndexUtility.hpp"

namespace Hubbard::Helper {
	class TermOnSquare : protected DetailModelConstructor<Hubbard::SquareLattice::UsingBroyden>
	{
	protected:
		complex_prec getExpectationValue(const SymbolicOperators::WickOperator& op, const Eigen::Vector2i& momentum_value) const;

		Eigen::Vector2i computeMomentum(const SymbolicOperators::Momentum& momentum, const std::vector<Eigen::Vector2i>& indizes, const std::vector<char>& momenta) const;

		complex_prec computeTerm(const SymbolicOperators::WickTerm& term, int l, int k) const;

	public:
		TermOnSquare(Utility::InputFileReader& input) : DetailModelConstructor(input) {};
	};
}