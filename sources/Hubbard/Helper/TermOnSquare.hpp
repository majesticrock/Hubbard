#pragma once
#include "../GlobalDefinitions.hpp"
#include "../Constants.hpp"
#include "../Models/SquareLattice/UsingBroyden.hpp"
#include "DetailModelConstructor.hpp"
#include "../MomentumIndexUtility.hpp"

namespace Hubbard::Helper {
	class TermOnSquare : protected DetailModelConstructor<Hubbard::Models::SquareLattice::UsingBroyden>
	{
	protected:
		complex_prec getExpectationValue(const SymbolicOperators::WickOperator& op, const Eigen::Vector2i& momentum_value) const;

		Eigen::Vector2i computeMomentum(const SymbolicOperators::MomentumList& momentum, const std::vector<Eigen::Vector2i>& indizes,
			const std::vector<char>& momenta) const;

		complex_prec computeTerm(const SymbolicOperators::WickTerm& term, int l, int k) const;

	public:
		Eigen::Vector2i mode_momentum;

		TermOnSquare(Utility::InputFileReader& input, const Models::ModelParameters& modelParameters, const Eigen::Vector2i& _mode_momentum = { 0, 0 })
			: DetailModelConstructor(input, modelParameters), mode_momentum(_mode_momentum) {};
		TermOnSquare(std::unique_ptr<Hubbard::Models::SquareLattice::UsingBroyden>&& model_ptr)
			: DetailModelConstructor(std::move(model_ptr)) {};
	};
}