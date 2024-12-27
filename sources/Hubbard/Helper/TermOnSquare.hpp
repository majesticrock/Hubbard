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
		global_floating_type getExpectationValue(const mrock::SymbolicOperators::WickOperator& op, const Eigen::Vector2i& momentum_value) const;

		Eigen::Vector2i compute_momentum_list(const mrock::SymbolicOperators::MomentumList& momentum, const Eigen::Vector2i& k, const Eigen::Vector2i& l) const;
		Eigen::Vector2i compute_momentum_no_q(const mrock::SymbolicOperators::Momentum& momentum, const Eigen::Vector2i& k, const Eigen::Vector2i& l) const;

		global_floating_type compute_sum(const mrock::SymbolicOperators::WickTerm& term, const Eigen::Vector2i& k, const Eigen::Vector2i& l) const;

		global_floating_type computeTerm(const mrock::SymbolicOperators::WickTerm& term, const int l, const int k) const;

	public:
		Eigen::Vector2i mode_momentum;

		TermOnSquare(mrock::Utility::InputFileReader& input, const Models::ModelParameters& modelParameters, const Eigen::Vector2i& _mode_momentum = { 0, 0 })
			: DetailModelConstructor(input, modelParameters), mode_momentum(_mode_momentum) {};
		TermOnSquare(std::unique_ptr<Hubbard::Models::SquareLattice::UsingBroyden>&& model_ptr)
			: DetailModelConstructor(std::move(model_ptr)) {};
	};
}