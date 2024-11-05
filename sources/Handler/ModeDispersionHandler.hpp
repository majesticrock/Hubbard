#pragma once

#include "ModeHandler.hpp"
#include "../Hubbard/Constants.hpp"
#include <Eigen/Dense>

class ModeDispersionHandler : public ModeHandler {
private:
	static inline Eigen::Vector2i path_Gamma_to_X(int i) {
		return { 0, i };
	}
	static inline Eigen::Vector2i path_X_to_R(int i) {
		return { i, Hubbard::Constants::K_DISCRETIZATION };
	}
	static inline Eigen::Vector2i path_R_to_Gamma(int i) {
		return { Hubbard::Constants::K_DISCRETIZATION - i, Hubbard::Constants::K_DISCRETIZATION - i };
	}
	static inline Eigen::Vector2i eval_point(int i) {
		switch (i / Hubbard::Constants::K_DISCRETIZATION) {
		case 0:
			return path_Gamma_to_X(i);
			break;
		case 1:
			return path_X_to_R(i - Hubbard::Constants::K_DISCRETIZATION);
			break;
		case 2:
			return path_R_to_Gamma(i - 2 * Hubbard::Constants::K_DISCRETIZATION);
			break;
		default:
			throw std::runtime_error("Could not find eval point!");
		}
	}

	const int eval_index{};
	inline std::string name_prefix() const { return (eval_index < 0 ? "" : (std::to_string(eval_index) + "_"));}
public:
	ModeDispersionHandler(Utility::InputFileReader& input, int _rank, int _numberOfRanks, int _eval_index)
		: HandlerBase(input, _rank, _numberOfRanks), ModeHandler(input, _rank, _numberOfRanks), eval_index{_eval_index} {};
	virtual void execute(Utility::InputFileReader& input) const override;
};