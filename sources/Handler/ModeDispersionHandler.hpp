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
		return { i, Hubbard::Constants::K_DISCRETIZATION - 1 };
	}
	static inline Eigen::Vector2i path_R_to_Gamma(int i) {
		return { Hubbard::Constants::K_DISCRETIZATION - 1 - i, Hubbard::Constants::K_DISCRETIZATION - 1 - i };
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
public:
	ModeDispersionHandler(Utility::InputFileReader& input, int _rank, int _numberOfRanks)
		: HandlerBase(input, _rank, _numberOfRanks), ModeHandler(input, _rank, _numberOfRanks) {};
	virtual void execute(Utility::InputFileReader& input) const override;
};