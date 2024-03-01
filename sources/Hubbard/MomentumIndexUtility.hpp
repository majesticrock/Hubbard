#include "Constants.hpp"

///////////////////////
// Utility functions //
///////////////////////

namespace Hubbard {
	// Computes the respective x or y component from a given input index
	inline int x(int idx) {
		return idx / (2 * Constants::K_DISCRETIZATION) - Constants::K_DISCRETIZATION;
	};
	inline int y(int idx) {
		return idx % (2 * Constants::K_DISCRETIZATION) - Constants::K_DISCRETIZATION;
	};
	inline global_floating_type gammaFromIndex(int index) {
		return cos(Constants::PI_DIV_DISCRETIZATION * x(index)) + cos(Constants::PI_DIV_DISCRETIZATION * y(index));
	};
	inline int equal_up_to_Q(const Eigen::Vector2i& l, const Eigen::Vector2i& r) {
		if (l == r) return 0;
		if (l(0) == r(0) + Constants::K_DISCRETIZATION || l(0) == r(0) - Constants::K_DISCRETIZATION) {
			if (l(1) == r(1) + Constants::K_DISCRETIZATION || l(1) == r(1) - Constants::K_DISCRETIZATION) {
				return 1;
			}
		}
		return -1;
	};
	// returns a value in [0, N_K), note that N_K = 2*constants::k_disc
	inline void clean_factor_2pi(Eigen::Vector2i& toClean) {
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
	inline int addQTo(int k) {
		Eigen::Vector2i l_buf_vec = { x(k), y(k) };
		l_buf_vec(0) += Constants::K_DISCRETIZATION;
		l_buf_vec(1) += Constants::K_DISCRETIZATION;
		clean_factor_2pi(l_buf_vec);
		return l_buf_vec(0) * 2 * Constants::K_DISCRETIZATION + l_buf_vec(1);
	};
}