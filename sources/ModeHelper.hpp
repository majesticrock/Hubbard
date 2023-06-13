#pragma once
#include <memory>

#include "Hubbard/Model.hpp"

class ModeHelper {
protected:
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Matrix_L;
	typedef Eigen::Vector<double, Eigen::Dynamic> Vector_L;
	typedef Utility::Resolvent<double> Resolvent_L;

    std::unique_ptr<Hubbard::Model> model;
    /*
	* 0 - n
	* 1 - g_up
	* 2 - sc
	* 3 - eta
    * 4 - g_down
	*/
	std::vector<Hubbard::Matrix_L> expecs;
	double sum_of_all[5] = { 0, 0, 0, 0, 0 };

	Hubbard::Matrix_L M, N;
	Hubbard::Matrix_L K_plus, K_minus, L;
	int number_of_basis_terms;
	int start_basis_at;

	const std::map<std::pair<std::string, bool>, int> wick_map = { {"n", 0}, {"g", 1}, {"f", 2}, {"\\eta", 3} };
	std::vector<std::vector<SymbolicOperators::WickTerm>> wicks_M, wicks_N;
    
    /////////////
    // methods //
    /////////////
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
			    toClean(i) = ((2 * Constants::K_DISCRETIZATION) - std::abs(toClean(i) % (2 * Constants::K_DISCRETIZATION))) % (2 * Constants::K_DISCRETIZATION);
		    }
		    else {
			    toClean(i) = (toClean(i) % (2 * Constants::K_DISCRETIZATION)) % (2 * Constants::K_DISCRETIZATION);
		    }
	    }
    };
    inline Eigen::Vector2i computeMomentum(const SymbolicOperators::Momentum& momentum,
	    const std::vector<Eigen::Vector2i>& indizes, const std::vector<char>& momenta) const {
	    Eigen::Vector2i buffer = { 0,0 };
	    for (int i = 0; i < momenta.size(); ++i)
	    {
	    	int mom_idx = momentum.isUsed(momenta[i]);
	    	if (mom_idx < 0) continue;
	    	buffer += momentum.momentum_list[mom_idx].first * indizes[i];
	    }
	    if (momentum.add_Q) {
	    	buffer(0) += Constants::K_DISCRETIZATION;
	    	buffer(1) += Constants::K_DISCRETIZATION;
	    }
	    clean_factor_2pi(buffer);
	    return buffer;
	};
    double computeTerm(const SymbolicOperators::WickTerm& term, int l, int k) const;
    void fill_M_N();
	void fill_M_N_xp_basis();
};