#pragma once
#include "../../../../FermionCommute/sources/WickTerm.hpp"
#include "../../Utility/InputFileReader.hpp"
#include "../SquareLattice/UsingBroyden.hpp"
#include <map>

// Both methods yield precisely the same data!
#define _PSEUDO_INVERSE

namespace Hubbard::Helper {
	constexpr double SQRT_SALT = 1e-5;
	constexpr double SALT = SQRT_SALT * SQRT_SALT;
	constexpr double ERROR_MARGIN = 1e-10;

	constexpr int OPERATION_NONE = 0;
	constexpr int OPERATION_INVERSE = 1;
	constexpr int OPERATION_SQRT = 2;
	constexpr int OPERATION_INVERSE_SQRT = 3;

	class ModeHelper {
	protected:
		std::vector<MatrixCL> expecs{};
		std::vector<complex_prec> sum_of_all{};

		std::vector<std::vector<SymbolicOperators::WickTerm>> wicks_M{}, wicks_N{};
		const std::map<std::string, int> wick_map = { {"n", 0}, {"g", 1}, {"f", 2}, {"\\eta", 3} };
		const std::map<std::string, int> wick_spin_offset = { {"\\uparrow", 0}, {"\\downarrow", 4}, {"\\sigma", 6} };

		std::unique_ptr<SquareLattice::UsingBroyden> model{};
		size_t TOTAL_BASIS{};
		/*
		* 0 - n
		* 1 - g_up
		* 2 - sc
		* 3 - eta
		* 4 - n_down
		* 5 - g_down
		* 6 - n_up + n_down
		* 7 - g_up + g_down
		*/

		int number_of_basis_terms{};
		int start_basis_at{};

		///////////////////////
		// Utility functions //
		///////////////////////
		// maps an index; [0, N_K) -> [-pi, pi)
		template <typename T>
		inline global_floating_type index_to_k_vector(const T index) const {
			return (((index * LONG_PI) / Constants::K_DISCRETIZATION) - LONG_PI);
		};
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

		/////////////
		// methods //
		/////////////
		void loadWick(const std::string& filename);
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

		complex_prec getExpectationValue(const SymbolicOperators::WickOperator& op, const Eigen::Vector2i& momentum_value) const;

		complex_prec getSumOfAll(const SymbolicOperators::WickOperator& op) const;

		complex_prec computeTerm(const SymbolicOperators::WickTerm& term, int l, int k) const;

		virtual void fillMatrices() = 0;

		/* Takes a positive semidefinite vector (the idea is that this contains eigenvalues) and applies an operation on it
		* 0: Correct for negative eigenvalues
		* 1: Compute the pseudoinverse
		* 2: Compute the square root
		* 3: Compute the pseudoinverse square root
		*/
		template<const int option>
		void applyMatrixOperation(Vector_L& evs) const {
			for (size_t i = 0; i < evs.size(); i++)
			{
				if (evs(i) < -SALT) {
					std::cerr << "M:   " << evs(i) << std::endl;
					throw std::invalid_argument("Matrix is not positive!  " + to_string(evs(i)));
				}
				if (evs(i) < SALT) {
#ifdef _PSEUDO_INVERSE
					evs(i) = 0;
#else
					switch (option) {
					default:
						break;
					case 1:
						evs(i) = 1 / SALT;
						break;

					case 2:
						evs(i) = SQRT_SALT;
						break;
					case 3:
						evs(i) = 1. / SQRT_SALT;
					}
#endif
				}
				else {
					switch (option) {
					default:
						break;
					case 1:
						evs(i) = 1. / evs(i);
						break;
					case 2:
						evs(i) = sqrt(evs(i));
						break;
					case 3:
						evs(i) = 1. / sqrt(evs(i));
					}
				}
			}
		};

	public:
		ModeHelper(Utility::InputFileReader& input);
		virtual ~ModeHelper() = default;

		SquareLattice::UsingBroyden& getModel() const {
			return *model;
		}

		virtual std::vector<Resolvent_L> computeCollectiveModes(std::vector<std::vector<global_floating_type>>& reciever) = 0;
	};
}