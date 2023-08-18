#pragma once
#include "BaseDOS.hpp"

namespace Hubbard::DensityOfStates {
	struct SimpleCubic : public BaseDOS {
		static constexpr double LOWER_BORDER = -3;
		constexpr static int num_positions = 300U;
		// Needs to be a multiple of three of optimal accuracy
		// This is due to the discontinuity in the first derivative at +/- 1
		// It should also be even to reflect the symmetry of the DOS
		static int n_splits;
		static dos_precision b_minus_a_halved;

		static std::vector<std::pair<dos_precision, dos_precision>> split_limits;
		static std::array<dos_precision, num_positions> abscissa;
		static std::array<dos_precision, num_positions> weights;

		static constexpr int DIMENSION = 3;

		inline static size_t n_abscissa() noexcept {
			return num_positions * n_splits;
		};
		// Returns the offset in the gauss quadrature (a+b)/2 depending on the interval number i
		inline static dos_precision functionQuadratureOffset(int i) {
			return 0.5 * (split_limits[i].first + split_limits[i].second);
		}
		virtual void computeValues() override;

		template <class ResultType>
		class DOSIntegrator {
		private:
			ResultType result;
			ResultType buffer;

			template <bool byValue, class UnaryFunction>
			const ResultType& _internal_integrate(const UnaryFunction& F) {
				result *= dos_precision{};

				for (int i = 0; i < n_splits; ++i)
				{
					for (size_t j = 0U; j < num_positions; ++j)
					{
						if constexpr (byValue) {
							result += weights[j] * values[i * num_positions + j] * F(functionQuadratureOffset(i) + abscissa[j]);
						}
						else {
							F(functionQuadratureOffset(i) + abscissa[j], buffer);
							result += weights[j] * values[i * num_positions + j] * buffer;
						}
					}
				}

				return result;
			};
		public:
			// This function passes the result of F by reference,
			// i.e. F(gamma, result) and expects F to fill result accordingly.
			template <class UnaryFunction>
			inline const ResultType& integrate_by_reference(const UnaryFunction& F) {
				return _internal_integrate<false, UnaryFunction>(F);
			};

			// This function assumes that F returns its result by value.
			template <class UnaryFunction>
			inline const ResultType& integrate_by_value(const UnaryFunction& F) {
				return _internal_integrate<true, UnaryFunction>(F);
			};

			DOSIntegrator() = default;
			// If it is neccessary for the ResultType to be initaliazed,
			// e.g. give a vector a certain size. ResultType needs to have a copy constructor though
			DOSIntegrator(const ResultType& initialize_result)
				: result(initialize_result), buffer(initialize_result)
			{};
		};
	};
}