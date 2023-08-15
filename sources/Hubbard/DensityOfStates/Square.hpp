#pragma once
#include <boost/multiprecision/float128.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include "BaseDOS.hpp"
#include <cmath>

namespace Hubbard::DensityOfStates {
	//typedef boost::multiprecision::float128 abscissa_t;
	typedef boost::multiprecision::cpp_bin_float_100 abscissa_t;
	struct Square : public BaseDOS {
		static std::vector<abscissa_t> abscissa;
		static std::vector<abscissa_t> upper_border_to_abscissa;
		static std::vector<dos_precision> weights;

		static dos_precision LOWER_BORDER;
		static dos_precision b_minus_a_halved;

		virtual void computeValues() override;
		inline static size_t n_abscissa() noexcept {
			return abscissa.size();
		};

		template <class ResultType>
		class DOSIntegrator {
		private:
			ResultType result;
			ResultType buffer;

			template <bool byValue, class UnaryFunction>
			const ResultType& _internal_integrate(const UnaryFunction& F) {
				if constexpr (byValue) {
					result = values[0] * static_cast<dos_precision>(weights[0]) * F(static_cast<dos_precision>(abscissa.front()));
				}
				else {
					F(static_cast<dos_precision>(abscissa.front()), buffer);
					result = values[0] * static_cast<dos_precision>(weights[0]) * buffer;
				}
				for (size_t i = 1U; i < values.size(); ++i)
				{
					if constexpr (byValue) {
						result += values[i] * static_cast<dos_precision>(weights[i]) * F(static_cast<dos_precision>(abscissa[i]));
					}
					else {
						F(static_cast<dos_precision>(abscissa[i]), buffer);
						result += values[i] * static_cast<dos_precision>(weights[i]) * buffer;
					}
				}
				result *= 2;
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