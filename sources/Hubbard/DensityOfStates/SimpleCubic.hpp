#pragma once
#include "BaseDOS.hpp"

namespace Hubbard::DensityOfStates {
	struct SimpleCubic : public BaseDOS {
		static std::vector<abscissa_t> upper_border_to_abscissa;
		static dos_precision b_minus_a_halved;

		static constexpr double LOWER_BORDER = -3;
		static constexpr int DIMENSION = 3;
		static constexpr int COORDINATION_NUMBER = 6;
		static constexpr int n_splits = 3;

		inline static size_t n_abscissa() noexcept {
			return abscissa.size();
		};
		static global_floating_type computeValue(const global_floating_type& gamma);
		virtual void computeValues() override;

		template <class T>
		using Integrator = DOSIntegrator<T, SimpleCubic>;
	};
}